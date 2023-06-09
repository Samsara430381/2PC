#include "MyProtocol/MyProtocol.h"

using namespace std;
using namespace sci;

#define KKOT_LIMIT 8
#define SQRT_LOOKUP_SCALE 2

#define exp_TR
#define ADD_PRECISION 9

MyProtocol::MyProtocol(int party, IOPack *iopack, OTPack *otpack)
{
  this->party = party;
  this->iopack = iopack;
  this->otpack = otpack;
  this->aux = new AuxProtocols(party, iopack, otpack);
  this->xt = new XTProtocol(party, iopack, otpack);
  this->trunc = new Truncation(party, iopack, otpack);
  this->mult = new LinearOT(party, iopack, otpack);
}

MyProtocol::~MyProtocol()
{
  delete aux;
  delete xt;
  delete trunc;
  delete mult;
}

uint64_t Fix(double x, int32_t s_x)
{
  uint64_t y;
  y = (x * (1 << s_x));
  return y;
}

double urt(uint64_t x, int32_t s_x)
{
  double y;
  y = ((double)(uint64_t)(x) / (double)(uint64_t)(1 << s_x));
  return y;
}
// Protocol 1
void MyProtocol::secExp(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                        int32_t bw_y, int32_t s_x, int32_t s_y)
{
  int s_tmp = s_x + ADD_PRECISION;
  assert((bw_y - s_y + 2 * s_tmp) <= 64 && "bw_y - s_y + 2 * s_tmp (s_tmp = s_x + ADD_PRECISION) is not allowed to be greater than 64");
  double L = exp2(bw_x - s_x);
  uint64_t masky = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));
  uint64_t maskx = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  // bw_exp_x is selected according to the bit-width of the integer part of x, e.g. 4bit-->12   3bit-->6  2bit-->3  1bit-->2
  int32_t bw_exp_x;
  uint32_t choice = bw_x - s_x;
  switch (choice)
  {
  case 1:
    bw_exp_x = 2;
    break;
  case 2:
    bw_exp_x = 3;
    break;
  case 3:
    bw_exp_x = 6;
    break;
  default:
    bw_exp_x = 12;
  }
  int32_t bw_exp_mul = bw_exp_x + s_tmp;
  // 1.Compute MSB(x)
  uint8_t *msb_x = new uint8_t[dim];
  // uint64_t num_rounds = iopack->get_rounds();
  aux->MSB(x, msb_x, dim, bw_x);
  // num_rounds = iopack->get_rounds() - num_rounds;
  // cout << "MSB num_rounds:" << num_rounds << endl;
  // 2.Compute Wrap(x)
  uint8_t *wrap_x = new uint8_t[dim];
  aux->MSB_to_Wrap(x, msb_x, wrap_x, dim, bw_x);
  // 3.local Compute exp(x)
  double *x_exp_double1 = new double[dim];
  uint64_t *x_exp_mul1 = new uint64_t[dim];
  uint64_t *x_exp_mul2 = new uint64_t[dim];
  uint64_t *x_exp_mul3 = new uint64_t[dim];

  for (int32_t i = 0; i < dim; i++)
  {
    x_exp_double1[i] = exp(urt(x[i], s_x));                    // exp(Uint64_t2Float(x[i]));
    x_exp_mul1[i] = Fix(x_exp_double1[i], s_tmp);              // Float2Uint64(x_exp_double1[i]);
    x_exp_mul2[i] = Fix(x_exp_double1[i] / exp(L / 2), s_tmp); // Float2Uint64(x_exp_double1[i] / exp(8.0));
    x_exp_mul3[i] = Fix(x_exp_double1[i] / exp(L), s_tmp);     // Float2Uint64(x_exp_double1[i] / exp(16.0));
  }
  // 4.Compute secExp(x)
  //    1) Compute exp(x-wL)
  uint64_t *x_exp_add1 = new uint64_t[dim];
  uint64_t *x_exp_add2 = new uint64_t[dim];
  uint64_t *x_exp_add3 = new uint64_t[dim];

  MultMode mode = MultMode::Alice_has_A;
  //    get [exp(x)] form <exp(x)>
  mult->hadamard_cross_terms(dim, x_exp_mul1, x_exp_mul1, x_exp_add1, bw_exp_mul, bw_exp_mul, bw_y - s_y + 2 * s_tmp, mode);
  //    get [exp(x-L)] form <exp(x-L)>
  mult->hadamard_cross_terms(dim, x_exp_mul2, x_exp_mul2, x_exp_add2, bw_exp_mul, bw_exp_mul, bw_y - s_y + 2 * s_tmp, mode);
  //    get [exp(x-2L)] form <exp(x-2L)>
  mult->hadamard_cross_terms(dim, x_exp_mul3, x_exp_mul3, x_exp_add3, bw_exp_mul, bw_exp_mul, bw_y - s_y + 2 * s_tmp, mode);
  int32_t shift = 2 * s_tmp - s_y;
  //    secret select share
  //      1) Compute msb & wrap
  uint8_t *mw_x = new uint8_t[dim];
  aux->AND(msb_x, wrap_x, mw_x, dim);
  //      2) Compute msb ^ wrap
  uint8_t *mw1_x = new uint8_t[dim];
  for (int32_t i = 0; i < dim; i++)
  {
    mw1_x[i] = msb_x[i] ^ wrap_x[i];
    x_exp_add3[i] = x_exp_add3[i] - x_exp_add1[i];
    // x_exp_add3[i] = ((x_exp_add3[i] - x_exp_add1[i]) >> shift) & masky;
    // x_exp_add1[i] = (x_exp_add1[i] >> shift) & masky;
    // x_exp_add2[i] = (x_exp_add2[i] >> shift) & masky;
  }
  //    secret select share
  uint64_t *x_exp_add_sel = new uint64_t[dim];

  aux->multiplexer(mw_x, x_exp_add3, x_exp_add_sel, dim, bw_y + shift, bw_y + shift);
  uint64_t mask_out = ((bw_y + shift) == 64 ? -1 : ((1ULL << (bw_y + shift)) - 1));
  for (int32_t i = 0; i < dim; i++)
  {
    x_exp_add_sel[i] += x_exp_add1[i];
    x_exp_add_sel[i] &= mask_out;

    x_exp_add2[i] -= x_exp_add_sel[i];
    x_exp_add2[i] &= mask_out;
  }
  uint64_t *x_exp_add_sel2 = new uint64_t[dim];
  aux->multiplexer(mw1_x, x_exp_add2, x_exp_add_sel2, dim, bw_y + shift, bw_y + shift);
  for (int32_t i = 0; i < dim; i++)
  {
    y[i] = (x_exp_add_sel2[i] + x_exp_add_sel[i]) & mask_out;
  }
  trunc->truncate_and_reduce(dim, y, y, shift, bw_y - s_y + 2 * s_tmp);
  // for (int32_t i = 0; i < dim; i++)
  // {
  //   y[i] = (x_exp_add_sel2[i] + x_exp_add_sel[i]) & masky;
  // }
  delete[] msb_x;
  delete[] wrap_x;
  delete[] x_exp_double1;
  delete[] x_exp_mul1;
  delete[] x_exp_mul2;
  delete[] x_exp_mul3;
  delete[] x_exp_add1;
  delete[] x_exp_add2;
  delete[] x_exp_add3;
  delete[] mw_x;
  delete[] mw1_x;
  delete[] x_exp_add_sel;
  delete[] x_exp_add_sel2;
}

#ifdef exp_TR
// // Protocol 2 (Faithful  Truncate-and-Reduce)
// // Assumes x is always negative
void MyProtocol::secExp_N(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                          int32_t bw_y, int32_t s_x, int32_t s_y)
{
  int s_tmp = s_x + ADD_PRECISION;
  assert((bw_y - s_y + 2 * s_tmp) <= 64 && "bw_y - s_y + 2 * s_tmp (s_tmp = s_x + ADD_PRECISION) is not allowed to be greater than 64");
  double L = exp2(bw_x - s_x);
  uint64_t masky = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));
  uint64_t maskx = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  // bw_exp_x is selected according to the bit-width of the integer part of x, e.g. 4bit-->12   3bit-->6  2bit-->3  1bit-->2
  int32_t bw_exp_x;
  uint32_t choice = bw_x - s_x;
  switch (choice)
  {
  case 1:
    bw_exp_x = 2;
    break;
  case 2:
    bw_exp_x = 3;
    break;
  case 3:
    bw_exp_x = 6;
    break;
  default:
    bw_exp_x = 12;
  }
  int32_t bw_exp_mul = bw_exp_x + s_tmp;

  // 1.local Compute exp(x)
  double *x_exp_double1 = new double[dim];
  uint64_t *x_exp_mul = new uint64_t[dim * 2];

  for (int32_t i = 0; i < dim; i++)
  {
    x_exp_double1[i] = exp(urt(x[i], s_x));
    x_exp_mul[i] = Fix(x_exp_double1[i] / exp(L / 2), s_tmp);   // s1 = exp2(bw_x - s_x -1), exp(s1)=exp(8.0)
    x_exp_mul[i + dim] = Fix(x_exp_double1[i] / exp(L), s_tmp); // s1 = exp2(bw_x - s_x -1), exp(2.0 * s1)=exp(16.0)
  }

  // 2.Set msb(x)=1
  uint8_t *msb_x = new uint8_t[dim];
  if (party == ALICE)
  {
    for (int32_t i; i < dim; i++)
    {
      msb_x[i] = 0;
    }
  }
  else
  {
    for (int32_t i; i < dim; i++)
    {
      msb_x[i] = 1;
    }
  }
  // 3.Compute wrap(x)
  uint8_t *wrap_x = new uint8_t[dim];
  aux->MSB_to_Wrap(x, msb_x, wrap_x, dim, bw_x);

  // 4.Compute secExp(x)
  //    1) Compute exp(x-wL)
  uint64_t *x_exp_add = new uint64_t[dim * 2];
  uint64_t *x_exp_add_2 = x_exp_add + dim;
  // int32_t bw_exp_x = (bw_x - s_x - 1) * 3;
  MultMode mode = MultMode::Alice_has_A;
  mult->hadamard_cross_terms(dim * 2, x_exp_mul, x_exp_mul, x_exp_add, bw_exp_mul, bw_exp_mul, bw_y - s_y + 2 * s_tmp, mode);
  // cout << "bw_exp_x: " << bw_exp_x << " \tbw_exp_mul: " << bw_exp_mul << "  \tbw_y: " << bw_y << endl;
  int32_t shift = 2 * s_tmp - s_y;
  for (int32_t i; i < dim; i++)
  {
    x_exp_add[i + dim] = x_exp_add[i + dim] - x_exp_add[i];
    // x_exp_add[i + dim] = ((x_exp_add[i + dim] - x_exp_add[i]) >> shift) & masky;
    // x_exp_add[i] = (x_exp_add[i] >> shift) & masky;
  }

  aux->multiplexer(wrap_x, x_exp_add_2, x_exp_add_2, dim, bw_y - s_y + 2 * s_tmp, bw_y - s_y + 2 * s_tmp);
  // aux->multiplexer(wrap_x, x_exp_add_2, x_exp_add_2, dim, bw_y, bw_y);

  for (int32_t i; i < dim; i++)
  {
    y[i] = x_exp_add[i] + x_exp_add[i + dim];
    // y[i] = (x_exp_add[i] + x_exp_add[i + dim]) & masky;
  }

  trunc->truncate_and_reduce(dim, y, y, shift, bw_y - s_y + 2 * s_tmp);

  delete[] msb_x;
  delete[] wrap_x;
  delete[] x_exp_double1;
  delete[] x_exp_mul;
  delete[] x_exp_add;
}
#else
// // Protocol 2 (Local Truncate-and-Reduce)
// // Assumes x is always negative
void MyProtocol::secExp_N(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                          int32_t bw_y, int32_t s_x, int32_t s_y)
{
  assert((bw_y - s_y + 2 * s_tmp) <= 64 && "bw_y - s_y + 2 * s_tmp (s_tmp = s_x + ADD_PRECISION) is not allowed to be greater than 64");
  double L = exp2(bw_x - s_x);
  uint64_t masky = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));
  uint64_t maskx = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  // bw_exp_x is selected according to the bit-width of the integer part of x, e.g. 4bit-->12   3bit-->6  2bit-->3  1bit-->2
  int32_t bw_exp_x;
  uint32_t choice = bw_x - s_x;
  switch (choice)
  {
  case 1:
    bw_exp_x = 2;
    break;
  case 2:
    bw_exp_x = 3;
    break;
  case 3:
    bw_exp_x = 6;
    break;
  default:
    bw_exp_x = 12;
  }
  int32_t bw_exp_mul = bw_exp_x + s_tmp;

  // 1.local Compute exp(x)
  double *x_exp_double1 = new double[dim];
  uint64_t *x_exp_mul = new uint64_t[dim * 2];

  for (int32_t i = 0; i < dim; i++)
  {
    x_exp_double1[i] = exp(urt(x[i], s_x));
    x_exp_mul[i] = Fix(x_exp_double1[i] / exp(L / 2), s_tmp);   // s1 = exp2(bw_x - s_x -1), exp(s1)=exp(8.0)
    x_exp_mul[i + dim] = Fix(x_exp_double1[i] / exp(L), s_tmp); // s1 = exp2(bw_x - s_x -1), exp(2.0 * s1)=exp(16.0)
  }

  // 2.Set msb(x)=1
  uint8_t *msb_x = new uint8_t[dim];
  if (party == ALICE)
  {
    for (int32_t i; i < dim; i++)
    {
      msb_x[i] = 0;
    }
  }
  else
  {
    for (int32_t i; i < dim; i++)
    {
      msb_x[i] = 1;
    }
  }
  // 3.Compute wrap(x)
  uint8_t *wrap_x = new uint8_t[dim];
  aux->MSB_to_Wrap(x, msb_x, wrap_x, dim, bw_x);

  // 4.Compute secExp(x)
  //    1) Compute exp(x-wL)
  uint64_t *x_exp_add = new uint64_t[dim * 2];
  uint64_t *x_exp_add_2 = x_exp_add + dim;
  MultMode mode = MultMode::Alice_has_A;
  mult->hadamard_cross_terms(dim * 2, x_exp_mul, x_exp_mul, x_exp_add, bw_exp_mul, bw_exp_mul, bw_y - s_y + 2 * s_tmp, mode);
  int32_t shift = 2 * s_tmp - s_y;
  for (int32_t i; i < dim; i++)
  {
    // x_exp_add[i + dim] = x_exp_add[i + dim] - x_exp_add[i];
    x_exp_add[i + dim] = ((x_exp_add[i + dim] - x_exp_add[i]) >> shift) & masky;
    x_exp_add[i] = (x_exp_add[i] >> shift) & masky;
  }

  // aux->multiplexer(wrap_x, x_exp_add_2, x_exp_add_2, dim, bw_y - s_y + 2 * s_tmp, bw_y - s_y + 2 * s_tmp);
  aux->multiplexer(wrap_x, x_exp_add_2, x_exp_add_2, dim, bw_y, bw_y);

  for (int32_t i; i < dim; i++)
  {
    // y[i] = x_exp_add[i] + x_exp_add[i + dim];
    y[i] = (x_exp_add[i] + x_exp_add[i + dim]) & masky;
  }

  // trunc->truncate_and_reduce(dim, y, y, shift, bw_y - s_y + 2 * s_tmp);

  delete[] msb_x;
  delete[] wrap_x;
  delete[] x_exp_double1;
  delete[] x_exp_mul;
  delete[] x_exp_add;
}
#endif

uint64_t lookup_neg_exp(uint64_t val_in, int32_t s_in, int32_t s_out)
{
  if (s_in < 0)
  {
    s_in *= -1;
    val_in *= (1ULL << (s_in));
    s_in = 0;
  }
  uint64_t res_val =
      exp(-1.0 * (val_in / double(1ULL << s_in))) * (1ULL << s_out);
  return res_val;
}

void MyProtocol::lookup_table_exp(int32_t dim, uint64_t *x, uint64_t *y,
                                  int32_t bw_x, int32_t bw_y, int32_t s_x,
                                  int32_t s_y)
{
  assert(bw_y >= (s_y + 2));

  int LUT_size = KKOT_LIMIT;

  uint64_t s_tmpmask = (bw_x == 64 ? -1 : (1ULL << bw_x) - 1);
  uint64_t LUT_out_mask = ((s_y + 2) == 64 ? -1 : (1ULL << (s_y + 2)) - 1);

  uint64_t *tmp_1 = new uint64_t[dim];
  for (int i = 0; i < dim; i++)
  {
    tmp_1[i] = (-1 * x[i]) & s_tmpmask;
  }
  int digit_size = LUT_size;
  int num_digits = ceil(double(bw_x) / digit_size);
  int last_digit_size = bw_x - (num_digits - 1) * digit_size;
  uint64_t *x_digits = new uint64_t[num_digits * dim];

  aux->digit_decomposition_sci(dim, tmp_1, x_digits, bw_x, digit_size);

  uint64_t digit_mask = (digit_size == 64 ? -1 : (1ULL << digit_size) - 1);
  uint64_t last_digit_mask =
      (last_digit_size == 64 ? -1 : (1ULL << last_digit_size) - 1);
  int N = (1ULL << digit_size);
  int last_N = (1ULL << last_digit_size);
  int N_digits = (digit_size == last_digit_size ? num_digits : num_digits - 1);
  uint64_t *digits_exp = new uint64_t[num_digits * dim];
  if (party == ALICE)
  {
    uint64_t **spec;
    spec = new uint64_t *[num_digits * dim];
    PRG128 prg;
    prg.random_data(digits_exp, num_digits * dim * sizeof(uint64_t));
    for (int i = 0; i < N_digits * dim; i++)
    {
      int digit_idx = i / dim;
      spec[i] = new uint64_t[N];
      digits_exp[i] &= LUT_out_mask;
      for (int j = 0; j < N; j++)
      {
        int idx = (x_digits[i] + j) & digit_mask;
        spec[i][j] = (lookup_neg_exp(idx, s_x - digit_size * digit_idx, s_y) -
                      digits_exp[i]) &
                     LUT_out_mask;
      }
    }
    aux->lookup_table<uint64_t>(spec, nullptr, nullptr, N_digits * dim,
                                digit_size, s_y + 2);

    if (digit_size != last_digit_size)
    {
      int offset = N_digits * dim;
      int digit_idx = N_digits;
      for (int i = offset; i < num_digits * dim; i++)
      {
        spec[i] = new uint64_t[last_N];
        digits_exp[i] &= LUT_out_mask;
        for (int j = 0; j < last_N; j++)
        {
          int idx = (x_digits[i] + j) & last_digit_mask;
          spec[i][j] = (lookup_neg_exp(idx, s_x - digit_size * digit_idx, s_y) -
                        digits_exp[i]) &
                       LUT_out_mask;
        }
      }
      aux->lookup_table<uint64_t>(spec + offset, nullptr, nullptr, dim,
                                  last_digit_size, s_y + 2);
    }

    for (int i = 0; i < num_digits * dim; i++)
      delete[] spec[i];
    delete[] spec;
  }
  else
  {
    aux->lookup_table<uint64_t>(nullptr, x_digits, digits_exp, N_digits * dim,
                                digit_size, s_y + 2);
    if (digit_size != last_digit_size)
    {
      int offset = N_digits * dim;
      aux->lookup_table<uint64_t>(nullptr, x_digits + offset,
                                  digits_exp + offset, dim, last_digit_size,
                                  s_y + 2);
    }
    for (int i = 0; i < num_digits * dim; i++)
    {
      digits_exp[i] &= LUT_out_mask;
    }
  }

  uint8_t *zero_shares = new uint8_t[dim];
  for (int i = 0; i < dim; i++)
  {
    zero_shares[i] = 0;
  }
  for (int i = 1; i < num_digits; i *= 2)
  {
    for (int j = 0; j < num_digits and j + i < num_digits; j += 2 * i)
    {
      mult->hadamard_product(dim, digits_exp + j * dim,
                             digits_exp + (j + i) * dim, digits_exp + j * dim,
                             s_y + 2, s_y + 2, 2 * s_y + 2, false, false,
                             MultMode::None, zero_shares, zero_shares);
      trunc->truncate_and_reduce(dim, digits_exp + j * dim,
                                 digits_exp + j * dim, s_y, 2 * s_y + 2);
    }
  }
  xt->z_extend(dim, digits_exp, y, s_y + 2, bw_y, zero_shares);

  delete[] x_digits;
  delete[] tmp_1;
  delete[] digits_exp;
  delete[] zero_shares;
}

#ifndef secExpn
#define secExpn
void MyProtocol::secExp_N2(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                           int32_t bw_y, int32_t s_x, int32_t s_y)
{
  int s_tmp = s_x + ADD_PRECISION;
  assert((bw_y - s_y + 2 * s_tmp) <= 64 && "bw_y - s_y + 2 * s_tmp (s_tmp = s_x + ADD_PRECISION) is not allowed to be greater than 64");
  uint64_t masky = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));
  uint64_t maskx = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  if (bw_x - s_x > 4)
  {
    // x = u || v
    uint64_t mask_u = ((bw_x - s_x) == 64 ? -1 : ((1ULL << (bw_x - s_x)) - 1));
    uint64_t mask_v = (s_x == 64 ? -1 : ((1ULL << s_x) - 1));
    uint64_t *u = new uint64_t[dim];
    uint64_t *v = new uint64_t[dim];
    uint64_t *exp_u = new uint64_t[dim];

    // trunc->truncate_and_reduce(dim, x, u, s_x, bw_x);
    uint64_t *inA_lower = new uint64_t[dim];
    uint8_t *wrap_v = new uint8_t[dim];
    for (int i = 0; i < dim; i++)
    {
      inA_lower[i] = x[i] & mask_v;
    }
    this->aux->wrap_computation(inA_lower, wrap_v, dim, s_x);
    uint64_t *arith_wrap = new uint64_t[dim];
    this->aux->B2A(wrap_v, arith_wrap, dim, (bw_x - s_x));
    for (int i = 0; i < dim; i++)
    {
      u[i] = ((x[i] >> s_x) + arith_wrap[i]) & mask_u;
    }
    delete[] inA_lower;
    delete[] arith_wrap;

    for (int32_t i = 0; i < dim; i++)
    {
      // u[i] = (x[i] >> s_x) & mask_u;
      v[i] = x[i] & mask_v;
    }
    // Compute exp(u)
    lookup_table_exp(dim, u, exp_u, bw_x - s_x, s_y + 2, 0, s_y);
    // Compute exp(v)
    // 1.local Compute exp(x)
    double *v_exp_double = new double[dim];
    uint64_t *v_exp_mul = new uint64_t[2 * dim];
    double L = 1.0;
    for (int32_t i = 0; i < dim; i++)
    {
      v_exp_double[i] = exp(urt(v[i], s_x));
      v_exp_mul[i] = Fix(v_exp_double[i], s_tmp);
      v_exp_mul[i + dim] = Fix(v_exp_double[i] / exp(L / 2), s_tmp);
    }
    int32_t bw_exp_mul = 2 + s_tmp;
    uint64_t *v_exp_add = new uint64_t[2 * dim];
    MultMode mode = MultMode::Alice_has_A;
    mult->hadamard_cross_terms(dim * 2, v_exp_mul, v_exp_mul, v_exp_add, bw_exp_mul, bw_exp_mul, 2 * bw_exp_mul, mode);
    for (int32_t i; i < dim; i++)
    {
      v_exp_add[i + dim] = v_exp_add[i + dim] - v_exp_add[i];
    }
    uint64_t *v_exp_add_2 = v_exp_add + dim;
    aux->multiplexer(wrap_v, v_exp_add_2, v_exp_add_2, dim, 2 * bw_exp_mul, 2 * bw_exp_mul);
    for (int32_t i; i < dim; i++)
    {
      y[i] = v_exp_add[i] + v_exp_add[i + dim];
    }

    uint8_t *zero_shares = new uint8_t[dim];
    for (int i = 0; i < dim; i++)
    {
      zero_shares[i] = 0;
    }
    mult->hadamard_product(dim, exp_u, y, y,
                           s_y + 2, 2 * bw_exp_mul, s_y + 2 + 2 * s_tmp, false, false,
                           MultMode::None, zero_shares, zero_shares);
    trunc->truncate_and_reduce(dim, y, y, 2 * s_tmp, s_y + 2 + 2 * s_tmp);
    xt->z_extend(dim, y, y, s_y + 2, bw_y, zero_shares);
    delete[] u;
    delete[] v;
    delete[] exp_u;
    delete[] v_exp_double;
    delete[] v_exp_mul;
    delete[] v_exp_add;
    delete[] wrap_v;
    delete[] zero_shares;
  }
  else
  {
    secExp_N(dim, x, y, bw_x, bw_y, s_x, s_y);
  }
}
#endif

// A0 \in (1/4, 1)
uint64_t lookup_A0(uint64_t index, int m)
{
  uint64_t k = 1ULL << m;
  double p = 1 + (double(index) / double(k));
  double A1 = 1.0 / (p * (p + 1.0 / double(k)));
  int32_t scale = m + 3;
  uint64_t mask = (1ULL << scale) - 1;
  uint64_t val = uint64_t(A1 * (1ULL << scale)) & mask;
  return val;
}

// A1 \in (1/2, 1)
uint64_t lookup_A1(uint64_t index, int m)
{
  uint64_t k = 1ULL << m;
  double p = 1 + (double(index) / double(k));
  double z = (p * (p + (1.0 / double(k))));
  double A1 = ((1.0 / double(k * 2)) + sqrt(z)) / z;
  int32_t scale = 2 * m + 2;
  uint64_t mask = (1ULL << scale) - 1;
  uint64_t val = uint64_t(A1 * (1ULL << scale)) & mask;
  return val;
}

void MyProtocol::reciprocal_approximation(int32_t dim, int32_t m,
                                          uint64_t *dn, uint64_t *out,
                                          int32_t bw_dn, int32_t bw_out,
                                          int32_t s_dn, int32_t s_out)
{
  assert(bw_out == m + s_dn + 4);
  assert(s_out == m + s_dn + 4);

  uint64_t s_dn_mask = (1ULL << s_dn) - 1;
  uint64_t m_mask = (1ULL << m) - 1;
  uint64_t s_min_m_mask = (1ULL << (s_dn - m)) - 1;

  uint64_t *tmp_1 = new uint64_t[dim];
  uint64_t *tmp_2 = new uint64_t[dim];

  for (int i = 0; i < dim; i++)
  {
    tmp_1[i] = dn[i] & s_dn_mask;
  }
  trunc->truncate_and_reduce(dim, tmp_1, tmp_2, s_dn - m, s_dn);

  int M = (1ULL << m);
  uint64_t c0_mask = (1ULL << (m + 4)) - 1;
  uint64_t c1_mask = (1ULL << (2 * m + 3)) - 1;
  uint64_t *c0 = new uint64_t[dim];
  uint64_t *c1 = new uint64_t[dim];
  if (party == ALICE)
  {
    uint64_t **spec;
    spec = new uint64_t *[dim];
    PRG128 prg;
    prg.random_data(c0, dim * sizeof(uint64_t));
    prg.random_data(c1, dim * sizeof(uint64_t));
    for (int i = 0; i < dim; i++)
    {
      spec[i] = new uint64_t[M];
      c0[i] &= c0_mask;
      c1[i] &= c1_mask;
      for (int j = 0; j < M; j++)
      {
        int idx = (tmp_2[i] + j) & m_mask;
        spec[i][j] = (lookup_A0(idx, m) - c0[i]) & c0_mask;
        spec[i][j] <<= (2 * m + 3);
        spec[i][j] |= (lookup_A1(idx, m) - c1[i]) & c1_mask;
      }
    }
    aux->lookup_table<uint64_t>(spec, nullptr, nullptr, dim, m, 3 * m + 7);

    for (int i = 0; i < dim; i++)
      delete[] spec[i];
    delete[] spec;
  }
  else
  {
    aux->lookup_table<uint64_t>(nullptr, tmp_2, c1, dim, m, 3 * m + 7);

    for (int i = 0; i < dim; i++)
    {
      c0[i] = (c1[i] >> (2 * m + 3)) & c0_mask;
      c1[i] = c1[i] & c1_mask;
    }
  }
  for (int i = 0; i < dim; i++)
  {
    tmp_1[i] = dn[i] & s_min_m_mask;
  }
  uint8_t *zero_shares = new uint8_t[dim];
  for (int i = 0; i < dim; i++)
  {
    zero_shares[i] = 0;
  }

  // Unsigned mult
  mult->hadamard_product(dim, c0, tmp_1, tmp_2, m + 4, s_dn - m, s_dn + 4,
                         false, false, MultMode::None, zero_shares, nullptr);

  xt->z_extend(dim, tmp_2, tmp_1, s_dn + 4, s_dn + m + 4, zero_shares);

  uint64_t out_mask = (1ULL << (s_dn + m + 4)) - 1;
  uint64_t scale_up = (1ULL << (s_dn - m + 1));
  for (int i = 0; i < dim; i++)
  {
    out[i] = ((c1[i] * scale_up) - tmp_1[i]) & out_mask;
  }

  delete[] tmp_1;
  delete[] tmp_2;
  delete[] c0;
  delete[] c1;
  delete[] zero_shares;
}

void MyProtocol::div(int32_t dim, uint64_t *nm, uint64_t *dn, uint64_t *out,
                     int32_t bw_nm, int32_t bw_dn, int32_t bw_out,
                     int32_t s_nm, int32_t s_dn, int32_t s_out,
                     bool signed_nm, bool compute_msnzb)
{
  assert(s_out <= s_dn);

  // out_precision = iters * (2*m + 2)
  int32_t m, iters;
  m = (s_out <= 18 ? ceil((s_out - 2) / 2.0)
                   : ceil((ceil(s_out / 2.0) - 2) / 2.0));
  iters = (s_out <= 18 ? 0 : 1);

  int32_t s_tmp_dn;
  int32_t bw_adjust;
  int32_t s_adjust;
  uint64_t *tmp_dn;
  uint64_t *adjust;
  if (compute_msnzb)
  {
    s_tmp_dn = bw_dn - 1;
    bw_adjust = bw_dn + 1;
    s_adjust = bw_dn - 1 - s_dn;
    uint64_t mask_adjust = (bw_adjust == 64 ? -1 : ((1ULL << bw_adjust) - 1));
    // MSB is always 0, thus, not including it
    uint8_t *msnzb_vector_bool = new uint8_t[dim * bw_dn];
    uint64_t *msnzb_vector = new uint64_t[dim * bw_dn];
    aux->msnzb_one_hot(dn, msnzb_vector_bool, bw_dn, dim);
    aux->B2A(msnzb_vector_bool, msnzb_vector, dim * bw_dn, bw_adjust);
    // adjust: bw = bw_dn, scale = bw_dn - 1 - s_dn
    adjust = new uint64_t[dim];
    for (int i = 0; i < dim; i++)
    {
      adjust[i] = 0;
      for (int j = 0; j < bw_dn; j++)
      {
        adjust[i] += (1ULL << (bw_dn - 1 - j)) * msnzb_vector[i * bw_dn + j];
      }
      adjust[i] &= mask_adjust;
    }
    // tmp_dn: bw = bw_dn, scale = bw_dn - 1
    tmp_dn = new uint64_t[dim];
    mult->hadamard_product(dim, dn, adjust, tmp_dn, bw_dn, bw_dn + 1, bw_dn + 1,
                           false, false, MultMode::None);

    delete[] msnzb_vector_bool;
    delete[] msnzb_vector;
  }
  else
  {
    // tmp_dn: bw = s_dn + 1, scale = s_dn
    s_tmp_dn = s_dn;
    tmp_dn = dn;
  }

  uint64_t *tmp_1 = new uint64_t[dim];
  uint64_t *tmp_2 = new uint64_t[dim];
  // tmp_1: bw = s_tmp_dn + m + 4, scale = s_tmp_dn + m + 3
  reciprocal_approximation(dim, m, tmp_dn, tmp_1, bw_dn, s_tmp_dn + m + 4,
                           s_tmp_dn, s_tmp_dn + m + 4);

  uint64_t *w0 = new uint64_t[dim];
  // w0: bw = s_out + 1, scale = s_out
  trunc->truncate_and_reduce(dim, tmp_1, w0, s_tmp_dn + m + 3 - s_out,
                             s_tmp_dn + m + 4);

  uint8_t *msb_nm = new uint8_t[dim];
  aux->MSB(nm, msb_nm, dim, bw_nm);
  uint8_t *zero_shares = new uint8_t[dim];
  for (int i = 0; i < dim; i++)
  {
    zero_shares[i] = 0;
  }

  // a0: bw = bw_out, scale = s_out
  uint64_t *a0 = new uint64_t[dim];
  // Mixed mult with w0 unsigned
  mult->hadamard_product(dim, nm, w0, tmp_1, bw_nm, s_out + 1, s_out + bw_nm,
                         signed_nm, false, MultMode::None, msb_nm, zero_shares);
  trunc->truncate_and_reduce(dim, tmp_1, tmp_2, s_nm, s_out + bw_nm);
  if ((bw_nm - s_nm) >= (bw_out - s_out))
  {
    aux->reduce(dim, tmp_2, a0, bw_nm - s_nm + s_out, bw_out);
  }
  else
  {
    if (signed_nm)
    {
      xt->s_extend(dim, tmp_2, a0, s_out + bw_nm - s_nm, bw_out, msb_nm);
    }
    else
    {
      xt->z_extend(dim, tmp_2, a0, s_out + bw_nm - s_nm, bw_out, nullptr);
    }
  }

  if (compute_msnzb)
  {
    int32_t bw_tmp1 =
        (bw_out + s_adjust < bw_adjust ? bw_adjust : bw_out + s_adjust);
    // tmp_1: bw = bw_tmp1, scale = s_out + s_adjust
    mult->hadamard_product(dim, a0, adjust, tmp_1, bw_out, bw_adjust, bw_tmp1,
                           signed_nm, false, MultMode::None,
                           (signed_nm ? msb_nm : nullptr), zero_shares);
    // a0: bw = bw_out, scale = s_out
    trunc->truncate_and_reduce(dim, tmp_1, a0, s_adjust, bw_out + s_adjust);
  }

  // tmp_1: bw = s_tmp_dn + 2, scale = s_tmp_dn
  uint64_t s_plus_2_mask = (1ULL << (s_tmp_dn + 2)) - 1;
  for (int i = 0; i < dim; i++)
  {
    tmp_1[i] = tmp_dn[i] & s_plus_2_mask;
  }

  if (iters > 0)
  {
    // d0: bw = s_out + 2, scale = s_out
    uint64_t *d0 = new uint64_t[dim];
    mult->hadamard_product(dim, w0, tmp_1, tmp_2, s_out + 1, s_tmp_dn + 2,
                           s_out + s_tmp_dn + 2, false, false, MultMode::None,
                           zero_shares, zero_shares);
    trunc->truncate_and_reduce(dim, tmp_2, d0, s_tmp_dn, s_out + s_tmp_dn + 2);

    // e0: bw = s_out + 2, scale = s_out
    // e0 = 1 - d0
    uint64_t *e0 = new uint64_t[dim];
    for (int i = 0; i < dim; i++)
    {
      e0[i] = (party == ALICE ? (1ULL << (s_out)) : 0) - d0[i];
    }

    uint64_t e_mask = (1ULL << (s_out + 2)) - 1;
    uint64_t *a_curr = new uint64_t[dim];
    uint64_t *e_curr = new uint64_t[dim];
    uint64_t *a_prev = a0;
    uint64_t *e_prev = e0;
    for (int i = 0; i < iters - 1; i++)
    {
      // tmp_1: bw = 2*s_out+2, scale: 2*s_out
      mult->hadamard_product(dim, e_prev, e_prev, tmp_1, s_out + 2, s_out + 2,
                             2 * s_out + 2, true, true, MultMode::None,
                             zero_shares, zero_shares);
      // e_curr: bw = s_out + 2, scale: s_out
      trunc->truncate_and_reduce(dim, tmp_1, e_curr, s_out, 2 * s_out + 2);
      // e_prev = 1 + e_prev
      for (int j = 0; j < dim; j++)
      {
        e_prev[j] =
            ((party == ALICE ? (1ULL << (s_out)) : 0) + e_prev[j]) & e_mask;
      }
      // tmp_1: bw = bw_out + s_out, scale: 2*s_out
      mult->hadamard_product(dim, a_prev, e_prev, tmp_1, bw_out, s_out + 2,
                             bw_out + s_out, signed_nm, true, MultMode::None,
                             (signed_nm ? msb_nm : nullptr), zero_shares);
      // a_curr: bw = bw_out, scale: s_out
      trunc->truncate_and_reduce(dim, tmp_1, a_curr, s_out, bw_out + s_out);
      memcpy(a_prev, a_curr, dim * sizeof(uint64_t));
      memcpy(e_prev, e_curr, dim * sizeof(uint64_t));
    }
    // e_prev = 1 + e_prev
    for (int j = 0; j < dim; j++)
    {
      e_prev[j] =
          ((party == ALICE ? (1ULL << (s_out)) : 0) + e_prev[j]) & e_mask;
    }
    // out: bw = bw_out, scale: s_out
    // Mixed mult with e_prev unsigned
    mult->hadamard_product(dim, a_prev, e_prev, tmp_1, bw_out, s_out + 2,
                           bw_out + s_out, signed_nm, false, MultMode::None,
                           (signed_nm ? msb_nm : nullptr), zero_shares);
    trunc->truncate_and_reduce(dim, tmp_1, out, s_out, bw_out + s_out);

    delete[] d0;
    delete[] e0;
    delete[] a_curr;
    delete[] e_curr;
  }
  else
  {
    memcpy(out, a0, dim * sizeof(uint64_t));
  }

  delete[] tmp_1;
  delete[] tmp_2;
  delete[] w0;
  delete[] a0;
  delete[] msb_nm;
  delete[] zero_shares;
  if (compute_msnzb)
  {
    delete[] tmp_dn;
    delete[] adjust;
  }
}

// #undef secExpn
// Protocol 3
void MyProtocol::secSigmoid(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                            int32_t bw_y, int32_t s_x, int32_t s_y)
{
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));
  uint64_t mask_exp_out = ((s_y + 2) == 64 ? -1 : ((1ULL << (s_y + 2)) - 1));
  uint8_t *zero_shares = new uint8_t[dim];

  for (int i = 0; i < dim; i++)
  {
    zero_shares[i] = 0;
  }
  uint8_t *msb_x = new uint8_t[dim];
  aux->MSB(x, msb_x, dim, bw_x);
  uint64_t *tmp_1 = new uint64_t[dim];
  uint64_t *tmp_2 = new uint64_t[dim];
  for (int i = 0; i < dim; i++)
  {
    tmp_1[i] = (-1 * x[i] - (party == ALICE ? 1 : 0)) & mask_x;
    tmp_2[i] = (2 * x[i] + (party == ALICE ? 1 : 0)) & mask_x;
    // tmp_1[i] = (-1 * x[i]) & mask_x;
    // tmp_2[i] = (2 * x[i]) & mask_x;
  }
  // neg_x = msb_x*(2x+1)-x-1 , neg_x is always negative
  uint64_t *neg_x = new uint64_t[dim];
  aux->multiplexer(msb_x, tmp_2, neg_x, dim, bw_x, bw_x);
  for (int i = 0; i < dim; i++)
  {
    neg_x[i] = (neg_x[i] + tmp_1[i]) & mask_x;
  }
  uint64_t *exp_neg_x = new uint64_t[dim];

#ifndef secExpn
  lookup_table_exp(dim, neg_x, exp_neg_x, bw_x, s_y + 2, s_x, s_y);
#else
  // assert((bw_x - s_x) <= 4 && "secExp_N can only support the integer part of up to 4bits");
  secExp_N2(dim, neg_x, exp_neg_x, bw_x, s_y + 2, s_x, s_y);
#endif

  for (int i = 0; i < dim; i++)
  {
    tmp_1[i] = ((party == ALICE ? (1ULL << s_y) : 0) + exp_neg_x[i]) & mask_exp_out;
    tmp_2[i] = (party == ALICE ? 1 : 0);
  }
  uint64_t *sig_neg_x = new uint64_t[dim];
  // sig_neg_x = 1/(1 + exp_neg_x)
  div(dim, tmp_2, tmp_1, sig_neg_x, 2, s_y + 2, s_y + 2, 0, s_y, s_y, true, false);
  // tmp_1 = (1 - 2 * sig_neg_x)
  for (int i = 0; i < dim; i++)
  {
    tmp_1[i] = ((party == ALICE ? 1ULL << s_y : 0) - 2 * sig_neg_x[i]) & mask_exp_out;
  }
  aux->multiplexer(msb_x, tmp_1, tmp_2, dim, s_y + 2, s_y + 2);
  for (int i = 0; i < dim; i++)
  {
    tmp_2[i] = (sig_neg_x[i] + tmp_2[i]) & mask_exp_out;
  }
  if (bw_y <= (s_y + 2))
  {
    for (int i = 0; i < dim; i++)
    {
      y[i] = tmp_2[i] & mask_y;
    }
  }
  else
  {
    xt->z_extend(dim, tmp_2, y, s_y + 2, bw_y, zero_shares);
  }

  delete[] zero_shares;
  delete[] msb_x;
  delete[] tmp_1;
  delete[] tmp_2;
  delete[] neg_x;
  delete[] exp_neg_x;
  delete[] sig_neg_x;
}

void MyProtocol::secTanh(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                         int32_t bw_y, int32_t s_x, int32_t s_y)
{

  // Compute Sigmoid(2x)
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));
  uint64_t mask_exp_out = ((s_y + 2) == 64 ? -1 : ((1ULL << (s_y + 2)) - 1));

  uint8_t *msb_x = new uint8_t[dim];
  aux->MSB(x, msb_x, dim, bw_x);
  // neg_x = -x-1 + msb_x * (2x + 1) (neg_x is always negative)
  uint64_t *tmp_1 = new uint64_t[dim];
  uint64_t *tmp_2 = new uint64_t[dim];
  for (int i = 0; i < dim; i++)
  {
    tmp_1[i] = (-1 * x[i] - (party == ALICE ? 1 : 0)) & mask_x;
    tmp_2[i] = (2 * x[i] + (party == ALICE ? 1 : 0)) & mask_x;
    // tmp_1[i] = (-1 * x[i]) & mask_x;
    // tmp_2[i] = (2 * x[i]) & mask_x;
  }
  uint64_t *neg_x = new uint64_t[dim];
  aux->multiplexer(msb_x, tmp_2, neg_x, dim, bw_x, bw_x);

  for (int i = 0; i < dim; i++)
  {
    neg_x[i] = (neg_x[i] + tmp_1[i]) & mask_x;
  }
  uint64_t *exp_neg_2x = new uint64_t[dim];
#ifndef secExpn
  lookup_table_exp(dim, neg_x, exp_neg_2x, bw_x, s_y + 2, s_x - 1, s_y);
#else
  // assert((bw_x - s_x) <= 3 && "secExp_N can only support the integer part of up to 4bits");
  secExp_N2(dim, neg_x, exp_neg_2x, bw_x, s_y + 2, s_x - 1, s_y);
#endif
  // tmp_1 = 1 + exp_neg_2x
  for (int i = 0; i < dim; i++)
  {
    tmp_1[i] = ((party == ALICE ? (1ULL << s_y) : 0) + exp_neg_2x[i]) & mask_exp_out;
    tmp_2[i] = (party == ALICE ? 1 : 0);
  }
  uint64_t *sig_neg_2x = new uint64_t[dim];
  // sig_neg_2x = 1/(1 + exp_neg_2x)
  div(dim, tmp_2, tmp_1, sig_neg_2x, 2, s_y + 2, s_y + 2, 0, s_y, s_y, true, false);
  for (int i = 0; i < dim; i++)
  {
    tmp_1[i] = ((party == ALICE ? 1ULL << s_y : 0) - 2 * sig_neg_2x[i]) & mask_exp_out;
  }
  // tmp_2 = msb_x * (1 - 2 * sig_neg_2x)
  aux->multiplexer(msb_x, tmp_1, tmp_2, dim, s_y + 2, s_y + 2);
  // tmp_1 = sig_neg2_x + msb_x * (1 - 2 * sig_neg_2x)
  for (int i = 0; i < dim; i++)
  {
    tmp_2[i] = (sig_neg_2x[i] + tmp_2[i]) & mask_exp_out;
    // Compute Tan_x = 2 * tmp_2 -1
    tmp_2[i] = (2 * tmp_2[i] - (party == ALICE ? 1ULL << s_y : 0)) & mask_exp_out;
  }

  if (bw_y <= (s_y + 2))
  {
    for (int i = 0; i < dim; i++)
    {
      y[i] = tmp_2[i] & mask_y;
    }
  }
  else
  {
    xt->s_extend(dim, tmp_2, y, s_y + 2, bw_y, msb_x);
  }

  delete[] msb_x;
  delete[] tmp_1;
  delete[] tmp_2;
  delete[] neg_x;
  delete[] exp_neg_2x;
}

//  0.5 < x =< 1
void MyProtocol::secNewtonRaphson(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x, int32_t bw_y, int32_t s_x, int32_t s_y)
{
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t *w0 = new uint64_t[dim];
  uint64_t *tmp = new uint64_t[dim]; // tmp = x * w0
  uint64_t *e0 = new uint64_t[dim];

  uint8_t *zero_shares = new uint8_t[dim];
  for (int i = 0; i < dim; i++)
  {
    zero_shares[i] = 0;
  }
  for (int i = 0; i < dim; i++)
  {
    w0[i] = ((party == ALICE ? Fix(2.9192, s_x) : 0)) - 2 * x[i];
  }
  mult->hadamard_product(dim, x, w0, tmp, bw_x, bw_x, 2 * s_x + 2, false, false, MultMode::None, zero_shares, zero_shares);
  for (int i = 0; i < dim; i++)
  {
    tmp[i] = (tmp[i] >> s_x) & mask_x;
    e0[i] = (Fix(1.0, s_x) - tmp[i]) & mask_x;
  }
  mult->hadamard_product(dim, w0, e0, y, bw_x, bw_x, 2 * s_x + 2, false, false, MultMode::None, zero_shares, zero_shares);
  for (int i = 0; i < dim; i++)
  {
    y[i] = (y[i] >> s_x) & mask_x;
  }
  delete[] w0;
  delete[] tmp;
  delete[] e0;
  delete[] zero_shares;
}