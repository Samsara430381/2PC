#ifndef MY_PROTOCOL_H__
#define MY_PROTOCOL_H__

#include "BuildingBlocks/aux-protocols.h"
#include "BuildingBlocks/truncation.h"
#include "BuildingBlocks/value-extension.h"
#include "LinearOT/linear-ot.h"

class MyProtocol
{
public:
    int party;
    sci::IOPack *iopack;
    sci::OTPack *otpack;
    AuxProtocols *aux;
    XTProtocol *xt;
    Truncation *trunc;
    LinearOT *mult;

    MyProtocol(int party, sci::IOPack *iopack, sci::OTPack *otpack);

    ~MyProtocol();

    void secExp(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                int32_t bw_y, int32_t s_x, int32_t s_y); // int32_t s_y = s_x

    // Assumes x is always negative
    void lookup_table_exp(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                          int32_t bw_y, int32_t s_x, int32_t s_y);

    // Assumes x is always negative
    void secExp_N(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                  int32_t bw_y, int32_t s_x, int32_t s_y); // int32_t s_y = s_x
    void secExp_N2(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                   int32_t bw_y, int32_t s_x, int32_t s_y); // int32_t s_y = s_x
    // Assumes x is always positive
    void secExp_P(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                  int32_t bw_y, int32_t s_x); // int32_t s_y = s_x

    // Current implementation assumes that dn is always of the form 1.y1y2y3..yn
    void reciprocal_approximation(int32_t dim, int32_t m, uint64_t *dn,
                                  uint64_t *out, int32_t bw_dn, int32_t bw_out,
                                  int32_t s_dn, int32_t s_out);

    // If compute_msnzb = false, dn = 1.y1y2y3....
    // Else if compute_msnzb = true, dn is always positive
    void div(int32_t dim,
             // numerator
             uint64_t *nm,
             // denominator
             uint64_t *dn,
             // output
             uint64_t *out,
             // bitwidths
             int32_t bw_nm, int32_t bw_dn, int32_t bw_out,
             // scales
             int32_t s_nm, int32_t s_dn, int32_t s_out, bool signed_nm = true,
             bool compute_msnzb = false);

    void secSigmoid(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                    int32_t bw_y, int32_t s_x, int32_t s_y); // int32_t s_y = s_x

    void secTanh(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                 int32_t bw_y, int32_t s_x, int32_t s_y); // int32_t s_y = s_x

    void secSoftMax(int32_t dim, int32_t n_class, uint64_t *x, uint64_t *y, int32_t bw_x,
                    int32_t bw_y, int32_t s_x, int32_t s_y); // int32_t s_y = s_x

    void secNewtonRaphson(int32_t dim,
                          // Current implementation assumes that x is always of the form 1.y1y2y3..yn
                          uint64_t *x,
                          uint64_t *y,
                          // bitwidths
                          int32_t bw_x, int32_t bw_y,
                          // scales
                          int32_t s_x, int32_t s_y);

    // void secSigmoid2(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
    //                  int32_t bw_y, int32_t s_x);
};

#endif
