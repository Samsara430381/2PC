#include "MyProtocol/MyProtocol.h"
#include <fstream>
#include <iostream>
#include <thread>

using namespace sci;
using namespace std;

#define MAX_THREADS 4
#define n 6

int party, port = 32000;
int num_threads = 4;
string address = "172.17.0.3";

int dim = exp10(n);
int bw_x = 14;
int bw_y = 14;
int s_x = 12;
int s_y = 12; //(bwC <= (bwA + bwB)) && (bwC >= bwA) && (bwC >= bwB)

uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));

IOPack *iopackArr[MAX_THREADS];
OTPack *otpackArr[MAX_THREADS];

uint64_t computeULPErr(double calc, double actual, int SCALE)
{
    int64_t calc_fixed = (double(calc) * (1ULL << SCALE));
    int64_t actual_fixed = (double(actual) * (1ULL << SCALE));
    uint64_t ulp_err = (calc_fixed - actual_fixed) > 0
                           ? (calc_fixed - actual_fixed)
                           : (actual_fixed - calc_fixed);
    if (ulp_err > 20)
        cout << "  ulp:" << ulp_err << "  calc:" << calc << "   actual:" << actual << endl;
    return ulp_err;
}

void exp_thread(int tid, uint64_t *tmp1, uint64_t *x, uint64_t *y, int num_exp)
{
    MyProtocol *math;
    if (tid & 1)
    {
        math = new MyProtocol(3 - party, iopackArr[tid], otpackArr[tid]);
    }
    else
    {
        math = new MyProtocol(party, iopackArr[tid], otpackArr[tid]);
    }
    // math->div(num_exp, tmp1, x, y, 2, bw_x, bw_y, 0, s_x, s_y, false, false);
    math->secNewtonRaphson(num_exp, x, y,  bw_x,  bw_y, s_x, s_y);
    delete math;
}

int main(int argc, char **argv)
{
    /************* Argument Parsing  ************/
    /********************************************/
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("N", dim, "Number of exponentiations");
    amap.arg("nt", num_threads, "Number of threads");
    amap.arg("ip", address, "IP Address of server (ALICE)");

    amap.parse(argc, argv);

    assert(num_threads <= MAX_THREADS);

    /********** Setup IO and Base OTs ***********/
    /********************************************/
    for (int i = 0; i < num_threads; i++)
    {
        iopackArr[i] = new IOPack(party, port + i, address);
        if (i & 1)
        {
            otpackArr[i] = new OTPack(iopackArr[i], 3 - party);
        }
        else
        {
            otpackArr[i] = new OTPack(iopackArr[i], party);
        }
    }
    std::cout << "All Base OTs Done" << std::endl;

    /************ Generate Test Data ************/
    /********************************************/
    PRG128 prg;

    uint64_t *tmp1 = new uint64_t[dim];
    uint64_t *x = new uint64_t[dim];
    uint64_t *y = new uint64_t[dim];

    prg.random_data(x, dim * sizeof(uint64_t));
    for (int i = 0; i < dim; i++)
    {
        tmp1[i] = (party == ALICE ? 1 : 0);
    }
    if (party == ALICE)
    {
        iopackArr[0]->io->send_data(x, dim * sizeof(uint64_t));
    }
    else
    {
        uint64_t *x0 = new uint64_t[dim];
        iopackArr[0]->io->recv_data(x0, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; i++)
        {
            // x is always negative
            x[i] = ((1ULL << (bw_x - 3)) + (x[i] & (mask_x >> 3))) - x0[i];
            // x is always 0
            // x[i] = - x0[i];
        }
        delete[] x0;
    }
    for (int i = 0; i < dim; i++)
    {
        x[i] &= mask_x;
    }

    /************** Fork Threads ****************/
    /********************************************/
    uint64_t total_comm = 0;
    uint64_t thread_comm[num_threads];

    uint64_t thread_num_rounds[num_threads];
    for (int i = 0; i < num_threads; i++)
    {
        thread_comm[i] = iopackArr[i]->get_comm();
        thread_num_rounds[i] = iopackArr[i]->get_rounds();
    }

    auto start = clock_start();
    std::thread exp_threads[num_threads];
    int chunk_size = dim / num_threads;
    for (int i = 0; i < num_threads; ++i)
    {
        int offset = i * chunk_size;
        int lnum_exp;
        if (i == (num_threads - 1))
        {
            lnum_exp = dim - offset;
        }
        else
        {
            lnum_exp = chunk_size;
        }
        exp_threads[i] =
            std::thread(exp_thread, i, tmp1 + offset, x + offset, y + offset, lnum_exp);
    }
    for (int i = 0; i < num_threads; ++i)
    {
        exp_threads[i].join();
    }
    long long t = time_from(start);

    for (int i = 0; i < num_threads; i++)
    {
        thread_comm[i] = iopackArr[i]->get_comm() - thread_comm[i];
        thread_num_rounds[i] = iopackArr[i]->get_rounds() - thread_num_rounds[i];
        total_comm += thread_comm[i];
    }

    /************** Verification ****************/
    /********************************************/
    if (party == ALICE)
    {
        iopackArr[0]->io->send_data(x, dim * sizeof(uint64_t));
        iopackArr[0]->io->send_data(y, dim * sizeof(uint64_t));
    }
    else
    { // party == BOB
        uint64_t *x0 = new uint64_t[dim];
        uint64_t *y0 = new uint64_t[dim];
        iopackArr[0]->io->recv_data(x0, dim * sizeof(uint64_t));
        iopackArr[0]->io->recv_data(y0, dim * sizeof(uint64_t));

        uint64_t total_err = 0;
        uint64_t max_ULP_err = 0;
        for (int i = 0; i < dim; i++)
        {
            double dbl_x = (unsigned_val(x0[i] + x[i], bw_x)) / double(1LL << s_x);
            double dbl_y = (unsigned_val(y0[i] + y[i], bw_y)) / double(1ULL << s_y);
            double rec_x = 1 / dbl_x;

            // cout << "\tdbl_x=" << dbl_x << "\n\tdbl_y=" << dbl_y << endl;
            // cout << "\trec(x)=" << rec_x << endl;

            uint64_t err = computeULPErr(dbl_y, rec_x, s_y);
            total_err += err;
            max_ULP_err = std::max(max_ULP_err, err);
        }

        cerr << "Average ULP error: " << total_err / dim << endl;
        cerr << "Max ULP error: " << max_ULP_err << endl;
        cerr << "Number of tests: " << dim << endl;

        delete[] x0;
        delete[] y0;
    }
    cout << "Number of Exp/s:\t" << (double(dim) / t) * 1e6 << std::endl;
    cout << "Exp Time\t" << t / (1000.0) << " ms" << endl;
    cout << "Exp Bytes Sent\t" << total_comm << " bytes" << endl;
    cout << "Num rounds: " << thread_num_rounds[0] << endl;
    /******************* Cleanup ****************/
    /********************************************/
    delete[] x;
    delete[] y;
    for (int i = 0; i < num_threads; i++)
    {
        delete iopackArr[i];
        delete otpackArr[i];
    }
}
