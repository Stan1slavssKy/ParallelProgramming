#include <iostream>
#include <cmath>
#include <mpi.h>
#include <fstream>
#include <cassert>

#include "cross.h"

int main(int argc, char** argv) {
    MPI::Init(argc, argv);

    auto rank = MPI::COMM_WORLD.Get_rank();
    auto commsize = MPI::COMM_WORLD.Get_size();

    int distance = M / commsize;
    int remainder = M % commsize;
    int start = 1 + distance * rank;
    int end = 1 + distance * (rank + 1);

    // If the remainder is not zero, we add 1 term of the sum to each process starting from zero
    if (remainder != 0) {
        if (remainder - rank > 0) {
            start += rank;
            end += rank + 1;
        }
        else {
            start += remainder;
            end += remainder;
        }
    }

    // made constructor that allocate memory for different rank
    Solution result {};
    
    FillInitialConditions(&result);
    CalculateFirstLineByRectangleMethod(&result);
    CalculateOtherLines(&result);

    if (rank == 0) {
        std::string name = "res.txt";
        std::ofstream f(name);
        PrintRes(f, &result);
    }

    MPI::Finalize();
    return 0;
}

void FillInitialConditions(Solution* result)
{
    for (int i = 0; i < K; ++i) {
        result->t[i] = tau * i;
        result->values[0 + i * M] = Equation::psi(tau * i);
    }
    for (int i = 0; i < M; ++i) {
        result->x[i] = h * i;
        result->values[i] = Equation::phi(h * i);
    }
}

void CalculateFirstLineByRectangleMethod(Solution* result)
{
    double* values = result->values;

    for (int m = 1; m < M; ++m) {
        double first_part = (values[M + m - 1] - values[m - 1] - values[M + m]) / (2 * tau);
        double second_part = Equation::a * (-1.0 * values[M + m - 1] + values[m]  - values[m - 1]) / (2 * h);
        double func_part = Equation::func((m + 0.5) * h, (1 + 0.5) * tau);
        values[M + m] = (func_part - first_part - second_part) * 2 * tau * h / (tau + h);
    }
}

void CalculateOtherLines(Solution* result)
{
    double* values = result->values;

    for (int k = 1; k < K - 1; ++k) {
        for (int m = 1; m < M - 1; ++m) {
            double first_part = (-1.0 * values[M * (k - 1) + m]) / (2 * tau);
            double second_part = Equation::a * (values[M * k + m + 1] - values[M * k + m - 1]) / (2 * h);
            double func_part = Equation::func(m * h, k * tau);
            values[M * (k + 1) + m] = (func_part - first_part - second_part) * 2 * tau;
        }
        // code for m = M - 1
        double first_part = (values[M * (k + 1) + M - 1 - 1] - values[M * k + M - 1 - 1] - values[M * k + M - 1]) / (2 * tau);
        double second_part = Equation::a * (-1.0 * values[M * (k + 1) + M - 1 - 1] + values[M * k + M - 1] - values[M * k + M - 1 - 1]) / (2 * h);
        double func_part = Equation::func((M - 1 + 0.5) * h, (k + 0.5) * tau);
        values[M * (k + 1) + M - 1] = (func_part - first_part - second_part) * 2 * tau * h / (tau + h);
    }
}

void PrintRes(std::ostream &ost, Solution* result)
{
    ost << "tau: " << tau << std::endl;
    ost << "h: " << h << std::endl;
    ost << "T: " << Equation::T << std::endl;
    ost << "X: " << Equation::X << std::endl;
    ost << "M: " << M << std::endl;
    ost << "K: " << K << std::endl;

    for (std::size_t i = 0, sz = K * M; i < sz; ++i)
        ost << result->values[i] << std::endl;
}
