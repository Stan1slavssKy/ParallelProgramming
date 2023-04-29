#ifndef CROSS_H
#define CROSS_H

#include <array>
#include <vector>

namespace Equation {

constexpr double a = 1.0;
constexpr double X = 10.0;
constexpr double T = 10.0;

double func(double x, double t)
{ 
    return x * t;
}

double phi(double x)
{
    return x;
}

double psi(double t)
{
    return -1.0 * t;
}

}; // namespace Equation

constexpr double h = 0.05;
constexpr double tau = 0.001;
constexpr int K = Equation::T / tau + 1;
constexpr int M = Equation::X / h + 1;

struct Solution {
    Solution()
    {
        values = new double[K * M];
        t = new double[K];
        x = new double[M];
    }
    ~Solution()
    {
        delete[] values;
        delete[] t;
        delete[] x;
    }

    double* values {nullptr};
    double* t {nullptr};
    double* x {nullptr};
};

void FillInitialConditions(Solution* result);
void CalculateFirstLineByRectangleMethod(Solution* result);
void CalculateOtherLines(Solution* result);
void PrintRes(std::ostream &ost, Solution* result);

#endif  // CROSS_H