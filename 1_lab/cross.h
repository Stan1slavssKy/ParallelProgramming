#ifndef CROSS_H
#define CROSS_H

#include <array>
#include <vector>
#include <cmath>

namespace Equation {

constexpr double a = 1.0;
constexpr double X = 5000.0;
constexpr double T = 1.0;

double func(double x, double t)
{ 
    return std::exp(std::sin(x * t));
}

double phi(double x)
{
    return std::cos(M_PI * x);
}

double psi(double t)
{
    return std::exp(-t);
}

}; // namespace Equation

// For the method to be stable, Courant number must be less than 1
// a * tau / h < 1
constexpr double h = 0.05;
constexpr double tau = 0.001;
constexpr long int K = Equation::T / tau;
constexpr long int M = Equation::X / h;

class Solution {
public:
    Solution(int rank, int commsize)
    {
        rank_ = rank;
        commsize_ = commsize;

        CalculateStartAndEndPos();
        cols_ = end_pos_ - start_pos_;
        values_ = new double[cols_ * K_] {};

        if (rank_ == 0) {
           result_values_ = new double[M_ * K_] {}; 
        }
    }
    ~Solution()
    {
        delete[] values_;
        if (rank_ == 0) {
            delete[] result_values_;
        }
    }

    void FillInitialConditions();
    void CalculateFirstLineByRectangleMethod();
    void CalculateOtherLines();
    void CalculateOtherLinesSingleWorker();
    void CalculateOtherLinesMultiWorkers();
    void GatheringAllWork();
    void PrintRes(std::ostream &ost);

private:
    int rank_ {0};
    int commsize_ {0};
    int start_pos_ {0};
    int end_pos_ {0};
    int cols_ {0};

    int M_ {M};
    int K_ {K};

    double* values_ {nullptr};
    double* result_values_ {nullptr};

    void CalculateStartAndEndPos();
    /**   Example for M = 5
     *
     *    [0, 2), [2, 4), [4, 5)
     *     |                  |
     *     start_pos          end_pos
     *
     *
     *    0 + + + +
     *    0 0 0 0 +
     *    0 0 0 0 0
     *    0 0 0 0 0
     *    ---------
     *    0 1 2 3 4
     *            |
     *            If pos == M - 1 we should use Rectangle method
     *            So when end_pos == M, we can be sure that it is the last worker
     *   
    **/     
};

#endif  // CROSS_H
