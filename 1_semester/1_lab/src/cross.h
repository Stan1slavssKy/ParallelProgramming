#ifndef CROSS_H
#define CROSS_H

#include <array>
#include <vector>
#include <cmath>

namespace Equation {

constexpr double a = 1.0;
constexpr double X = 50.0;
constexpr double T = 1.0;

double func(double x, double t)
{ 
    return std::exp(std::sin(x * t / X));
}

double phi(double x)
{
    return std::cos(M_PI * x / X);
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
constexpr size_t K = Equation::T / tau;
constexpr size_t M = Equation::X / h;

/**
 * Worker created for each process
 * 
 * For each process allocated array named values_, that aim is to store its part of the calculations
 * 
 * In root_process result_values_ is allocated to gather all work (in the form of values_) 
 * from every proccess after main algorithms is ends
 * 
 * Description of the algorithm:
 * 
 * Firstly we need to fill initial conditions using functions psi and phi
 * The 2 step with the use of the rectangle method is to fill in the second layer
 * We use the rectangle in order not to lose the accuracy of the Cross algorithm
 * 
 * 0 0 0 0 0   (1)   + 0 0 0 0   (2)   + 0 0 0 0
 * 0 0 0 0 0  ---->  + 0 0 0 0  ---->  + 0 0 0 0 
 * 0 0 0 0 0         + 0 0 0 0         + + + + +
 * 0 0 0 0 0         + + + + +         + + + + +
 * 
 * Now let's take a closer look at the parallelization algorithm:
 * 
 * The algorithm consists in splitting the field into n equal parts, where n is the number of processes
 * The problem of Cross is that in boundary parts we cant calculate values. Suggestion sulution is to 
 * use Bsend to send left and rigth boundary value of process k. And in the last worker last boundary 
 * value is calculating by rectangle method
*/

class Worker {
public:
    Worker(int rank, int commsize)
    {
        rank_ = rank;
        commsize_ = commsize;

        CalculateStartAndEndPos();
        values_ = new double[cols_ * K_] {};

        if (rank_ == 0) {
            if (commsize_ > 1) {
                result_values_ = new double[M_ * K_] {}; 
            } else {
                result_values_ = values_;
            }
        }
    }

    ~Worker()
    {
        delete[] values_;
        values_ = nullptr;

        if (rank_ == 0 && commsize_ > 1) {
            delete[] result_values_;
            result_values_ = nullptr;
        }
    }

    void FillInitialConditions();
    void CalculateFirstLineByRectangleMethod();
    void CalculateOtherLinesByCrossMethod();
    void GatheringAllWork();
    void PrintResult();

private:
    size_t rank_ {0};
    size_t commsize_ {0};
    size_t start_pos_ {0};
    size_t end_pos_ {0};
    size_t cols_ {0};

    size_t M_ {M};
    size_t K_ {K};

    double* values_ {nullptr};
    double* result_values_ {nullptr};

    void CalculateStartAndEndPos();
    /**   Example for M = 5
     *
     *    [0, 2), [2, 4), [4, 5)
     *     |                  |
     *     start_pos          end_pos
     *   
    **/
};

#endif  // CROSS_H
