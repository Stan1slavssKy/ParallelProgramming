#include <iostream>
#include <cmath>
#include <mpi.h>
#include <fstream>
#include <cassert>
#include <memory.h>
#include <unistd.h>

#include "cross.h"

#define ENABLE_SAVING_PICTURE_OFF

int main(int argc, char** argv) {
    MPI::Init(argc, argv);

    auto rank = MPI::COMM_WORLD.Get_rank();
    auto commsize = MPI::COMM_WORLD.Get_size();

    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(2 * rank, &mask);
    sched_setaffinity(getpid(), sizeof(cpu_set_t), &mask);

    double start_wtime = MPI::Wtime();

    Worker result(rank, commsize);
    
    result.FillInitialConditions();
    result.CalculateFirstLineByRectangleMethod();
    result.CalculateOtherLinesByCrossMethod();

    if (commsize > 1) {
        result.GatheringAllWork();
    }

    double end_wtime = MPI::Wtime();

    if (rank == 0) {
        std::cout << "exec_time: " << end_wtime - start_wtime << std::endl;
        #ifdef ENABLE_SAVING_PICTURE
        result.PrintResult();
        #endif
    }

    MPI::Finalize();
    return 0;
}

void Worker::FillInitialConditions()
{
    if (rank_ == 0) {
        for (size_t i = 0; i < K_; ++i) {
            values_[0 + i * cols_] = Equation::psi(tau * i);
        }
    }
    for (size_t i = 0; i < cols_; ++i) {
        values_[i] = Equation::phi(h * (i + start_pos_));
    }
}

void Worker::CalculateFirstLineByRectangleMethod()
{
    if (rank_ == 0) {
        for (size_t m = 1; m < cols_; ++m) {
            double first_part = (values_[cols_ + m - 1] - values_[m - 1] - values_[m]) / (2 * tau);
            double second_part = Equation::a * (-1.0 * values_[cols_ + m - 1] + values_[m]  - values_[m - 1]) / (2 * h);
            double func_part = Equation::func((m + 0.5) * h, (1 + 0.5) * tau);
            values_[cols_ + m] = (func_part - first_part - second_part) * 2 * tau * h / (h + Equation::a * tau);
        }
        if (commsize_ > 1) {
            MPI_Send(values_ + cols_ + cols_ - 1, 1, MPI::DOUBLE, rank_ + 1, 0, MPI_COMM_WORLD);
            MPI_Send(values_ + cols_ - 1, 1, MPI::DOUBLE, rank_ + 1, 0, MPI_COMM_WORLD);
        }
        return;
    }

    double recv_up_value = 0;
    double recv_down_value = 0;
    MPI_Status status;
    MPI_Recv(&recv_up_value, 1, MPI::DOUBLE, rank_ - 1, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&recv_down_value, 1, MPI::DOUBLE, rank_ - 1, 0, MPI_COMM_WORLD, &status);
    
    for (size_t m = 0; m < cols_; ++m) {
        if (m != 0) {
            recv_up_value = values_[cols_ + m - 1];
            recv_down_value = values_[m - 1];
        }
        double first_part = (recv_up_value - recv_down_value - values_[m]) / (2 * tau);
        double second_part = Equation::a * (-1.0 * recv_up_value + values_[m]  - recv_down_value) / (2 * h);
        double func_part = Equation::func((start_pos_ + m + 0.5) * h, (1 + 0.5) * tau);
        values_[cols_ + m] = (func_part - first_part - second_part) * 2 * tau * h / (h + Equation::a * tau);
    }
    if (rank_ != commsize_ - 1) {
        MPI_Send(values_ + cols_ + cols_ - 1, 1, MPI::DOUBLE, rank_ + 1, 0, MPI_COMM_WORLD);
        MPI_Send(values_ + cols_ - 1, 1, MPI::DOUBLE, rank_ + 1, 0, MPI_COMM_WORLD);
    }
}

void Worker::CalculateOtherLinesByCrossMethod()
{
    for (size_t k = 1; k < K_ - 1; ++k) {
        double recv_value = 0; 
        if (rank_ != 0) {
            MPI_Status status;
            MPI_Recv(&recv_value, 1, MPI::DOUBLE, rank_ - 1, 0, MPI_COMM_WORLD, &status);
        }
        if (rank_ != commsize_ - 1) {
            // send last value in already filled layer
            MPI_Send(values_ + cols_ * k + cols_ - 1, 1, MPI::DOUBLE, rank_ + 1, 0, MPI_COMM_WORLD);
        }

        // if rank = 0 we dont need to calculate m = 0
        for (size_t m = (rank_ == 0); m < cols_ - 1; ++m) {
            if (m != 0) {
                recv_value = values_[cols_ * k + m - 1];
            }
            double first_part = (-1.0 * values_[cols_ * (k - 1) + m]) / (2 * tau);
            double second_part = Equation::a * (values_[cols_ * k + m + 1] - recv_value) / (2 * h);
            double func_part = Equation::func((start_pos_ + m) * h, k * tau);
            values_[cols_ * (k + 1) + m] = (func_part - first_part - second_part) * 2 * tau;
        }
        size_t m = cols_ - 1;
        double first_part = (values_[cols_ * (k + 1) + m - 1] - values_[cols_ * k + m - 1] - values_[cols_ * k + m]) / (2 * tau);
        double second_part = Equation::a * (-1.0 * values_[cols_ * (k + 1) + m - 1] + values_[cols_ * k + m] - values_[cols_ * k + m - 1]) / (2 * h);
        double func_part = Equation::func((start_pos_ + m + 0.5) * h, (k + 0.5) * tau);
        values_[cols_ * (k + 1) + m] = (func_part - first_part - second_part) * 2 * tau * h / (h + Equation::a * tau);
    }
}

void Worker::PrintResult()
{
    std::string name = "./build/res.txt";
    std::ofstream ost(name);

    ost << "tau: " << tau << std::endl;
    ost << "h: " << h << std::endl;
    ost << "T: " << Equation::T << std::endl;
    ost << "X: " << Equation::X << std::endl;
    ost << "M: " << M_ << std::endl;
    ost << "K: " << K_ << std::endl;

    double* array_to_print = nullptr;
    if (commsize_ == 1) {
        array_to_print = values_;
    } else {
        array_to_print = result_values_;
    }

    for (size_t i = 0, size = K_ * M_; i < size; ++i) {
        ost << array_to_print[i] << std::endl;
    }
}

void Worker::CalculateStartAndEndPos()
{
    size_t distance = M_ / commsize_;
    size_t remainder = M_ % commsize_;
    if (remainder != 0) {
        M_ = (distance + 1) * commsize_;
        ++distance;
    }
    start_pos_ = distance * rank_;
    end_pos_ = distance * (rank_ + 1);

    cols_ = end_pos_ - start_pos_;
}

void Worker::GatheringAllWork()
{    
    for (size_t k = 0; k < K_; ++k) {
        MPI_Gather(values_ + k * cols_, cols_, MPI::DOUBLE, result_values_ + k * M_, cols_, MPI::DOUBLE, 0, MPI_COMM_WORLD);
    }
}
