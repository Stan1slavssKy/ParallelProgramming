#ifndef LAB_2_SRC_INTEGRATION_H
#define LAB_2_SRC_INTEGRATION_H

#include <cstdio>
#include <chrono>
#include <iostream>
#include <cassert>
#include <pthread.h>
#include <vector>

struct CLIArgs {
    size_t thread_number {0};
    double accuracy {0};
};

class WorkerThread {
public:
    static constexpr double interval_begin = 0;
    static constexpr double interval_end = 0.99;
    static constexpr double integration_segment = interval_end - interval_begin;

    static double function(double x)
    { 
        return 1.0 / ((x - 1) * (x - 1));
    }

    WorkerThread() = default;
    ~WorkerThread() = default;

    explicit WorkerThread(size_t task_id, size_t thread_number, double accuracy)
    {
        task_id_ = task_id;
        thread_number_ = thread_number;
        accuracy_ = accuracy;
        CalculateStartAndEndPos();
    }

    void CalculateStartAndEndPos()
    {
        if (thread_number_ == 1) {
            double distance = integration_segment / static_cast<double>(thread_number_);
            start_pos_ = interval_begin + distance * task_id_;
            end_pos_ = interval_begin + distance * (task_id_ + 1);
            return;
        }
        if (task_id_ == 0) {
            start_pos_ = interval_begin;
            end_pos_ = interval_end / 1.8;
        } else {
            double distance = (integration_segment - interval_end / 1.8) / static_cast<double>(thread_number_ - 1);
            start_pos_ = interval_begin + interval_end / 1.8 + distance * (task_id_ - 1);
            end_pos_ = interval_begin + interval_end / 1.8 + distance * (task_id_ );
        }
    }

    static void * Integrate(void *wrkr);
    static double DoIntegrate(double a, double b, double accuracy);
    static double DoSimpsonIteration(double a, double b);
    static double RunSimpsonMethod(double from, double to, double step);

    size_t task_id_ {0};
    size_t thread_number_ {0};

    double accuracy_ {0};
    double start_pos_ {0};
    double end_pos_ {0};
    double part_result_value_ {0};
};

class MainThread {
public:
    explicit MainThread(CLIArgs *cli_args)
    {
        assert(cli_args != nullptr);

        thread_number_ = cli_args->thread_number;
        accuracy_ = cli_args->accuracy;

        threads_.resize(thread_number_);
        workers_.resize(thread_number_);
    }
    ~MainThread() {}

    bool CreateThreads();
    void JoinThreads();

    double GetResult() 
    {
        return result_value_;
    }

    static int ParseCLIArgs(int argc, char **argv, CLIArgs* cli_args);

private:
    size_t thread_number_ {0};
    double accuracy_ {0};
    double result_value_ {0};

    std::vector<pthread_t> threads_;
    std::vector<WorkerThread> workers_;
};

using chr_time_p = std::chrono::time_point<std::chrono::high_resolution_clock>;

class ScopedTimeMeasure {
public:
    explicit ScopedTimeMeasure()
    {
        start_ = std::chrono::high_resolution_clock::now();
    }
    ~ScopedTimeMeasure()
    {
        end_ = std::chrono::high_resolution_clock::now();
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(end_ - start_).count();
        std::cout << "Elapsed time " << us << "us" << std::endl;
    }

private:
    chr_time_p start_;
    chr_time_p end_;
};

double runSimpson(double from, double to, double step);
double doSimpsonIter(double a, double b);
double f(double x);

#endif  // LAB_2_SRC_INTEGRATION_H
