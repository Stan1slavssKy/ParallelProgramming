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

struct IntegrateTaskInfo {
    double start;
    double end;

    double f_start;
    double f_end;

    double area;
};

class WorkerThread {
public:
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
    }

    static void * Integrate(void *wrkr);

    size_t task_id_ {0};
    size_t thread_number_ {0};

    double accuracy_ {0};
    double part_result_value_ {0};

    std::stack<IntegrateTaskInfo> local_stack_;
};

class MainThread {
public:
    static constexpr double interval_begin = 0;
    static constexpr double interval_end = 0.99;
    static constexpr double integration_segment = interval_end - interval_begin;

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

class ScopedTimeMeasure {
public:
    using chr_time_p = std::chrono::time_point<std::chrono::high_resolution_clock>;

    explicit ScopedTimeMeasure()
    {
        start_ = std::chrono::high_resolution_clock::now();
    }
    ~ScopedTimeMeasure()
    {
        end_ = std::chrono::high_resolution_clock::now();
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(end_ - start_).count();
        std::cout << us << " us" << std::endl;
    }

private:
    chr_time_p start_;
    chr_time_p end_;
};

#endif  // LAB_2_SRC_INTEGRATION_H
