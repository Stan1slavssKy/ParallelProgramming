#include <cmath>
#include <stack>
#include <unistd.h>
#include "integration.h"

namespace Global {

    static std::stack<IntegrateTaskInfo> stack;

    static size_t active_threads = 0;

    // lock if global stack is empty
    static pthread_mutex_t stack_present_lock;
    static pthread_mutex_t stack_access_lock;

    static constexpr int local_to_global_stack_count = 15;
    static constexpr int max_stack_size = 50;

}  // namespace Global

#define PRINT_RESULT_INTEGRAL_OFF

int main(int argc, char **argv)
{
    CLIArgs cli_args;
    if (MainThread::ParseCLIArgs(argc, argv, &cli_args) < 0) {
        std::cerr << "FAILED TO PARSE CLI ARGS" << std::endl;
        return 1;
    }
    
    MainThread main_thread(&cli_args);
    {
        ScopedTimeMeasure stm;
        main_thread.CreateThreads();
        main_thread.JoinThreads();
    }

#ifdef PRINT_RESULT_INTEGRAL_ON
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    std::cout << "I = " << main_thread.GetResult() << std::endl;
#endif

    return 0;
}

int MainThread::ParseCLIArgs(int argc, char **argv, CLIArgs* cli_args)
{
    if (argc != 3) {
        std::cerr << "USAGE: [THREAD_NUM] [ACCURACY]" << std::endl;
        return -1;
    }

    cli_args->thread_number = std::atoi(argv[1]);
    if (cli_args->thread_number < 1) {
        std::cerr << "Incorrect threads amount: " << cli_args->thread_number << "thread_number should be >= 1" << std::endl;
        return -1;
    }

    cli_args->accuracy = std::abs(std::atof(argv[2]));

    return 0;
}

bool MainThread::CreateThreads()
{
    double f_a = WorkerThread::function(interval_begin);
    double f_b = WorkerThread::function(interval_end);

    Global::stack.push(IntegrateTaskInfo{interval_begin, interval_end, f_a, f_b, (f_a + f_b) * (interval_end - interval_begin) / 2});

    for (size_t i = 0; i < thread_number_; ++i) {
        pthread_t thread;
        WorkerThread worker(i, thread_number_, accuracy_);
        workers_[i] = worker;

        if (pthread_create(&thread, nullptr, WorkerThread::Integrate, &workers_[i])) {
            std::cerr << "[ERROR:] Can't create thread" << std::endl;
            return false;
        }
        threads_[i] = thread;
    }
    return true;
}

void MainThread::JoinThreads()
{
    void *status = nullptr;
    for (size_t i = 0; i < thread_number_; ++i) {
        pthread_join(threads_[i], &status);
        WorkerThread *worker = reinterpret_cast<WorkerThread *>(status);
        result_value_ += worker->part_result_value_;
    }
}

/* static */
void *WorkerThread::Integrate(void *wrkr)
{
    WorkerThread *worker = (WorkerThread *)wrkr;

    // Architecture specific mapping on cpu's
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(worker->task_id_, &cpuset);
    pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);

    double a = 0;
    double b = 0;
    double f_a = 0;
    double f_b = 0;
    double int_ab = 0;

    while (true) {
        pthread_mutex_lock(&Global::stack_present_lock);
        pthread_mutex_lock(&Global::stack_access_lock);

        IntegrateTaskInfo new_task = Global::stack.top();

        a = new_task.start;
        b = new_task.end;
        f_a = new_task.f_start;
        f_b = new_task.f_end;
        int_ab = new_task.area;

        Global::stack.pop();

        if (!Global::stack.empty()) {
            pthread_mutex_unlock(&Global::stack_present_lock);
        }
        if (a <= b) {
            ++Global::active_threads;
        }

        pthread_mutex_unlock(&Global::stack_access_lock);

        if (a > b) {
            // at the end we push a special task in which a > b
            break;
        }

        while (true) {
            const double c = (a + b) / 2;
            const double f_c = WorkerThread::function(c);

            const double int_ac = (f_a + f_c) * (c - a) / 2;
            const double int_cb = (f_c + f_b) * (b - c) / 2;
            const double int_acb = int_ac + int_cb;

            if (std::abs(int_ab - int_acb) >= worker->accuracy_ * std::abs(int_acb)) {
                // push in stack left part
                worker->local_stack_.push(IntegrateTaskInfo{a, c, f_a, f_c, int_ac});
                // start next iteration with right part
                a = c;
                f_a = f_c;
                int_ab = int_cb;
            } else {
                worker->part_result_value_ += int_acb;
                if (worker->local_stack_.empty()) {
                    break;
                }

                IntegrateTaskInfo new_task = worker->local_stack_.top();

                a = new_task.start;
                b = new_task.end;
                f_a = new_task.f_start;
                f_b = new_task.f_end;
                int_ab = new_task.area;

                worker->local_stack_.pop();
            }

            if (worker->local_stack_.size() > Global::local_to_global_stack_count) {
                pthread_mutex_lock(&Global::stack_access_lock);
                if (Global::stack.empty()) {
                    while (worker->local_stack_.size() > 1 && Global::stack.size() < Global::max_stack_size) {
                        Global::stack.push(worker->local_stack_.top());
                        worker->local_stack_.pop();
                    }
                    pthread_mutex_unlock(&Global::stack_present_lock);
                }
                pthread_mutex_unlock(&Global::stack_access_lock);
            }
        }
        
        pthread_mutex_lock(&Global::stack_access_lock);
        --Global::active_threads;

        if (Global::active_threads == 0 && Global::stack.empty()) {
            for(size_t i = 0; i < worker->thread_number_; ++i) {
                // push a special task in which a > b for break out of the loop
                Global::stack.push(IntegrateTaskInfo{2.0, 1.0, 0.0, 0.0, 0.0});
            }
            pthread_mutex_unlock(&Global::stack_present_lock);
        }

        pthread_mutex_unlock(&Global::stack_access_lock);
    }
    pthread_exit(worker);
}
