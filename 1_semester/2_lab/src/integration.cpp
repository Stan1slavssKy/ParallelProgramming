#include <cmath>
#include "integration.h"

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

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    std::cout << "I = " << main_thread.GetResult() << std::endl;

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
static pthread_mutex_t mutex;
/* static */
void *WorkerThread::Integrate(void *wrkr)
{
    WorkerThread *worker = (WorkerThread *)wrkr;
        // pthread_mutex_lock(&mutex);
        // std::cout << worker->start_pos_ << " " << worker->end_pos_ << std::endl;
        // pthread_mutex_unlock(&mutex);
    worker->part_result_value_ = WorkerThread::DoIntegrate(worker->start_pos_, worker->end_pos_, worker->accuracy_);
    pthread_exit(worker);
}

/* static */
double WorkerThread::DoIntegrate(double a, double b, double accuracy)
{
    double res = 0;
    double c = (a + b) / 2;
    double f_a = WorkerThread::function(a);
    double f_b = WorkerThread::function(b);
    double f_c = WorkerThread::function(c);

    double int_ab = (f_a + f_b) * (b - a) / 2;
    double int_ac = (f_a + f_c) * (c - a) / 2;
    double int_cb = (f_c + f_b) * (b - c) / 2;
    
    double int_acb = int_ac + int_cb;
    
    if (std::abs(int_ab - int_acb) >= accuracy * std::abs(int_acb)) {
        res = DoIntegrate(a, c, accuracy) + DoIntegrate(c, b, accuracy);
    } else {
        res = int_acb;
    }
    return res;
}

/* static */
double WorkerThread::DoSimpsonIteration(double a, double b)
{
    return (b - a) / 6 * (WorkerThread::function(a) + 4 * WorkerThread::function((a + b) / 2) + WorkerThread::function(b));
}

/* static */
double WorkerThread::RunSimpsonMethod(double from, double to, double step)
{
    double result = 0;
    for (double point = from; point < to; point += step) {
        result += DoSimpsonIteration(point, point + step);
    }
    return result;
}
