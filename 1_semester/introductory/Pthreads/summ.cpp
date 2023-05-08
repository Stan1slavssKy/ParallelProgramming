#include <pthread.h>
#include <iostream>
#include <vector>

constexpr size_t THREAD_NUMBER = 5;

using task_id_t = size_t;

struct Worker {
    int start {0};
    int end {0};
    double result {0};
};

void *CalculatePartSumm(void *arg)
{
    Worker *pos = (Worker *)arg;
    for (int i = pos->start; i < pos->end; ++i) {
        pos->result = pos->result + 1.0 / i;
    }
    if (pos->start == 1) {
        // zero (master) worker
        return nullptr;
    } else {
        pthread_exit(pos);
    }
}

int GetSumElementsNumber(int argc, char** argv)
{
    if (argc != 2) {
        std::cout << "Error. Please pass the parameter at startup" << std::endl;
        return -1;
    }

    int N = std::atoi(argv[argc - 1]);
    if (N == 0) {
        std::cout << "You pass into arg N = 0 or atoi conversion can't be performed" << std::endl;
        return -1;
    }

    return N;
}

int main(int argc, char **argv)
{
    int N = GetSumElementsNumber(argc, argv);
    if (N < 0) {
        return 1;
    }
    int distance = N / THREAD_NUMBER;
    if (distance == 0) {
        std::cerr << "Too many threads or N too low" << std::endl;
        return 1;
    }
    int remainder = N % THREAD_NUMBER;
    
    std::vector<pthread_t> threads(THREAD_NUMBER);
    std::vector<Worker> workers(THREAD_NUMBER);

    for (size_t i = 0; i < THREAD_NUMBER; ++i) {
        pthread_t thread;
        Worker worker;

        int start = 1 + distance * i;
        int end = 1 + distance * (i + 1);

        // If the remainder is not zero, we add 1 term of the sum to each process starting from zero
        if (remainder != 0) {
            if (remainder - i > 0) {
                start += i;
                end += i + 1;
            }
            else {
                start += remainder;
                end += remainder;
            }
        }

        worker.start = start;
        worker.end = end;
        workers[i] = worker;

        if (i == 0) {
            CalculatePartSumm(&workers[i]);
        } else {
            if (pthread_create(&thread, nullptr, CalculatePartSumm, &workers[i])) {
                std::cerr << "[ERROR:] Can't create thread" << std::endl;
            }
        }
        
        threads[i] = thread;
    }

    double result = workers[0].result;
    void *status = nullptr;
    for (size_t i = 1; i < THREAD_NUMBER; ++i) {
        pthread_join(threads[i], &status);
        result += ((Worker *)status)->result;
    }

    std::cout << "[SUMM = " << result << "]"<< std::endl;

    return 0;
}