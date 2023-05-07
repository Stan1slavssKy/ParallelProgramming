#include <pthread.h>
#include <iostream>
#include <vector>

constexpr size_t THREAD_NUMBER = 10;

using task_id_t = size_t;

void *hello_message(void *arg)
{
    printf("Hello world! [task_id = %ld out of %ld processors]\n", *((task_id_t*)arg), THREAD_NUMBER);
    return nullptr;
}

int main(int argc, char **argv)
{
    std::vector<pthread_t> threads(THREAD_NUMBER);
    std::vector<task_id_t> tasks_id(THREAD_NUMBER);

    for (size_t i = 0; i < THREAD_NUMBER; ++i) {
        pthread_t thread;
        tasks_id[i] = i;
        if (pthread_create(&thread, nullptr, hello_message, &tasks_id[i])) {
            std::cerr << "[ERROR:] Can't create thread" << std::endl;
        }

        threads[i] = thread;
    }

    for (size_t i = 0; i < THREAD_NUMBER; ++i) {
        pthread_join(threads[i], nullptr);
    }

    return 0;
}