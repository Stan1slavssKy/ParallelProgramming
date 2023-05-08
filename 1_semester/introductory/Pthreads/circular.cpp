#include <pthread.h>
#include <iostream>
#include <vector>

constexpr size_t THREAD_NUMBER = 10;

using task_id_t = size_t;

static int value_to_increment = 0;
static task_id_t current_tid = 0;
static pthread_mutex_t mutex;
static pthread_cond_t cv;

void *increment(void *arg)
{
    task_id_t task_id = *((task_id_t *)arg);
    
    pthread_mutex_lock(&mutex);

    while (task_id != current_tid) {
        pthread_cond_wait(&cv, &mutex);
    }
    ++value_to_increment;
    ++current_tid;

    std::cout << "thread " << task_id << " increment value, " << "current_value = " << value_to_increment << std::endl;
    
    pthread_cond_broadcast(&cv);
    pthread_mutex_unlock(&mutex);
    pthread_exit(nullptr);
}

int main(int argc, char **argv)
{
    std::vector<pthread_t> threads(THREAD_NUMBER);
    std::vector<task_id_t> tasks_id(THREAD_NUMBER);

    for (size_t i = 0; i < THREAD_NUMBER; ++i) {
        pthread_t thread;
        tasks_id[i] = i;
        if (pthread_create(&thread, nullptr, increment, &tasks_id[i])) {
            std::cerr << "[ERROR:] Can't create thread" << std::endl;
        }

        threads[i] = thread;
    }

    for (size_t i = 0; i < THREAD_NUMBER; ++i) {
        pthread_join(threads[i], nullptr);
    }

    return 0;
}
