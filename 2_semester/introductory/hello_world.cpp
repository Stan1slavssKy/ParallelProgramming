#include <omp.h>
#include <stdio.h>

int main() {
    omp_set_num_threads(16);

#pragma omp parallel
{
    printf("Hello World from thread = %d; Thread number = %d.\n", omp_get_thread_num(), omp_get_num_threads());
}

    return 0;
}
