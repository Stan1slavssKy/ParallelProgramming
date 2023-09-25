#include <omp.h>
#include <iostream>

int main() {
    omp_set_num_threads(16);
    
    int counter = 0;

    #pragma omp parallel shared(counter)
    {
        #pragma omp for ordered
        for (int i = 0; i < omp_get_num_threads(); ++i) {
            #pragma omp ordered
            {
                std::cout << "Incrementing from thread = " << omp_get_thread_num() << "; Counter value = " 
                                                                                   << ++counter << std::endl;
            }
        }
    }

    return 0;
}
