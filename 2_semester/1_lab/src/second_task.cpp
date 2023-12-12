#include <omp.h>
#include <cmath>
#include <unistd.h>
#include <string>
#include <cstring>

#include "array_file_output.h"
#include "scoped_time_measure.h"

static constexpr size_t ISIZE = 15000;
static constexpr size_t JSIZE = ISIZE;

void DumpToFile(const char *filename, double *array);
int GetThreadNumber(int argc, char** argv);

int main(int argc, char **argv) {
    int thread_number = GetThreadNumber(argc, argv);
    if (thread_number == -1) {
        return 1;
    }
    omp_set_num_threads(thread_number);

    double *main_array = new double[ISIZE * JSIZE];

    for (size_t i = 0; i < ISIZE; ++i) {
        for (size_t j = 0; j < JSIZE; ++j) {
            main_array[JSIZE * i + j] = 10 * i + j;
        }
    }

    double start_wtime = omp_get_wtime();

    size_t distance = ISIZE / thread_number;
    size_t remainder = ISIZE % thread_number;

#pragma omp parallel
    {
        int rank = omp_get_thread_num();
        size_t start = distance * rank;
        size_t end = distance * (rank + 1);

        if (remainder != 0) {
            if (remainder > static_cast<size_t>(rank)) {
                start += rank;
                end += rank + 1;
            }
            else {
                start += remainder;
                end += remainder;
            }
        }

        double *array = nullptr;
        if (rank != thread_number - 1) {
            array = new double[4 * JSIZE];

            size_t fetch_idx = end;
            size_t fetch_size = std::min(static_cast<size_t>(4), (ISIZE - fetch_idx + 1));

            std::memcpy(array, main_array + JSIZE * fetch_idx, fetch_size * JSIZE * sizeof(double));
        }
#pragma omp barrier

        if (rank != thread_number - 1) {
            for (size_t i = start; i < end - 4; ++i) {
                for (size_t j = 0; j < JSIZE - 2; ++j) {
                    main_array[JSIZE * i + j] = std::sin(0.1 * main_array[JSIZE * (i + 4) + (j + 2)]);
                }
            }
            for (size_t i = end - 4; i < end; ++i) {
                for (size_t j = 0; j < JSIZE - 2; ++j) {
                    main_array[JSIZE * i + j] = std::sin(0.1 * array[JSIZE * (i - end + 4) + (j + 2)]);
                }
            }
        } else {
            for (size_t i = start; i < ISIZE - 4; ++i) {
                for (size_t j = 0; j < JSIZE - 2; ++j) {
                    main_array[JSIZE * i + j] = std::sin(0.1 * main_array[JSIZE * (i + 4) + (j + 2)]);
                }
            }
        }
        delete[] array;
    }

    double end_wtime = omp_get_wtime();

    std::cout << "exec_time: " << end_wtime - start_wtime << std::endl;
    // DumpToFile(("second_task_openmp_" + std::to_string(thread_number) + ".txt").c_str(), main_array);

    delete[] main_array;

    return 0;
}

void DumpToFile(const char *filename, double *array)
{
    FILE *ff = nullptr;
    ff = fopen(filename, "w");

    for (size_t i = 0; i < ISIZE; ++i) {
        for (size_t j = 0; j < JSIZE; ++j) {
            fprintf(ff, "%f ", array[JSIZE * i + j]);
        }
        fprintf(ff, "\n");
    }
    fclose(ff);
}

int GetThreadNumber(int argc, char** argv)
{
    if (argc > 2) {
        std::cerr << "Too many arguments" << std::endl;
        return -1;
    }
    if (argc < 2) {
        std::cerr << "Too few arguments" << std::endl;
        return -1;
    }

    int thread_number = std::strtol(argv[argc - 1], nullptr, 10);
    if (errno == ERANGE) {
        std::cerr << "errno = STRTOL_ERROR" << std::endl;
        return -1;
    }
    if (thread_number < 1) {
        std::cerr << "Invalid thread_number = " << thread_number << ". Thread number should be >= 1" << std::endl;
        return -1;
    }

    return thread_number;
}
