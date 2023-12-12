#include <omp.h>
#include <cmath>
#include <unistd.h>
#include <string>

#include "array_file_output.h"
#include "scoped_time_measure.h"

static constexpr size_t ISIZE = 15000;
static constexpr size_t JSIZE = ISIZE;

int GetThreadNumber(int argc, char** argv);

int main(int argc, char **argv) {
    int thread_number = GetThreadNumber(argc, argv);
    if (thread_number == -1) {
        return 1;
    }
    omp_set_num_threads(thread_number);

    ArrayHandler array_handler(ISIZE, JSIZE);
    array_handler.DefaultFillIn();

    auto *array = array_handler.GetArray();

    double start_wtime = omp_get_wtime();

#pragma omp parallel for
    for (size_t i = 0; i < ISIZE; ++i) {
        for (size_t j = 0; j < JSIZE; ++j) {
            array[JSIZE * i + j] = std::sin(2 * array[JSIZE * i + j]);
        }
    }

    double end_wtime = omp_get_wtime();

    std::cout << "exec_time: " << end_wtime - start_wtime << std::endl;
    // array_handler.DumpToFile(("baseline_openmp_" + std::to_string(thread_number) + ".txt").c_str());

    return 0;
}

void DumpToFile(const char *filename, double *array)
{
    FILE *ff = nullptr;
    ff = fopen(filename, "w");

    for (size_t i = 0; i < ISIZE; ++i) {
        for (size_t j = 0; j < JSIZE; ++j) {
            fprintf(ff, "%f", array[JSIZE * i + j]);
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