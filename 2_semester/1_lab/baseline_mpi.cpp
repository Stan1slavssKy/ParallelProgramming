#include <mpi.h>
#include <cmath>
#include <unistd.h>

#include "array_file_output.h"
#include "scoped_time_measure.h"

static constexpr size_t ISIZE = 5000;
static constexpr size_t JSIZE = ISIZE;

void MPIParallelization(double *current_array, int number_rows);
void GatheringAllWork(double *current_array, double *main_array, int number_rows, int rank);
void DumpToFile(const char *filename, double *array);

int main(int argc, char **argv) {
    MPI::Init(argc, argv);

    auto rank = MPI::COMM_WORLD.Get_rank();
    auto commsize = MPI::COMM_WORLD.Get_size();

    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(2 * rank, &mask);
    sched_setaffinity(getpid(), sizeof(cpu_set_t), &mask);
    
    double *main_array = nullptr;
    if (rank == 0) {
        main_array = new double[ISIZE * JSIZE];
    }

    int distance = ISIZE / commsize;
    int remainder = ISIZE % commsize;
    int start = distance * rank;
    int end = distance * (rank + 1);

    // If the remainder is not zero, we add 1 term of the sum to each process starting from zero
    if (remainder != 0) {
        if (remainder - rank > 0) {
            start += rank;
            end += rank + 1;
        }
        else {
            start += remainder;
            end += remainder;
        }
    }

    int number_rows = end - start;

    // std::cout << "[" << start << " " << end << ")" << std::endl;

    double *current_array = new double[number_rows * JSIZE];
    for (int i = 0; i < number_rows; ++i) {
        for (int j = 0; j < JSIZE; ++j) {
            current_array[JSIZE * i + j] = 10 * (start + i) + j;
        }
    }

    double start_wtime = MPI::Wtime();
   
    MPIParallelization(current_array, number_rows);

    if (commsize > 1) {
        GatheringAllWork(current_array, main_array, number_rows, rank);
    }

    double end_wtime = MPI::Wtime();
    if (rank == 0) {
        std::cout << "exec_time: " << end_wtime - start_wtime << std::endl;
        DumpToFile("baseline_MYYYY.txt", main_array);
        delete[] main_array;
    }
    
    delete[] current_array;

    MPI::Finalize();
    return 0;
}

void MPIParallelization(double *current_array, int number_rows)
{
    for (int i = 0; i < number_rows; ++i) {
        for (int j = 0; j < JSIZE; ++j) {
            current_array[JSIZE * i + j] = std::sin(2 * current_array[JSIZE * i + j]);
        }
    }
}

void GatheringAllWork(double *current_array, double *main_array, int number_rows, int rank)
{
    int *lengths = nullptr;
    int *put_positions = nullptr;

    if (rank == 0) {
        int commsize = MPI::COMM_WORLD.Get_size();
        lengths = new int[commsize];
        put_positions = new int[commsize];

        int distance = ISIZE / commsize;
        int remainder = ISIZE % commsize;

        for (int i = 0; i < commsize; ++i) {
            int start = distance * i;
            int end = distance * (i + 1);

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
            lengths[i] = (end - start) * JSIZE;
            put_positions[i] = start * JSIZE;
            // std::cout << lengths[i] << " " << put_positions[i] << std::endl;
        }
    }

    MPI_Gatherv(current_array, number_rows * JSIZE, MPI::DOUBLE, main_array, lengths, put_positions, MPI::DOUBLE, 0, MPI_COMM_WORLD);
    
    delete[] lengths;
    delete[] put_positions;
}

void DumpToFile(const char *filename, double *array)
{
    FILE *ff = nullptr;
    ff = fopen(filename, "w");

    for (int i = 0; i < ISIZE; ++i) {
        for (int j = 0; j < JSIZE; ++j) {
            fprintf(ff, "%f", array[JSIZE * i + j]);
        }
        fprintf(ff, "\n");
    }
    fclose(ff);
}

