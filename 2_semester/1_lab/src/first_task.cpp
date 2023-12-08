#include <mpi.h>
#include <cmath>
#include <unistd.h>
#include <string>
#include <algorithm>
#include <cstring>

#include "array_file_output.h"
#include "scoped_time_measure.h"

static constexpr size_t ISIZE = 150;
static constexpr size_t JSIZE = ISIZE;

void MPIParallelization(double **subarrays, int nmb_orig_to_calc, int nmb_calcs, int rank);
void GatheringAllWork(double **subarrays, double *main_array, int nmb_orig_to_calc, int rank, int start, int i_elems_per_array, int *length_array, int nmb_calc);
void DumpToFile(const char *filename, double *array);

int main(int argc, char **argv) {
    MPI::Init(argc, argv);

    auto rank = MPI::COMM_WORLD.Get_rank();
    auto commsize = MPI::COMM_WORLD.Get_size();
    if (commsize > 8) {
        MPI::Finalize();
        std::cerr << "Commsize should be <= 8" << std::endl;
        return 0;
    }

    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(2 * rank, &mask);
    sched_setaffinity(getpid(), sizeof(cpu_set_t), &mask);
    
    double *main_array = nullptr;
    if (rank == 0) {
        main_array = new double[ISIZE * JSIZE];
    }

    int number_original_values = 8;

    int distance = number_original_values / commsize;
    int remainder = number_original_values % commsize;
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

    int nmb_orig_to_calc = end - start;
    int nmb_calc = (ISIZE - rank) / number_original_values;
    int i_elems_per_array = JSIZE * nmb_calc;

    // std::cout << "[" << start << " " << end << ")" << std::endl;
    // std::cout << "nmb_elems = " << i_elems_per_array << std::endl;

    double **subarrays = new double*[nmb_orig_to_calc];

    // for (int i = 0; i < nmb_orig_to_calc; ++i) {
    //     double *current_array = new double[i_elems_per_array + 1];
    //     subarrays[i] = current_array;

    //     for (size_t j = 0; j < JSIZE; ++j) {
    //         current_array[j] = 10 * (i + start) + j;
    //     }
    // }

    for (int i = 0; i < nmb_orig_to_calc; ++i) {
        double *current_array = new double[i_elems_per_array + 1];
        subarrays[i] = current_array;

        for (int counter = 0; counter < nmb_calc; ++counter) {
            for (size_t j = 0; j < JSIZE; ++j) {
                current_array[JSIZE * counter + j] = 10 * (i + start + 8 * counter) + j;
            }
        }
    }

    int *length_array = nullptr;

    if (rank == 0) {
        length_array = new int[number_original_values];
        for (int i = 0; i < number_original_values; ++i) {
            length_array[i] = (ISIZE - i) / number_original_values;
            std::cout << "============================" << length_array[i] << std::endl;
        }
    }

    double start_wtime = MPI::Wtime();
   
    MPIParallelization(subarrays, nmb_orig_to_calc, nmb_calc, rank);

    if (commsize > 1) {
        GatheringAllWork(subarrays, main_array, nmb_orig_to_calc, rank, start, i_elems_per_array, length_array, nmb_calc);
    }

    double end_wtime = MPI::Wtime();
    if (rank == 0) {
        std::cout << "exec_time: " << end_wtime - start_wtime << std::endl;
        // if (commsize > 1) {
        DumpToFile(("first_task_" + std::to_string(commsize) + ".txt").c_str(), main_array);
        // } else {
        //     DumpToFile(("baseline_" + std::to_string(commsize) + ".txt").c_str(), current_array);
        // }
        delete[] main_array;
    }

    for (int i = 0; i < nmb_orig_to_calc; ++i) {
        delete[] subarrays[i];
    }

    delete[] subarrays;
    delete[] length_array;

    MPI::Finalize();
    return 0;
}

void MPIParallelization(double **subarrays, int nmb_orig_to_calc, int nmb_calcs, int rank)
{
    for (int arr_idx = 0; arr_idx < nmb_orig_to_calc; ++arr_idx) {
        double *current_array = subarrays[arr_idx];

        for (int counter = 0; counter < nmb_calcs - 1; ++counter) {
            for (size_t j = 0; j < JSIZE - 3; ++j) {
                current_array[JSIZE * (counter + 1) + j] = std::sin(4 * current_array[JSIZE * counter + (j + 3)]);
            }
        }
        // if (rank == 0) {
        //     for (size_t j = 0; j < JSIZE; ++j) {
        //         std::cout << std::sin(4 * current_array[JSIZE * 0 + (j + 3)]) << std::endl;
        //     }
        //     DumpToFile("debug.txt", current_array);
        // }
    }
}

void GatheringAllWork(double **subarrays, double *main_array, int nmb_orig_to_calc, int rank, int start, int i_elems_per_array, int *length_array, int nmb_calc)
{
    if (rank == 0) {
        auto commsize = MPI::COMM_WORLD.Get_size();

        for (int i = 0; i < nmb_orig_to_calc; ++i) {
            double *current_array = subarrays[i];
            int offset = start;

            for (size_t counter = 0; counter < nmb_calc; ++counter) {
                std::memcpy(main_array + offset * ISIZE, current_array + counter * JSIZE, JSIZE * sizeof(double));
                offset += 8;
            }
        }

        for (int cur_rank = 1; cur_rank < commsize; ++cur_rank) {
            for (int i = 0; i < length_array[cur_rank]; ++i) {
                // std::cout << "{============== wait for " << cur_rank << std::endl;

                double *temp_array = new double[JSIZE];
                int offset = 0;
                MPI_Status status;

                MPI_Recv(&offset,    1,                 MPI::INT,    cur_rank, 0, MPI_COMM_WORLD, &status);
                MPI_Recv(temp_array, JSIZE, MPI::DOUBLE, cur_rank, 0, MPI_COMM_WORLD, &status);
                
                // std::cout << offset << std::endl;

                // int send_rank = status.MPI_SOURCE;
                // std::cout << "{============== recv " << offset << " from " << send_rank << std::endl;
                
                std::memcpy(main_array + offset * ISIZE, temp_array, JSIZE * sizeof(double));
                delete[] temp_array;
            }
        }
    }
    else {
        for (int i = 0; i < nmb_orig_to_calc; ++i) {
            double *current_array = subarrays[i];
            
            // std::cout << "===================================================" << std::endl;
            int offset = start;

            for (int counter = 0; counter < nmb_calc; ++counter) {
                // std::cout << "send from " << rank << std::endl;
                // int offset = start + 8 * counter;
                // std::cout << offset << std::endl;
                MPI_Send(&offset,                         1,     MPI::INT,    0, 0, MPI_COMM_WORLD);
                MPI_Send(current_array + counter * JSIZE, JSIZE, MPI::DOUBLE, 0, 0, MPI_COMM_WORLD);
                offset += 8;
            }
        }
    }
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
