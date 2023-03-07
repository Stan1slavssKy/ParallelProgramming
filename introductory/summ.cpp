#include <mpi.h>
#include <iostream>
#include <cstdlib>

int GetSumElementsNumber(int argc, char** argv)
{
    if (argc != 2) {
        std::cout << "Error. Please pass the parameter at startup" << std::endl;
        return -1;
    }

    int N = std::atoi(argv[argc - 1]);
    if (N == 0) {
        std::cout << "You pass into arg N = 0 or atoi conversion can't be performed" << std::endl;
        return 0;
    }

    return N;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    
    int N = GetSumElementsNumber(argc, argv);

    const int root_task_id = 0;

    double result_sum = 0;
    double part_sum = 0;
    int commsize = 0;
    int task_id = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);

    int distance = N / commsize;
    int remainder = N % commsize;
    int start = 1 + distance * task_id;
    int end = 1 + distance * (task_id + 1);

    // If the remainder is not zero, we add 1 term of the sum to each process starting from zero
    if (remainder != 0) {
        if (remainder - task_id > 0) {
            start += task_id;
            end += task_id + 1;
        }
        else {
            start += remainder;
            end += remainder;
        }
    }

    for (int i = start; i < end; ++i) {
        part_sum += 1.0 / i;
    }

    MPI_Reduce(&part_sum, &result_sum, 1, MPI_DOUBLE, MPI_SUM, root_task_id, MPI_COMM_WORLD);

    if (task_id == root_task_id) {
        std::cout << "Result sum = " << result_sum << std::endl;
    }

    MPI_Finalize();

    return 0;
}