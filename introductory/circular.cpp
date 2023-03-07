#include <mpi.h>
#include <iostream>
#include <cstdlib>

int GetInputNumber(int argc, char** argv)
{
    if (argc != 2) {
        std::cout << "Error. Please pass the parameter at startup" << std::endl;
        return -1;
    }
    return std::atoi(argv[argc - 1]);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    
    int num = GetInputNumber(argc, argv);
    
    int result = num;
    int commsize = 0;
    int task_id = 0;

    const int root_task_id = 0;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);

    if (task_id != root_task_id) {
        MPI_Recv(&result, 1, MPI_INT, task_id - 1, 0, MPI_COMM_WORLD, &status);
    }

    ++result;

    std::cout << "Send from " << task_id << " to " << (task_id + 1) % commsize << std::endl;
    MPI_Send(&result, 1, MPI_INT, (task_id + 1) % commsize, 0, MPI_COMM_WORLD);

    if (task_id == root_task_id) {
        MPI_Recv(&result, 1, MPI_INT, commsize - 1, 0, MPI_COMM_WORLD, &status);
        std::cout << "Begin number = " << num << " Result number = " << result << std::endl;
    }

    MPI_Finalize();

    return 0;
}