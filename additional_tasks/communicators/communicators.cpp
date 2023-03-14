#include <mpi.h>
#include <iostream>
#include <cstdlib>

int GetSumElementsNumber(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "Error. Incorrect arguments" << std::endl;
        return -1;
    }

    int N = std::atoi(argv[argc - 1]);
    if (N == 0) {
        std::cerr << "You pass into arg N = 0 or atoi conversion can't be performed" << std::endl;
        return 0;
    }

    return N;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    
    const int root_rank = 0;
    int commsize = 0;
    int rank = 0;
    int N = GetSumElementsNumber(argc, argv);

    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Split communicator*/
    MPI_Comm new_comm;

    int color = -1;
    int key = -1;
    const int first_color = 1;
    const int second_color = 2;

    if (rank == root_rank) {
        color = first_color;
        key = rank;
    }
    else {
        color = second_color;
        key = rank - 1;
    }

    MPI_Comm_split(MPI_COMM_WORLD, color, key, &new_comm);

    /* Get rank and commsize in new communicator and print it */
    int new_rank = 0;
    int new_commsize = 0;

    MPI_Comm_rank(new_comm, &new_rank);
    MPI_Comm_size(new_comm, &new_commsize);

    std::cout << "[Base communicator] " << "rank = " << rank << " commsize = " << commsize << std::endl;
    std::cout << "[New  communicator] " << "new_rank = " << new_rank << " new_commsize = " << new_commsize << std::endl;
    std::cout << std::endl;

    /* Calculate summ in new communicator*/
    if (rank != root_rank) {
        double result_sum = 0;
        double part_sum = 0;

        int distance = N / new_commsize;
        int remainder = N % new_commsize;
        int start = 1 + distance * new_rank;
        int end = 1 + distance * (new_rank + 1);

        if (remainder != 0) { // If the remainder is not zero, we add 1 term of the sum to each process starting from zero
            if (remainder - new_rank > 0) {
                start += new_rank;
                end += new_rank + 1;
            }
            else {
                start += remainder;
                end += remainder;
            }
        }

        for (int i = start; i < end; ++i) {
            part_sum += 1.0 / i;
        }

        MPI_Reduce(&part_sum, &result_sum, 1, MPI_DOUBLE, MPI_SUM, root_rank, new_comm);

        if (new_rank == root_rank) {
            std::cout << "[Result summ = " << result_sum << "]" << std::endl;
        }
    }

    MPI_Comm_free(&new_comm);
    MPI_Finalize();

    return 0;
}