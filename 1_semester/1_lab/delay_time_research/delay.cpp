#include <iostream>
#include <cmath>
#include <cstddef>
#include <mpi.h>

int main(int argc, char** argv)
{
    MPI::Init(argc, argv);

    if (argc != 2) {
        std::cout << "Please pass arg with iterations number" << std::endl;
        MPI::Finalize();
        return 0;
    }

    int commsize = MPI::COMM_WORLD.Get_size();
    if (commsize < 2) {
        std::cout << "Need more than 1 process to test send/recv time" << std::endl;
        MPI::Finalize();
        return 0;
    }

    int rank = MPI::COMM_WORLD.Get_rank();

    size_t iter_nmb = std::atoi(argv[1]);
    double send_nmb = 0;

    if (rank == 0) {
        double start_time = MPI::Wtime();
        for (size_t i = 0; i < iter_nmb; ++i) {
            MPI::COMM_WORLD.Recv(&send_nmb, 1, MPI::DOUBLE, 1, 0);
        }
        double finish_time = MPI::Wtime();

        std::cout << "Single transfer time: " << (finish_time - start_time) / iter_nmb * 1e6 << " us" << std::endl;
    } else if (rank == 1) {
        for (size_t i = 0; i < iter_nmb; ++i) {
            MPI::COMM_WORLD.Send(&send_nmb, 1, MPI::DOUBLE, 0, 0);
        }
    }

    MPI::Finalize();
    return 0;
}
