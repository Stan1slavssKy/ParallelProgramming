#include <mpi.h>
#include <cmath>
#include <unistd.h>
#include <string>
#include <algorithm>
#include <cstring>

#include "array_file_output.h"
#include "scoped_time_measure.h"

static constexpr size_t ISIZE = 25;
static constexpr size_t JSIZE = ISIZE;

class Worker;

void GatheringAllWork(const Worker& worker);
void DumpToFile(const char *filename, double *array);

struct Worker {
    // Accepts a range of original values that will need to be calculated
    explicit Worker(int rank, size_t commsize, size_t start, size_t end) : 
        rank_(rank), commsize_(commsize), start_(start), end_(end)
    {
        if (rank_ == 0) {
            main_array_ = new double[ISIZE * JSIZE];
            std::fill_n(main_array_, ISIZE * JSIZE, 0);
        }

        if (static_cast<size_t>(rank_) == commsize_ - 1) {
            last_rank_ = true;
        }

        calc_size_ = end_ - start_;
        calc_array_ = new double[calc_size_ * JSIZE];

        first_array_ = new double[calc_size_ * JSIZE];
        
        for (size_t counter = 0; counter < calc_size_; ++counter) {
            for (size_t j = 0; j < JSIZE; ++j) {
                first_array_[JSIZE * counter + j] = 10 * (start_ + counter) + j;
            }
        }
        if (!last_rank_) {
            // Since i <- i+ 1 we need only one boundary element with idx = end
            // [start, end) => idx = end needed for calculate end - 1
            boundary_arrray_ = new double[1 * JSIZE];

            for (size_t j = 0; j < JSIZE; ++j) {
                boundary_arrray_[j] = 10 * end_ + j;
            }
        }
        // StartMeasusre()
        for (size_t counter = 0; counter < calc_size_; ++counter) {
            for (size_t j = 0; j < JSIZE; ++j) {
                first_array_[JSIZE * counter + j] = std::sin(0.1 * first_array_[JSIZE * counter + j]);
            }
        }

        if (!last_rank_) {
            for (size_t j = 0; j < JSIZE; ++j) {
                boundary_arrray_[j] = std::sin(0.1 * boundary_arrray_[j]);
            }
        }
    }

    ~Worker()
    {
        if (rank_ == 0) {
            delete[] main_array_;
        }
        
        if (!last_rank_) {
            delete[] boundary_arrray_;
        }

        delete[] calc_array_;
    }

    void Calculate()
    {
        for (size_t counter = 0; counter < calc_size_ - 1; ++counter) {
            for (size_t j = 2; j < JSIZE; ++j) {
                calc_array_[JSIZE * counter + j] = 1.5 * first_array_[JSIZE * (counter + 1) + (j - 2)];
            }
        }

        if (!last_rank_) {
            for (size_t j = 2; j < JSIZE; ++j) {
                calc_array_[JSIZE * (calc_size_ - 1) + j] = 1.5 * boundary_arrray_[JSIZE * ((calc_size_ - 1) + 1) + (j - 2)];
            }
        }
    }

    int rank_ {-1};
    size_t commsize_ {0};

    bool last_rank_ {false};

    size_t start_ {0};
    size_t end_ {0};

    size_t calc_size_ {0};
    // Array for calculation part of each worker
    double *calc_array_ {nullptr};
    // Part of array, from which calc array will be calculate
    double *first_array_ {nullptr};
    // Array for elements from different part. Actual for all workers except last
    double *boundary_arrray_ {nullptr};
    // Array with full calculations. Actual only for root (0) worker
    double *main_array_ {nullptr};
};

int main(int argc, char **argv) {
    MPI::Init(argc, argv);

    int rank = MPI::COMM_WORLD.Get_rank();
    size_t commsize = MPI::COMM_WORLD.Get_size();

    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(2 * rank, &mask);
    sched_setaffinity(getpid(), sizeof(cpu_set_t), &mask);

    size_t distance = ISIZE / commsize;
    size_t remainder = ISIZE % commsize;
    size_t start = distance * rank;
    size_t end = distance * (rank + 1);

    // If the remainder is not zero, we add 1 term of the sum to each process starting from zero
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

    printf("rank = %d, [%ld, %ld)\n", rank, start, end);

    Worker worker(rank, commsize, start, end);

    double start_wtime = MPI::Wtime();
    
    worker.Calculate();

    if (commsize > 1) {
        GatheringAllWork(worker);
    }
    
    if (rank == 0) {
        if (commsize == 1) {
    //         for (int i = 0; i < worker.number_init_states_; ++i) {
    //             double *current_array = worker.subarrays_[i]->Get();
    //             size_t transition_nmb = worker.subarrays_[i]->GetTransitionNmb();

    //             int offset = worker.start_ + i;
    //             // std::cout << transition_nmb << std::endl;
    //             for (size_t counter = 0; counter < transition_nmb + 1; ++counter) {
    //                 std::memcpy(worker.main_array_ + offset * ISIZE, current_array + counter * JSIZE, JSIZE * sizeof(double));
    //                 offset += 8;
    //             }
    //         }
        }
        double end_wtime = MPI::Wtime();
        std::cout << "exec_time: " << end_wtime - start_wtime << std::endl;
        DumpToFile(("third_task_" + std::to_string(commsize) + ".txt").c_str(), worker.main_array_);
    }

    MPI::Finalize();
    return 0;
}

void GatheringAllWork(const Worker& worker)
{
    if (worker.rank_ == 0) {
        double *main_array = worker.main_array_;
        size_t offset = 0;
        
        {
            size_t zero_array_size = worker.calc_size_ * JSIZE;
            std::memcpy(main_array, worker.calc_array_, zero_array_size * sizeof(double));
            offset += zero_array_size;
        }
        printf("offset = %ld\n", offset);

        for (int cur_rank = 1; cur_rank < static_cast<int>(worker.commsize_); ++cur_rank) {
            // printf("offset = %ld\n", offset);
            MPI_Status status;

            size_t recv_size = 0;
            printf("offset = %ld\n", offset);
            MPI_Recv(&recv_size, 1, MPI::LONG_LONG, cur_rank, 0, MPI_COMM_WORLD, &status);
            printf("recv_size = %ld, offset = %ld\n", recv_size, offset);

            double *recv_array = new double[recv_size];

            // printf("recv_size = %ld, offset = %ld\n", recv_size, offset);
            MPI_Recv(recv_array, recv_size, MPI::DOUBLE, cur_rank, 0, MPI_COMM_WORLD, &status);    
            std::memcpy(main_array + offset, recv_array, recv_size * sizeof(double));
            offset += recv_size;

            delete[] recv_array;
        }
    }
    else {
        size_t send_size = worker.calc_size_ * JSIZE;
        MPI_Send(&send_size, 1, MPI::LONG_LONG, 0, 0, MPI_COMM_WORLD);
        
        double *current_array = worker.calc_array_;
        MPI_Send(current_array, send_size, MPI::DOUBLE, 0, 0, MPI_COMM_WORLD);
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
