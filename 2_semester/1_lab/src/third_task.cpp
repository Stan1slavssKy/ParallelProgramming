#include <mpi.h>
#include <cmath>
#include <unistd.h>
#include <string>
#include <algorithm>
#include <cstring>

#include "array_file_output.h"
#include "scoped_time_measure.h"

static constexpr size_t ISIZE = 15000;
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

        first_array_ = new double[(calc_size_ + !last_rank_ - 1) * JSIZE];
        
        for (size_t counter = 0; counter < calc_size_ + (!last_rank_) - 1; ++counter) {
            for (size_t j = 0; j < JSIZE; ++j) {
                first_array_[JSIZE * counter + j] = 10 * (start_ + counter + 1) + j;
            }
        }
        
        StartMeasure();

        for (size_t counter = 0; counter < calc_size_ + (!last_rank_) - 1; ++counter) {
            for (size_t j = 0; j < JSIZE; ++j) {
                first_array_[JSIZE * counter + j] = std::sin(0.1 * first_array_[JSIZE * counter + j]);
            }
        }
    }

    ~Worker()
    {
        if (rank_ == 0) {
            delete[] main_array_;
        }
        
        delete[] calc_array_;
    }

    void Calculate()
    {
        for (size_t counter = 0; counter < calc_size_ + (!last_rank_) - 1; ++counter) {
            for (size_t j = 2; j < JSIZE; ++j) {
                calc_array_[JSIZE * counter + j] = 1.5 * first_array_[JSIZE * counter + (j - 2)];
            }
        }
    }

    void StartMeasure()
    {
        if (rank_ == 0) {
            time_start_ = MPI::Wtime();
        }
    }
    
    void EndMeasure()
    {
        if (rank_ == 0) {
            time_stop_ = MPI::Wtime();
        }
    }

    double GetElapsedTime() const
    {
        return time_stop_ - time_start_;
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
    // Array with full calculations. Actual only for root (0) worker
    double *main_array_ {nullptr};

    double time_start_ {0.0};
    double time_stop_  {0.0}; 
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

    Worker worker(rank, commsize, start, end);
    
    worker.Calculate();

    if (commsize > 1) {
        GatheringAllWork(worker);
    }
    
    if (rank == 0) {
        if (commsize == 1) {
            std::memcpy(worker.main_array_ , worker.calc_array_, worker.calc_size_ * JSIZE * sizeof(double));
        }
        worker.EndMeasure();
        double elapsed_time = worker.GetElapsedTime();

        std::cout << "exec_time: " << elapsed_time << std::endl;
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

        for (int cur_rank = 1; cur_rank < static_cast<int>(worker.commsize_); ++cur_rank) {
            MPI_Status status;

            size_t recv_size = 0;
            MPI_Recv(&recv_size, 1, MPI::LONG_LONG, cur_rank, 0, MPI_COMM_WORLD, &status);

            double *recv_array = new double[recv_size];

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
