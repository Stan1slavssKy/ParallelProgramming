#include <mpi.h>
#include <cmath>
#include <unistd.h>
#include <string>
#include <algorithm>
#include <cstring>

#include "array_file_output.h"
#include "scoped_time_measure.h"

static constexpr size_t ISIZE = 10000;
static constexpr size_t JSIZE = ISIZE;
static constexpr size_t ORIGINAL_VALUES_NMB = 8;

class Worker;

void GatheringAllWork(const Worker& worker);
void DumpToFile(const char *filename, double *array);

class TransitionArray {
public:
    explicit TransitionArray(size_t array_size, int transitions_nmb) : array_size_(array_size), transition_number_(transitions_nmb)
    {
        array_ = new double[array_size_];
    }

    ~TransitionArray()
    {
        delete[] array_;
    }

    double *Get() const
    {
        return array_;
    }

    size_t GetTransitionNmb() const
    {
        return transition_number_;
    }

private:
    double *array_ {nullptr};
    size_t array_size_ {0};
    size_t transition_number_ {0};
};

struct Worker {
    // Accepts a range of original values that will need to be calculated
    explicit Worker(int rank, int start, int end) : rank_(rank), start_(start), end_(end)
    {
        if (rank_ == 0) {
            main_array_ = new double[JSIZE * ISIZE];
        }
        // [start; end) => number of initial state to be calculated end - start
        number_init_states_ = end - start;
        // For each state we create array, subarrays_ saved pointers to this arrays
        subarrays_ = new TransitionArray*[number_init_states_];

        for (int i = 0; i < number_init_states_; ++i) {
            size_t transition_number = ((ISIZE - 1) - (start + i)) / ORIGINAL_VALUES_NMB;

            auto *transition_array = new TransitionArray(JSIZE * (transition_number + 1), transition_number);
            subarrays_[i] = transition_array;
            
            double *current_array = transition_array->Get();
            // Initialize array
            for (size_t counter = 0; counter < transition_number + 1; ++counter) {
                for (size_t j = 0; j < JSIZE; ++j) {
                    current_array[JSIZE * counter + j] = 10 * (start + i + 8 * counter) + j;
                }
            }
        }
    }
    ~Worker()
    {
        if (rank_ == 0) {
            delete[] main_array_;
        }
        for (int i = 0; i < number_init_states_; ++i) {
            delete subarrays_[i];
        }
        delete[] subarrays_;
    }
    void Calculate()
    {
        for (int arr_idx = 0; arr_idx < number_init_states_; ++arr_idx) {
            double *current_array = subarrays_[arr_idx]->Get();
            size_t transition_nmb = subarrays_[arr_idx]->GetTransitionNmb();

            for (size_t counter = 0; counter < transition_nmb; ++counter) {
                for (size_t j = 0; j < JSIZE - 3; ++j) {
                    current_array[JSIZE * (counter + 1) + j] = std::sin(4 * current_array[JSIZE * counter + (j + 3)]);
                }
            }
        }
    }

    int rank_ {-1};

    int start_ {0};
    int end_ {0};
    int number_init_states_ {0};
    
    TransitionArray **subarrays_ {nullptr};
    // actual only for root rank (0)
    double *main_array_ {nullptr};
};

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

    Worker worker(rank, start, end);

    double start_wtime = MPI::Wtime();
    
    worker.Calculate();

    if (commsize > 1) {
        GatheringAllWork(worker);
    }
    
    if (rank == 0) {
        if (commsize == 1) {
            for (int i = 0; i < worker.number_init_states_; ++i) {
                double *current_array = worker.subarrays_[i]->Get();
                size_t transition_nmb = worker.subarrays_[i]->GetTransitionNmb();

                int offset = worker.start_ + i;
                // std::cout << transition_nmb << std::endl;
                for (size_t counter = 0; counter < transition_nmb + 1; ++counter) {
                    std::memcpy(worker.main_array_ + offset * ISIZE, current_array + counter * JSIZE, JSIZE * sizeof(double));
                    offset += 8;
                }
            }
        }
        double end_wtime = MPI::Wtime();
        std::cout << "exec_time: " << end_wtime - start_wtime << std::endl;
        // DumpToFile(("first_task_" + std::to_string(commsize) + ".txt").c_str(), worker.main_array_);
    }

    MPI::Finalize();
    return 0;
}

void GatheringAllWork(const Worker& worker)
{
    if (worker.rank_ == 0) {
        double *main_array = worker.main_array_; 
        auto commsize = MPI::COMM_WORLD.Get_size();

        for (int i = 0; i < worker.number_init_states_; ++i) {
            double *current_array = worker.subarrays_[i]->Get();
            size_t transition_nmb = worker.subarrays_[i]->GetTransitionNmb();

            int offset = worker.start_ + i;

            for (size_t counter = 0; counter < transition_nmb + 1; ++counter) {
                std::memcpy(main_array + offset * ISIZE, current_array + counter * JSIZE, JSIZE * sizeof(double));
                offset += 8;
            }
        }

        for (int cur_rank = 1; cur_rank < commsize; ++cur_rank) {
            MPI_Status status;
            size_t number_recvs = 0;
            MPI_Recv(&number_recvs, 1, MPI::LONG_LONG, cur_rank, 0, MPI_COMM_WORLD, &status);

            for (size_t i = 0; i < number_recvs; ++i) {
                double *temp_array = new double[JSIZE];
                int offset = 0;

                MPI_Recv(&offset,    1,     MPI::INT,    cur_rank, 0, MPI_COMM_WORLD, &status);
                MPI_Recv(temp_array, JSIZE, MPI::DOUBLE, cur_rank, 0, MPI_COMM_WORLD, &status);
                
                std::memcpy(main_array + offset * JSIZE, temp_array, JSIZE * sizeof(double));
                delete[] temp_array;
            }
        }
    }
    else {
        size_t init_states_nmb = worker.number_init_states_;
        size_t number_of_sends = 0;

        for (size_t i = 0; i < init_states_nmb; ++i) {
            size_t transition_nmb = worker.subarrays_[i]->GetTransitionNmb();
            number_of_sends += (transition_nmb + 1);
        }
        
        MPI_Send(&number_of_sends, 1, MPI::LONG_LONG, 0, 0, MPI_COMM_WORLD);

        for (size_t i = 0; i < init_states_nmb; ++i) {
            double *current_array = worker.subarrays_[i]->Get();
            size_t transition_nmb = worker.subarrays_[i]->GetTransitionNmb();

            int offset = worker.start_ + i;

            for (size_t counter = 0; counter < transition_nmb + 1; ++counter) {
                MPI_Send(&offset                        , 1    , MPI::INT     , 0, 0, MPI_COMM_WORLD);
                MPI_Send(current_array + counter * JSIZE, JSIZE, MPI::DOUBLE  , 0, 0, MPI_COMM_WORLD);
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
