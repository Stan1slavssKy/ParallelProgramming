#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <cstring>
#include <algorithm>
#include <cmath>

#include <mpi.h>
#include <gmp.h>

enum Errors {
    FEW_ARGUMENTS = -1,
    MANY_ARGUMENTS = -2,
    STRTOL_ERROR = -3,
    INVALID_NUMBER_INPUT = -4,
};

int GetAccuracy(int argc, char** argv)
{
    int saved_errno = errno;

    if (argc > 2) {
        return MANY_ARGUMENTS;
    }
    if (argc < 2) {
        return FEW_ARGUMENTS;
    }

    int accuracy = std::strtol(argv[argc - 1], nullptr, 10);
    if (errno == ERANGE) {
        return STRTOL_ERROR;
    }
    if (accuracy <= 0) {
        return INVALID_NUMBER_INPUT;
    }

    errno = saved_errno;
    return accuracy;
}

void AccuracyErrorHandling(int accuracy, int errnum)
{
    switch (accuracy)
    {
    case FEW_ARGUMENTS:
        std::cerr << "Too few arguments in input!" << std::endl;
        break;
    case MANY_ARGUMENTS:
        std::cerr << "Too many arguments in input!" << std::endl;
        break;
    case STRTOL_ERROR:
        std::cerr << "Converted value falls out in std::strtol!" << std::strerror(errnum) << std::endl;
        break;
    case INVALID_NUMBER_INPUT:
        std::cerr << "Invalid number input!" << std::endl;
        break;
    default:
        std::cerr << "Unknown error in GetAccuracy!" << std::endl;
        break;
    }
}

/**
 * Newton's method for solving a nonlinear equation
 * We need n! > 10 ^ accuracy
 * n * ln(n) - n > accuracy * ln(10)
 * x * ln(x) - x > cnst
 * F(x_n) = x_n * ln(x_n) - x_n - cnst
 * x_{n+1} = x_n - F(x_n) / F'(x_n)
*/
int GetFactorialNumberNewtonMethod(int accuracy) 
{
    double x_0 = 2.0;
    double cnst = (double)accuracy * std::log(10.0);
    double epsilon = 1.0;

    double prev_x = x_0;
    double cur_x = 0;

    while (true) {
        cur_x = prev_x - ((prev_x * std::log(prev_x) - prev_x - cnst) / std::log(prev_x));
        if (abs(cur_x - prev_x) < epsilon) {
            return (int)ceil(cur_x);
        }
        prev_x = cur_x;
    }

    return -1;
}

void CalculatePartExponentSum(mpf_t part_sum, int factorial_number, int start, int end, int task_id, int commsize)
{
    mpz_t sum;
    mpz_init(sum);
	mpz_set_ui(sum, 1); // 1 is the last element in summ

    int current_num = end - 1;

    mpz_t reverse_factorial;
    mpz_init(reverse_factorial);
	mpz_set_ui(reverse_factorial, end - 1);
    if ((end - start) != 1) {
        mpz_add(sum, sum, reverse_factorial);
        --current_num;
    }

    for (int idx = 0; idx < end - start - 2; ++idx) {
        mpz_mul_ui(reverse_factorial, reverse_factorial, current_num);
        mpz_add(sum, sum, reverse_factorial);
        --current_num;
        if (current_num == 0) {
            break;
        }
    }

    if (commsize > 1) {
        if (task_id == 0) {
            char* factorial_str = mpz_get_str(nullptr, 10, reverse_factorial);
            MPI_Send(factorial_str, strlen(factorial_str) + 1, MPI_CHAR, task_id + 1, 0, MPI_COMM_WORLD);
        }
        if (task_id > 0) {
            MPI_Status status;

            int rcv_fact_len = 0;
            MPI_Probe(task_id - 1, 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_CHAR, &rcv_fact_len);

            char* rcv_factorial_str = (char*)calloc(rcv_fact_len + 1, sizeof(char));
            MPI_Recv(rcv_factorial_str, rcv_fact_len + 1, MPI_CHAR, task_id - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            mpz_t prev_factorial;
            mpz_init_set_str(prev_factorial, rcv_factorial_str, 10);

            mpz_mul_ui(reverse_factorial, reverse_factorial, start);
            mpz_mul(reverse_factorial, reverse_factorial, prev_factorial);

            if (task_id != (commsize - 1)) {
                char* factorial_str = mpz_get_str(nullptr, 10, reverse_factorial);
                int fact_len = std::strlen(factorial_str);
                MPI_Send(factorial_str, fact_len + 1, MPI_CHAR, task_id + 1, 0, MPI_COMM_WORLD);
            }

            free(rcv_factorial_str);
            mpz_clear(prev_factorial);
        }
    }

    mpf_t sum_f;
    mpf_init(sum_f);
    mpf_set_z(sum_f, sum);

    mpf_t reverse_factorial_f;
    mpf_init(reverse_factorial_f);
    mpf_set_z(reverse_factorial_f, reverse_factorial);

    mpf_div(part_sum, sum_f, reverse_factorial_f);
    
    mpz_clears(reverse_factorial, sum, nullptr);
    mpf_clears(reverse_factorial_f, sum_f, nullptr);
}

int main(int argc, char** argv)
{
    int mpi_status = MPI_Init(&argc, &argv);
    if (mpi_status != MPI_SUCCESS) {
        std::cerr << "MPI_Init failed, return value = " << mpi_status << std::endl;
        return EXIT_FAILURE;
    }

    int accuracy = GetAccuracy(argc, argv);
    if (accuracy < 0 || errno != 0) {
        AccuracyErrorHandling(accuracy, errno);
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    mpf_set_default_prec(64 + 4 * accuracy);

    int factorial_number = GetFactorialNumberNewtonMethod(accuracy);

    const int root_task_id = 0;
    int commsize = 0;
    int task_id = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);

    int distance = factorial_number / commsize;
    int remainder = factorial_number % commsize;
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

    mpf_t result_sum;
    mpf_init(result_sum);
    mpf_set_ui(result_sum, 1); // 0!

    mpf_t part_sum;
    mpf_init(part_sum);

    CalculatePartExponentSum(part_sum, factorial_number, start, end, task_id, commsize);

    MPI_Status status;

    if (task_id == root_task_id) {
        mpf_add(result_sum, result_sum, part_sum);

        char *rvc_buf = (char*)calloc(accuracy + 8, sizeof(char));

        for (int idx = 1; idx < commsize; ++idx) {
            MPI_Recv(rvc_buf, accuracy + 8, MPI_CHAR, idx, 0, MPI_COMM_WORLD, &status);
            mpf_set_str(part_sum, rvc_buf, 10);
            mpf_add(result_sum, result_sum, part_sum);
        }
        free(rvc_buf);
        gmp_printf ("[RESULT] Exp = %.*Ff\n", accuracy, result_sum);
    }
    else {
        char *mpi_buf = (char*)calloc(accuracy + 8, sizeof(char));
        char *format_str = (char*)calloc(accuracy + 12, sizeof(char));

        size_t format_str_size = accuracy + 12;

        snprintf(format_str, format_str_size, "%%.%dFf", accuracy);
        gmp_snprintf(mpi_buf, accuracy + 8, format_str, part_sum);
        MPI_Send(mpi_buf, accuracy + 8, MPI_CHAR, 0, 0, MPI_COMM_WORLD);

        free(format_str);
        free(mpi_buf);
    }

    mpf_clear(result_sum);
    mpf_clear(part_sum);

    MPI_Finalize();

    return 0;
}