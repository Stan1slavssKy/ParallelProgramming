#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <cstring>

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

/* Calculate the last factorial required for input accuracy to get the correct exponent */
int GetFactorialNumber(int accuracy) 
{
    mpz_t rop; // required_factorial_accuracy
    mpz_init(rop);

    mpz_t base;
    mpz_init(base);

    mpz_t factorial;
    mpz_init(factorial);
	mpz_set_ui(factorial, 2);

	mpz_set_ui(base, 10);
    mpz_pow_ui(rop, base, accuracy); // rop = 10 ^ accuracy

    int factorial_number = 2;

    while (mpz_cmp(rop, factorial) >= 0) {
        mpz_mul_si(factorial, factorial, ++factorial_number);
    }

    mpz_clear(factorial);
    mpz_clear(rop);
	mpz_clear(base);
    return factorial_number;
}

double CalculatePartExponentSum(int factorial_number, int start, int end)
{
    //std::cout << " [ " << start << ", " << end << " ] " << std::endl;
    mpz_t sum;
    mpz_init(sum);
    mpz_t factorial;
    mpz_init(factorial);
    mpz_fac_ui(factorial, factorial_number);

    for (int idx = 2; idx <= start; ++idx) {
        mpz_cdiv_ui(factorial, idx);
    }
    mpz_add(sum, sum, factorial);

    for (int idx = start + 1; idx < end; ++idx) {
        mpz_cdiv_ui(factorial, idx);
        mpz_add(sum, sum, factorial);
    }

    double part_summ = 

    mpz_clear(factorial);
    mpz_clear(sum);

    return part_summ;
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

    int factorial_number = GetFactorialNumber(accuracy);

    //std::cout << "factorial_number = " << factorial_number << std::endl;

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
    mpf_set_ui(result_sum, 1);

    double part_sum = CalculatePartExponentSum(factorial_number, start, end);

    mpf_clear(result_sum);

    MPI_Finalize();

    return 0;
}