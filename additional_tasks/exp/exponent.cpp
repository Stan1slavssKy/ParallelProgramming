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

    while (mpz_cmp(rop, factorial) >= 0) { // find n that n! > 10 ^ accuracy
        mpz_mul_si(factorial, factorial, ++factorial_number);
    }

    mpz_clear(factorial);
    mpz_clear(rop);
	mpz_clear(base);
    return factorial_number;
}

void CalculatePartExponentSum(mpf_t part_sum, int factorial_number, int start, int end)
{
   // std::cout << " [ " << start << ", " << end << " ] " << std::endl;
    mpz_t sum;
    mpz_init(sum);
    if ((end - start) != 1) {
	    mpz_set_ui(sum, 1); // 1 is the last element in summ
    }

    int current_num = end - 1;

    mpz_t reverse_factorial;
    mpz_init(reverse_factorial);
	mpz_set_ui(reverse_factorial, end - 1);
    mpz_add(sum, sum, reverse_factorial);
    --current_num;

    for (int idx = 0; idx < end - start - 2; ++idx) {
        mpz_mul_ui(reverse_factorial, reverse_factorial, current_num);
        mpz_add(sum, sum, reverse_factorial);
        // gmp_printf ("current_sum = %Zd\n", sum);
        --current_num;
        if (current_num == 0) {
            break;
        }
    }

    while (current_num > 1) {
        mpz_mul_ui(reverse_factorial, reverse_factorial, current_num);
        --current_num;
    }

    mpf_t sum_f;
    mpf_init(sum_f);
    mpf_set_z(sum_f, sum);

    mpf_t reverse_factorial_f;
    mpf_init(reverse_factorial_f);
    mpf_set_z(reverse_factorial_f, reverse_factorial);

    // gmp_printf ("sum_f____[D] = %.*Ff\n", 2000, sum_f);
    // gmp_printf ("reverse_f[D] = %.*Ff\n", 2000, reverse_factorial_f);

    mpf_div(part_sum, sum_f, reverse_factorial_f);
    
    mpz_clear(reverse_factorial);
    mpz_clear(sum);
    mpf_clear(reverse_factorial_f);
    mpf_clear(sum_f);
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

    mpf_set_default_prec(64 + ceil(3.33 * accuracy));

    int factorial_number = GetFactorialNumber(accuracy);

 //    std::cout << "factorial_number = " << factorial_number << std::endl;

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

    CalculatePartExponentSum(part_sum, factorial_number, start, end);
    // gmp_printf ("[Rank = %d] part_suma[D] = %.*Ff\n", task_id, 2000, part_sum);

    MPI_Status status;

    mp_exp_t exp = 0;
    char* part_sum_str = mpf_get_str(nullptr, &exp, 0, 0, part_sum);
    size_t part_sum_str_len = strlen(part_sum_str);
    // printf ("[Rank = %d] exp = %d, part_suma[D] = %s\n", task_id, exp, part_sum_str);
    if (task_id == root_task_id) {
        mpf_add(result_sum, result_sum, part_sum);
        //gmp_printf ("part_summmma( = %.*Ff\n", 2000, part_sum);

        for (int idx = 1; idx < commsize; ++idx) { // idx == task_id
            MPI_Recv(&exp, 1, MPI_LONG, idx, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&part_sum_str_len, 1, MPI_LONG, idx, 0, MPI_COMM_WORLD, &status);
            
            if (exp <= 0) {                                             // | 0 | . | exp | length |
                size_t recv_str_len = part_sum_str_len + abs(exp) + 2;  // |___|___|_____|________|
                                                                        //  |             |
                part_sum_str = new char[recv_str_len];                  //  part_s_str    receive here
                std::fill_n(part_sum_str, recv_str_len, '0');
                std::fill(part_sum_str + 1, part_sum_str + 2, '.');

                MPI_Recv(part_sum_str + abs(exp) + 2, part_sum_str_len, MPI_CHAR, idx, 0, MPI_COMM_WORLD, &status);
            } else {
                part_sum_str = new char[part_sum_str_len];
                MPI_Recv(part_sum_str, part_sum_str_len, MPI_CHAR, idx, 0, MPI_COMM_WORLD, &status); 
            }

            mpf_set_str(part_sum, part_sum_str, 10);
                    // std::cout << part_sum_str << "[TASK_ID = " << idx << "]" << "exp =" << exp << std::endl;
            delete[] part_sum_str;
            // gmp_printf ("part_summmma = %.*Ff\n", 2000, part_sum);
            mpf_add(result_sum, result_sum, part_sum);
        }      
        
        gmp_printf ("[RESULT] Exp = %.*Ff\n", accuracy, result_sum);
    } else {
        MPI_Send(&exp, 1, MPI_LONG, root_task_id, 0, MPI_COMM_WORLD);
        MPI_Send(&part_sum_str_len, 1, MPI_LONG, root_task_id, 0, MPI_COMM_WORLD);
        MPI_Send(part_sum_str, part_sum_str_len, MPI_CHAR, root_task_id, 0, MPI_COMM_WORLD);
    }

    mpf_clear(result_sum);
    mpf_clear(part_sum);

    MPI_Finalize();

    return 0;
}