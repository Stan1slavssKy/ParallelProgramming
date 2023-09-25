#include <omp.h>

#include <iostream>

int GetInputNumber(int argc, char** argv)
{
    if (argc != 2) {
        std::cout << "Error. Please pass num elemts in sum" << std::endl;
        return -1;
    }
    return std::atoi(argv[argc - 1]);
}

int main(int argc, char** argv)
{
    size_t num = GetInputNumber(argc, argv);
    float result = 0;

#pragma omp parallel for reduction(+:result)
    for (size_t i = 1; i <= num; ++i) {
        result = result + 1.0 / i;
    }

    std::cout << "Result = " << result << std::endl;

    return 0;
}
