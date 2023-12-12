#include <cmath>

#include "array_file_output.h"
#include "scoped_time_measure.h"

static constexpr size_t ISIZE = 15000;
static constexpr size_t JSIZE = ISIZE;

void InitVariantBaseline();
void InitVariantFirstTask();
void InitVariantSecondTask();
void InitVariantThirdTask();

int main() {
    // InitVariantBaseline();
    // InitVariantFirstTask();
    // InitVariantSecondTask();
    InitVariantThirdTask();
    return 0;
}

void InitVariantBaseline()
{
    ArrayHandler array_handler(ISIZE, JSIZE);
    array_handler.DefaultFillIn();

    double *array = array_handler.GetArray();
    {
        ScopedTimeMeasure stm;

        for (size_t i = 0; i < ISIZE; ++i) {
            for (size_t j = 0; j < JSIZE; ++j) {
                array[JSIZE * i + j] = std::sin(2 * array[JSIZE * i + j]);
            }
        }
    }

    array_handler.DumpToFile("file_to_compare_baseline.txt");
}

void InitVariantFirstTask()
{
    ArrayHandler array_handler(ISIZE, JSIZE);
    array_handler.DefaultFillIn();

    double *array = array_handler.GetArray();
    {
        ScopedTimeMeasure stm;

        for (size_t i = 8; i < ISIZE; ++i) {
            for (size_t j = 0; j < JSIZE - 3; ++j) {
                array[JSIZE * i + j] = std::sin(4 * array[JSIZE * (i - 8) + (j + 3)]);
            }
        }
    }

    array_handler.DumpToFile("mini_file_to_compare_1_task.txt");
}

void InitVariantSecondTask()
{
    ArrayHandler array_handler(ISIZE, JSIZE);
    array_handler.DefaultFillIn();

    double *array = array_handler.GetArray();
    {
        ScopedTimeMeasure stm;

        for (size_t i = 0; i < ISIZE - 4; ++i) {
            for (size_t j = 0; j < JSIZE - 2; ++j) {
                array[JSIZE * i + j] = std::sin(0.1 * array[JSIZE * (i + 4) + (j + 2)]);
            }
        }
    }

    array_handler.DumpToFile("file_to_compare_2_task.txt");  
}

void InitVariantThirdTask()
{
    ArrayHandler array_handler_a(ISIZE, JSIZE);
    ArrayHandler array_handler_b(ISIZE, JSIZE);

    array_handler_a.DefaultFillIn();
    array_handler_b.ZeroFillIn();

    double *a = array_handler_a.GetArray();
    double *b = array_handler_b.GetArray();

    {
        ScopedTimeMeasure stm;

        for (size_t i = 0; i < ISIZE; ++i) {
            for (size_t j = 0; j < JSIZE; ++j) {
                a[JSIZE * i + j] = std::sin(0.1 * a[JSIZE * i + j]);
            }
        }
        for (size_t i = 0; i < ISIZE - 1; ++i) {
            for (size_t j = 2; j < JSIZE; ++j) {
                b[JSIZE * i + j] = 1.5 * a[JSIZE * (i + 1) + (j - 2)];
            }
        }
    }

    array_handler_b.DumpToFile("mini_file_to_compare_3_task.txt");  
}
