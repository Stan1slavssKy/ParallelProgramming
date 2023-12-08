#include <cmath>

#include "array_file_output.h"
#include "scoped_time_measure.h"

static constexpr size_t ISIZE = 150;
static constexpr size_t JSIZE = ISIZE;

void InitVariantBaseline();
void InitVariantFirstTask();
void InitVariantSecondTask();
void InitVariantThirdTask();

int main() {
    // InitVariantBaseline();
    InitVariantFirstTask();
    // InitVariantSecondTask();
    // InitVariantThirdTask();
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
