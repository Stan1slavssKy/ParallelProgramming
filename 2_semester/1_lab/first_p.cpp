#include <cmath>

#include "array_file_output.h"
#include "scoped_time_measure.h"

static constexpr size_t ISIZE = 5000;
static constexpr size_t JSIZE = 5000;

int main() {
    ArrayHandler array_handler(ISIZE, JSIZE);
    array_handler.DefaultFillIn();

    double *array = array_handler.GetArray();
    {
        ScopedTimeMeasure stm;
        
        for (int i = 8; i < ISIZE; ++i) {
            for (int j = 0; j < JSIZE; ++j) {
                array[ISIZE * i + j] = std::sin(4 * array[ISIZE * (i - 8) + j + 3]);
            }
        }
    }

    array_handler.DumpToFile("first_task_result.txt");
    return 0;
}
