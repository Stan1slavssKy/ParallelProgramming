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
        
        for (int i = 0; i < ISIZE; ++i) {
            for (int j = 0; j < JSIZE; ++j) {
                array[ISIZE * i + j] = std::sin(2 * array[ISIZE * i + j]);
            }
        }
    }

    array_handler.DumpToFile("baseline.txt");
    return 0;
}
