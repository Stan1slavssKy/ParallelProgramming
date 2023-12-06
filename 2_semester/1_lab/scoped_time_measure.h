#ifndef PARALLELPROGRAMMING_2_SEMESTER_1_LAB_SCOPED_TIME_MEASURE_H
#define PARALLELPROGRAMMING_2_SEMESTER_1_LAB_SCOPED_TIME_MEASURE_H

#include <iostream>
#include <chrono>

class ScopedTimeMeasure {
public:
    using time_point_t = std::chrono::time_point<std::chrono::system_clock>;

    explicit ScopedTimeMeasure()
    {
        start_ = std::chrono::system_clock::now();
    }

    ~ScopedTimeMeasure()
    {
        end_ = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_ - start_;
        std::cout << "[Elapsed time: " << elapsed_seconds.count() << "s]" << std::endl;
    }

private:
    time_point_t start_;
    time_point_t end_;
};

#endif // PARALLELPROGRAMMING_2_SEMESTER_1_LAB_SCOPED_TIME_MEASURE_H
