#ifndef PARALLELPROGRAMMING_2_SEMESTER_1_LAB_ARRAY_HANDLER_H
#define PARALLELPROGRAMMING_2_SEMESTER_1_LAB_ARRAY_HANDLER_H

#include <cstddef>
#include <memory>

class ArrayHandler {
public:
    explicit ArrayHandler(size_t x, size_t y) : x_(x), y_(y)
    {
        array_ = std::make_unique<double[]>(x_ * y_);
    }

    ~ArrayHandler() = default;

    void DefaultFillIn()
    {
        for (size_t i = 0; i < x_; ++i) {
            for (size_t j = 0; j < y_; ++j) {
                array_[x_ * i + j] = 10 * i + j;
            }
        }
    }

    void DumpToFile(const char *filename) const
    {
        FILE *ff = nullptr;
        ff = fopen(filename, "w");

        for (size_t i = 0; i < x_; ++i) {
            for (size_t j = 0; j < y_; ++j) {
                fprintf(ff, "%f ", array_[x_ * i + j]);
            }
            fprintf(ff, "\n");
        }
        fclose(ff);
    }

    double *GetArray() const
    {
        return array_.get();
    }

private:
    std::unique_ptr<double[]> array_;
    size_t x_ {0};
    size_t y_ {0};
};

#endif // PARALLELPROGRAMMING_2_SEMESTER_1_LAB_ARRAY_HANDLER_H
