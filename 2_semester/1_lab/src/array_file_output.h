#ifndef PARALLELPROGRAMMING_2_SEMESTER_1_LAB_ARRAY_HANDLER_H
#define PARALLELPROGRAMMING_2_SEMESTER_1_LAB_ARRAY_HANDLER_H

#include <cstddef>
#include <memory>

class ArrayHandler {
public:
    explicit ArrayHandler(size_t x, size_t y) : isize_(x), jsize_(y)
    {
        array_ = std::make_unique<double[]>(isize_ * jsize_);
    }

    ~ArrayHandler() = default;

    void DefaultFillIn()
    {
        for (size_t i = 0; i < isize_; ++i) {
            for (size_t j = 0; j < jsize_; ++j) {
                array_[jsize_ * i + j] = 10 * i + j;
            }
        }
    }

    void ZeroFillIn()
    {
        for (size_t i = 0; i < isize_; ++i) {
            for (size_t j = 0; j < jsize_; ++j) {
                array_[jsize_ * i + j] = 0;
            }
        }
    }

    void DumpToFile(const char *filename) const
    {
        FILE *ff = nullptr;
        ff = fopen(filename, "w");

        for (size_t i = 0; i < isize_; ++i) {
            for (size_t j = 0; j < jsize_; ++j) {
                fprintf(ff, "%f ", array_[jsize_ * i + j]);
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
    size_t isize_ {0};
    size_t jsize_ {0};
};

#endif // PARALLELPROGRAMMING_2_SEMESTER_1_LAB_ARRAY_HANDLER_H
