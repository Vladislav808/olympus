#ifndef OLYMPUS_MATRIX_HPP
#define OLYMPUS_MATRIX_HPP

#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <vector>

namespace olympus {

class Matrix {
public:
    struct LUDecomposition {
        Matrix L;
        Matrix U;
        Matrix P;
    };

    class RowProxy {
    public:
        RowProxy(double* row_start, std::size_t cols);
        double& operator[](std::size_t col) noexcept;
        const double& operator[](std::size_t col) const noexcept;

    private:
        double* row_start_;
        std::size_t cols_;
    };

    class ConstRowProxy {
    public:
        ConstRowProxy(const double* row_start, std::size_t cols);
        const double& operator[](std::size_t col) const noexcept;

    private:
        const double* row_start_;
        std::size_t cols_;
    };

    Matrix() = default;
    Matrix(std::size_t rows, std::size_t cols, double value = 0.0);
    Matrix(std::initializer_list<std::initializer_list<double>> values);
    Matrix(std::size_t rows, std::size_t cols, const std::vector<double>& values);

    std::size_t rows() const noexcept;
    std::size_t cols() const noexcept;

    RowProxy operator[](std::size_t row) noexcept;
    ConstRowProxy operator[](std::size_t row) const noexcept;

    double& at(std::size_t row, std::size_t col);
    const double& at(std::size_t row, std::size_t col) const;

    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(double scalar) const;
    Matrix operator*(const Matrix& other) const;

    void addRow(const std::vector<double>& values);
    void addColumn(const std::vector<double>& values);
    void removeRow(std::size_t row);
    void removeColumn(std::size_t col);

    Matrix transpose() const;
    LUDecomposition luDecomposition() const;
    Matrix inverse() const;

private:
    std::size_t rows_ = 0;
    std::size_t cols_ = 0;
    std::vector<double> data_;

    std::size_t index(std::size_t row, std::size_t col) const;
};

Matrix operator*(double scalar, const Matrix& matrix);

} // namespace olympus

#endif // OLYMPUS_MATRIX_HPP
