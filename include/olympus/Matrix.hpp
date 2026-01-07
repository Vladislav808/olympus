#ifndef OLYMPUS_MATRIX_HPP
#define OLYMPUS_MATRIX_HPP

#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <vector>

namespace olympus {

/**
 * @brief Dense matrix of doubles with basic linear algebra operations.
 */
class Matrix {
public:
    /// @brief Result container for LU decomposition (L, U, and permutation P).
    struct LUDecomposition;

    /// @brief Proxy for mutable row access via operator[].
    class RowProxy {
    public:
        RowProxy(double* row_start, std::size_t cols);
        double& operator[](std::size_t col) noexcept;
        const double& operator[](std::size_t col) const noexcept;

    private:
        double* row_start_;
        std::size_t cols_;
    };

    /// @brief Proxy for const row access via operator[].
    class ConstRowProxy {
    public:
        ConstRowProxy(const double* row_start, std::size_t cols);
        const double& operator[](std::size_t col) const noexcept;

    private:
        const double* row_start_;
        std::size_t cols_;
    };

    /// @brief Create an empty matrix.
    Matrix() = default;
    /// @brief Create a matrix with fixed size and optional fill value.
    Matrix(std::size_t rows, std::size_t cols, double value = 0.0);
    /// @brief Create a matrix from nested initializer lists.
    Matrix(std::initializer_list<std::initializer_list<double>> values);
    /// @brief Create a matrix from a flat vector in row-major order.
    Matrix(std::size_t rows, std::size_t cols, const std::vector<double>& values);

    /// @brief Number of rows.
    std::size_t rows() const noexcept;
    /// @brief Number of columns.
    std::size_t cols() const noexcept;

    /// @brief Mutable row access (unchecked).
    RowProxy operator[](std::size_t row) noexcept;
    /// @brief Const row access (unchecked).
    ConstRowProxy operator[](std::size_t row) const noexcept;

    /// @brief Element access with bounds checking.
    double& at(std::size_t row, std::size_t col);
    /// @brief Element access with bounds checking.
    const double& at(std::size_t row, std::size_t col) const;

    /// @brief Element-wise addition.
    Matrix operator+(const Matrix& other) const;
    /// @brief Element-wise subtraction.
    Matrix operator-(const Matrix& other) const;
    /// @brief Scalar multiplication.
    Matrix operator*(double scalar) const;
    /// @brief Matrix multiplication.
    Matrix operator*(const Matrix& other) const;

    /// @brief Append a row at the bottom of the matrix.
    void addRow(const std::vector<double>& values);
    /// @brief Append a column at the right of the matrix.
    void addColumn(const std::vector<double>& values);
    /// @brief Remove a row by index.
    void removeRow(std::size_t row);
    /// @brief Remove a column by index.
    void removeColumn(std::size_t col);

    /// @brief Return the transposed matrix.
    Matrix transpose() const;
    /// @brief Compute LU decomposition with partial pivoting.
    LUDecomposition luDecomposition() const;
    /// @brief Compute matrix inverse (square matrices only).
    Matrix inverse() const;

private:
    std::size_t rows_ = 0;
    std::size_t cols_ = 0;
    std::vector<double> data_;

    std::size_t index(std::size_t row, std::size_t col) const;
};

struct Matrix::LUDecomposition {
    Matrix L;
    Matrix U;
    Matrix P;
};

/// @brief Scalar multiplication with scalar on the left.
Matrix operator*(double scalar, const Matrix& matrix);

} // namespace olympus

#endif // OLYMPUS_MATRIX_HPP
