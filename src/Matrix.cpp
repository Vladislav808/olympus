#include "olympus/Matrix.hpp"

#include <cmath>

namespace olympus {

Matrix::RowProxy::RowProxy(double* row_start, std::size_t cols)
    : row_start_(row_start), cols_(cols) {}

double& Matrix::RowProxy::operator[](std::size_t col) noexcept {
    return row_start_[col];
}

const double& Matrix::RowProxy::operator[](std::size_t col) const noexcept {
    return row_start_[col];
}

Matrix::ConstRowProxy::ConstRowProxy(const double* row_start, std::size_t cols)
    : row_start_(row_start), cols_(cols) {}

const double& Matrix::ConstRowProxy::operator[](std::size_t col) const noexcept {
    return row_start_[col];
}

Matrix::Matrix(std::size_t rows, std::size_t cols, double value)
    : rows_(rows), cols_(cols), data_(rows * cols, value) {}

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> values) {
    rows_ = values.size();
    cols_ = rows_ == 0 ? 0 : values.begin()->size();
    data_.reserve(rows_ * cols_);

    for (const auto& row : values) {
        if (row.size() != cols_) {
            throw std::invalid_argument("All rows must have the same length");
        }
        data_.insert(data_.end(), row.begin(), row.end());
    }
}

Matrix::Matrix(std::size_t rows, std::size_t cols, const std::vector<double>& values)
    : rows_(rows), cols_(cols), data_(values) {
    if (values.size() != rows * cols) {
        throw std::invalid_argument("Value count must match rows * cols");
    }
}

std::size_t Matrix::rows() const noexcept {
    return rows_;
}

std::size_t Matrix::cols() const noexcept {
    return cols_;
}

Matrix::RowProxy Matrix::operator[](std::size_t row) noexcept {
    return RowProxy(data_.data() + index(row, 0), cols_);
}

Matrix::ConstRowProxy Matrix::operator[](std::size_t row) const noexcept {
    return ConstRowProxy(data_.data() + index(row, 0), cols_);
}

double& Matrix::at(std::size_t row, std::size_t col) {
    if (row >= rows_ || col >= cols_) {
        throw std::out_of_range("Matrix indices are out of range");
    }
    return data_[index(row, col)];
}

const double& Matrix::at(std::size_t row, std::size_t col) const {
    if (row >= rows_ || col >= cols_) {
        throw std::out_of_range("Matrix indices are out of range");
    }
    return data_[index(row, col)];
}

Matrix Matrix::operator+(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::invalid_argument("Matrix sizes must match for addition");
    }

    Matrix result(rows_, cols_);
    for (std::size_t i = 0; i < data_.size(); ++i) {
        result.data_[i] = data_[i] + other.data_[i];
    }
    return result;
}

Matrix Matrix::operator-(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::invalid_argument("Matrix sizes must match for subtraction");
    }

    Matrix result(rows_, cols_);
    for (std::size_t i = 0; i < data_.size(); ++i) {
        result.data_[i] = data_[i] - other.data_[i];
    }
    return result;
}

Matrix Matrix::operator*(double scalar) const {
    Matrix result(rows_, cols_);
    for (std::size_t i = 0; i < data_.size(); ++i) {
        result.data_[i] = data_[i] * scalar;
    }
    return result;
}

Matrix Matrix::operator*(const Matrix& other) const {
    if (cols_ != other.rows_) {
        throw std::invalid_argument("Matrix sizes are incompatible for multiplication");
    }

    Matrix result(rows_, other.cols_, 0.0);
    for (std::size_t row = 0; row < rows_; ++row) {
        for (std::size_t col = 0; col < other.cols_; ++col) {
            double sum = 0.0;
            for (std::size_t k = 0; k < cols_; ++k) {
                sum += data_[index(row, k)] * other.data_[other.index(k, col)];
            }
            result.data_[result.index(row, col)] = sum;
        }
    }
    return result;
}

void Matrix::addRow(const std::vector<double>& values) {
    if (rows_ == 0 && cols_ == 0) {
        cols_ = values.size();
        rows_ = cols_ == 0 ? 0 : 1;
        data_ = values;
        return;
    }

    if (values.size() != cols_) {
        throw std::invalid_argument("Row size must match matrix column count");
    }

    data_.insert(data_.end(), values.begin(), values.end());
    ++rows_;
}

void Matrix::addColumn(const std::vector<double>& values) {
    if (rows_ == 0 && cols_ == 0) {
        rows_ = values.size();
        cols_ = rows_ == 0 ? 0 : 1;
        data_ = values;
        return;
    }

    if (values.size() != rows_) {
        throw std::invalid_argument("Column size must match matrix row count");
    }

    std::vector<double> updated;
    updated.reserve((cols_ + 1) * rows_);
    for (std::size_t row = 0; row < rows_; ++row) {
        const std::size_t row_start = row * cols_;
        updated.insert(updated.end(), data_.begin() + row_start, data_.begin() + row_start + cols_);
        updated.push_back(values[row]);
    }
    data_.swap(updated);
    ++cols_;
}

void Matrix::removeRow(std::size_t row) {
    if (row >= rows_) {
        throw std::out_of_range("Row index is out of range");
    }

    std::vector<double> updated;
    updated.reserve((rows_ - 1) * cols_);
    for (std::size_t current = 0; current < rows_; ++current) {
        if (current == row) {
            continue;
        }
        const std::size_t row_start = current * cols_;
        updated.insert(updated.end(), data_.begin() + row_start, data_.begin() + row_start + cols_);
    }
    data_.swap(updated);
    --rows_;
    if (rows_ == 0) {
        data_.clear();
    }
}

void Matrix::removeColumn(std::size_t col) {
    if (col >= cols_) {
        throw std::out_of_range("Column index is out of range");
    }

    if (cols_ == 1) {
        --cols_;
        data_.clear();
        return;
    }

    std::vector<double> updated;
    updated.reserve(rows_ * (cols_ - 1));
    for (std::size_t row = 0; row < rows_; ++row) {
        const std::size_t row_start = row * cols_;
        for (std::size_t current_col = 0; current_col < cols_; ++current_col) {
            if (current_col == col) {
                continue;
            }
            updated.push_back(data_[row_start + current_col]);
        }
    }
    data_.swap(updated);
    --cols_;
}

Matrix Matrix::transpose() const {
    Matrix result(cols_, rows_);
    for (std::size_t row = 0; row < rows_; ++row) {
        for (std::size_t col = 0; col < cols_; ++col) {
            result.data_[result.index(col, row)] = data_[index(row, col)];
        }
    }
    return result;
}

Matrix::LUDecomposition Matrix::luDecomposition() const {
    if (rows_ != cols_) {
        throw std::invalid_argument("LU decomposition requires a square matrix");
    }

    const std::size_t n = rows_;
    Matrix L(n, n, 0.0);
    Matrix U(*this);
    Matrix P(n, n, 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        L.data_[L.index(i, i)] = 1.0;
        P.data_[P.index(i, i)] = 1.0;
    }

    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot_row = k;
        double pivot_value = std::abs(U.data_[U.index(k, k)]);
        for (std::size_t row = k + 1; row < n; ++row) {
            double value = std::abs(U.data_[U.index(row, k)]);
            if (value > pivot_value) {
                pivot_value = value;
                pivot_row = row;
            }
        }

        if (pivot_value == 0.0) {
            throw std::runtime_error("Matrix is singular and cannot be decomposed");
        }

        if (pivot_row != k) {
            for (std::size_t col = 0; col < n; ++col) {
                std::swap(U.data_[U.index(k, col)], U.data_[U.index(pivot_row, col)]);
                std::swap(P.data_[P.index(k, col)], P.data_[P.index(pivot_row, col)]);
            }
            for (std::size_t col = 0; col < k; ++col) {
                std::swap(L.data_[L.index(k, col)], L.data_[L.index(pivot_row, col)]);
            }
        }

        for (std::size_t row = k + 1; row < n; ++row) {
            double factor = U.data_[U.index(row, k)] / U.data_[U.index(k, k)];
            L.data_[L.index(row, k)] = factor;
            for (std::size_t col = k; col < n; ++col) {
                U.data_[U.index(row, col)] -= factor * U.data_[U.index(k, col)];
            }
        }
    }

    return {L, U, P};
}

Matrix Matrix::inverse() const {
    if (rows_ != cols_) {
        throw std::invalid_argument("Matrix inverse requires a square matrix");
    }

    const std::size_t n = rows_;
    LUDecomposition decomposition = luDecomposition();
    const Matrix& L = decomposition.L;
    const Matrix& U = decomposition.U;
    const Matrix& P = decomposition.P;

    Matrix inverse_matrix(n, n, 0.0);
    std::vector<double> y(n);
    std::vector<double> x(n);

    for (std::size_t col = 0; col < n; ++col) {
        for (std::size_t row = 0; row < n; ++row) {
            double sum = P.data_[P.index(row, col)];
            for (std::size_t k = 0; k < row; ++k) {
                sum -= L.data_[L.index(row, k)] * y[k];
            }
            y[row] = sum;
        }

        for (std::size_t row = n; row-- > 0;) {
            double sum = y[row];
            for (std::size_t k = row + 1; k < n; ++k) {
                sum -= U.data_[U.index(row, k)] * x[k];
            }
            double pivot = U.data_[U.index(row, row)];
            if (pivot == 0.0) {
                throw std::runtime_error("Matrix is singular and cannot be inverted");
            }
            x[row] = sum / pivot;
        }

        for (std::size_t row = 0; row < n; ++row) {
            inverse_matrix.data_[inverse_matrix.index(row, col)] = x[row];
        }
    }

    return inverse_matrix;
}

std::size_t Matrix::index(std::size_t row, std::size_t col) const {
    return row * cols_ + col;
}

Matrix operator*(double scalar, const Matrix& matrix) {
    return matrix * scalar;
}

} // namespace olympus
