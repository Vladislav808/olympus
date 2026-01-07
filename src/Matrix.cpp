#include "olympus/Matrix.hpp"

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

Matrix Matrix::transpose() const {
    Matrix result(cols_, rows_);
    for (std::size_t row = 0; row < rows_; ++row) {
        for (std::size_t col = 0; col < cols_; ++col) {
            result.data_[result.index(col, row)] = data_[index(row, col)];
        }
    }
    return result;
}

std::size_t Matrix::index(std::size_t row, std::size_t col) const {
    return row * cols_ + col;
}

Matrix operator*(double scalar, const Matrix& matrix) {
    return matrix * scalar;
}

} // namespace olympus
