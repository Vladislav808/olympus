#include "olympus/Vector.hpp"

#include <cmath>

namespace olympus {

Vector::Vector(std::size_t size, double value) : data_(size, value) {}

Vector::Vector(std::initializer_list<double> values) : data_(values) {}

Vector::Vector(const std::vector<double>& values) : data_(values) {}

std::size_t Vector::size() const noexcept {
    return data_.size();
}

double& Vector::operator[](std::size_t index) noexcept {
    return data_[index];
}

const double& Vector::operator[](std::size_t index) const noexcept {
    return data_[index];
}

double& Vector::at(std::size_t index) {
    return data_.at(index);
}

const double& Vector::at(std::size_t index) const {
    return data_.at(index);
}

Vector Vector::operator+(const Vector& other) const {
    if (data_.size() != other.data_.size()) {
        throw std::invalid_argument("Vector sizes must match for addition");
    }

    Vector result(data_.size());
    for (std::size_t i = 0; i < data_.size(); ++i) {
        result.data_[i] = data_[i] + other.data_[i];
    }
    return result;
}

Vector Vector::operator-(const Vector& other) const {
    if (data_.size() != other.data_.size()) {
        throw std::invalid_argument("Vector sizes must match for subtraction");
    }

    Vector result(data_.size());
    for (std::size_t i = 0; i < data_.size(); ++i) {
        result.data_[i] = data_[i] - other.data_[i];
    }
    return result;
}

Vector Vector::operator*(double scalar) const {
    Vector result(data_.size());
    for (std::size_t i = 0; i < data_.size(); ++i) {
        result.data_[i] = data_[i] * scalar;
    }
    return result;
}

double Vector::dot(const Vector& other) const {
    if (data_.size() != other.data_.size()) {
        throw std::invalid_argument("Vector sizes must match for dot product");
    }

    double result = 0.0;
    for (std::size_t i = 0; i < data_.size(); ++i) {
        result += data_[i] * other.data_[i];
    }
    return result;
}

double Vector::norm() const {
    return std::sqrt(dot(*this));
}

Vector operator*(double scalar, const Vector& vector) {
    return vector * scalar;
}

} // namespace olympus
