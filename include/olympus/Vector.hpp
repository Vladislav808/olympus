#ifndef OLYMPUS_VECTOR_HPP
#define OLYMPUS_VECTOR_HPP

#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <vector>

namespace olympus {

class Vector {
public:
    Vector() = default;
    explicit Vector(std::size_t size, double value = 0.0);
    Vector(std::initializer_list<double> values);
    explicit Vector(const std::vector<double>& values);

    std::size_t size() const noexcept;

    double& operator[](std::size_t index) noexcept;
    const double& operator[](std::size_t index) const noexcept;

    double& at(std::size_t index);
    const double& at(std::size_t index) const;

    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;
    Vector operator*(double scalar) const;

    double dot(const Vector& other) const;
    double norm() const;

private:
    std::vector<double> data_;
};

Vector operator*(double scalar, const Vector& vector);

} // namespace olympus

#endif // OLYMPUS_VECTOR_HPP
