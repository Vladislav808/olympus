#ifndef OLYMPUS_VECTOR_HPP
#define OLYMPUS_VECTOR_HPP

#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <vector>

namespace olympus {

/**
 * @brief Dense vector of doubles with basic linear algebra operations.
 */
class Vector {
public:
    /// @brief Create an empty vector.
    Vector() = default;
    /// @brief Create a vector with a fixed size and optional fill value.
    explicit Vector(std::size_t size, double value = 0.0);
    /// @brief Create a vector from an initializer list.
    Vector(std::initializer_list<double> values);
    /// @brief Create a vector from a std::vector.
    explicit Vector(const std::vector<double>& values);

    /// @brief Return the number of elements.
    std::size_t size() const noexcept;

    /// @brief Access an element without bounds checking.
    double& operator[](std::size_t index) noexcept;
    /// @brief Access an element without bounds checking.
    const double& operator[](std::size_t index) const noexcept;

    /// @brief Access an element with bounds checking.
    double& at(std::size_t index);
    /// @brief Access an element with bounds checking.
    const double& at(std::size_t index) const;

    /// @brief Element-wise addition.
    Vector operator+(const Vector& other) const;
    /// @brief Element-wise subtraction.
    Vector operator-(const Vector& other) const;
    /// @brief Scalar multiplication.
    Vector operator*(double scalar) const;

    /// @brief Dot product with another vector.
    double dot(const Vector& other) const;
    /// @brief Euclidean norm (L2).
    double norm() const;

private:
    std::vector<double> data_;
};

/// @brief Scalar multiplication with scalar on the left.
Vector operator*(double scalar, const Vector& vector);

} // namespace olympus

#endif // OLYMPUS_VECTOR_HPP
