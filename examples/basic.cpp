#include "olympus/Matrix.hpp"
#include "olympus/Vector.hpp"

#include <iostream>

int main() {
    olympus::Vector a{1.0, 2.0, 3.0};
    olympus::Vector b{4.0, 5.0, 6.0};

    olympus::Vector sum = a + b;
    std::cout << "Dot product: " << a.dot(b) << "\n";
    std::cout << "Norm of sum: " << sum.norm() << "\n";

    olympus::Matrix left{{1.0, 2.0}, {3.0, 4.0}};
    olympus::Matrix right{{5.0, 6.0}, {7.0, 8.0}};
    olympus::Matrix product = left * right;
    olympus::Matrix transposed = product.transpose();

    std::cout << "Transposed(0, 1): " << transposed.at(0, 1) << "\n";
    return 0;
}
