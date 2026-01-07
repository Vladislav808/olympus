# olympus

Array (matrix and vector) library for C++.

## Build

```sh
cmake -S . -B build
cmake --build build
```

## Quick example

```cpp
#include "olympus/Matrix.hpp"
#include "olympus/Vector.hpp"

olympus::Vector v{1.0, 2.0, 3.0};
olympus::Vector w{4.0, 5.0, 6.0};

auto sum = v + w;
auto dot = v.dot(w);

auto m = olympus::Matrix{{1.0, 2.0}, {3.0, 4.0}};
auto t = m.transpose();
```
