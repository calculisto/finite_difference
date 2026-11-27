# calculisto::finite_difference
Central finite difference approximations of the first and second derivatives
of functions of several variables, to second or fourth order of accuracy.

## Usage
```cpp
#include <calculisto/finite_difference/finite_difference.hpp>
    using calculisto::finite_difference::central_finite_difference;

    const auto
f = [] (double x, double y) { return x * x * y * y; };
// d f / d x = 8 at (1, 2): the first derivative with respect to argument 0,
// to second order, with the step 1e-5.
    const auto
dfdx = central_finite_difference (f, 1e-5, 1.0, 2.0);
// d^2 f / d y^2 = 2 at (1, 2), to fourth order, with the step 1e-3.
    const auto
d2fdy2 = central_finite_difference <1, 2, 4> (f, 1e-3, 1.0, 2.0);
```
`central_finite_difference <Index, Order, Accuracy> (f, h, x...)` approximates
the derivative of order `Order` (1 or 2) of `f` with respect to its argument
number `Index`, at `x...`, with the step `h`, by a central scheme of accuracy
`Accuracy` (2 or 4): the error is proportional to `h` to the power
`Accuracy`. The defaults are `Index = 0`, `Order = 1` and `Accuracy = 2`.

The call is valid (and only then, so that it can be checked with concepts)
when:
- `Index` designates an argument, and the order and the accuracy are
  implemented;
- the shifted argument keeps its type: `x + k * h` has the type of `x`. With
  an `int` argument, a `double` step would be truncated, and with a `float`
  one it would be rounded: such calls are rejected;
- `f` can be called as an lvalue with the arguments as const lvalues, since it
  is called several times, and its values support the arithmetic of the
  scheme.

### Choosing the step
A smaller step reduces the truncation error, proportional to `h` to the power
`Accuracy`, but increases the rounding error, proportional to the machine
epsilon divided by `h` to the power `Order`. For arguments and derivatives of
order 1, the total error is smallest for a step of about epsilon to the power
1 / (`Order` + `Accuracy`); scale it with the magnitude of the argument. In
`double`, measured on `exp` and `sin` on [0.5, 2]:

| Order | Accuracy | Suggested step | Measured best step | Relative error |
|-------|----------|----------------|--------------------|----------------|
| 1     | 2        | 6e-6           | 8e-6               | 3e-11          |
| 1     | 4        | 7e-4           | 1e-3               | 3e-13          |
| 2     | 2        | 1e-4           | 2e-4               | 4e-9           |
| 2     | 4        | 3e-3           | 5e-3               | 8e-12          |

A step 10 times too small or too large typically costs one to three digits.

## Dependencies
- A C++20 compiler.
- [CMake](https://cmake.org) 3.25 or later, to install the library or run the
  tests.
- To run the tests, the [doctest](https://github.com/doctest/doctest) library,
  version 2.5 or later. If it isn't installed, CMake downloads it.

## Installation
This is a header-only library. It can be used in any of these ways:
- copy the `include` directory where your compiler can find it;
- from a CMake project, with `add_subdirectory` or `FetchContent`, then link to
  the `calculisto::finite_difference` target;
- install it with CMake:
  ```
  cmake --preset dev
  cmake --install build/dev --prefix <directory>
  ```
  then, from a CMake project:
  ```cmake
  find_package (calculisto_finite_difference 0.1 REQUIRED)
  target_link_libraries (<target> PRIVATE calculisto::finite_difference)
  ```
  adding `<directory>` to `CMAKE_PREFIX_PATH` if it isn't a standard location.

## Tests
To configure, build and run the tests with the default compiler, execute
```
cmake --workflow --preset dev
```
in the root directory of the project. The tests are compiled with warnings as
errors, and run with the address sanitizer and, except with MSVC and clang-cl,
the undefined behavior sanitizer. They also check that the installed library
can be used with `find_package`.

Other presets select a compiler: `gcc`, `clang`, `clang-libcxx`, and on
Windows `msvc` and `clang-cl`. `cmake --list-presets=workflow` lists them.

## License
SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

## Affiliation
This material is developed by the Numerical Modelling platform at the
Institut des Sciences de la Terre d'Orléans (https://www.isto-orleans.fr),
a joint laboratory of the University of Orléans (https://www.univ-orleans.fr),
the french National Center For Scientific Research (https://www.cnrs.fr) and
the french Geological Survey (https://www.brgm.fr).

[![logo ISTO][isto-logo]](https://www.isto-orleans.fr) &emsp;
[![logo CNRS][cnrs-logo]](https://www.cnrs.fr) &emsp;
[![logo UO][uo-logo]](https://www.univ-orleans.fr) &emsp;
[![logo BRGM][brgm-logo]](https://www.brgm.fr)

[isto-logo]: https://calcul-isto.cnrs-orleans.fr/logos/isto-156.png
[cnrs-logo]: https://calcul-isto.cnrs-orleans.fr/logos/cnrs-128.png
[uo-logo]: https://calcul-isto.cnrs-orleans.fr/logos/uo-128.png
[brgm-logo]: https://calcul-isto.cnrs-orleans.fr/logos/brgm-256.png
