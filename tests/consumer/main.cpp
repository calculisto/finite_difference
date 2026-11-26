#include <calculisto/finite_difference/finite_difference.hpp>
    using calculisto::finite_difference::central_finite_difference;

    int
main ()
{
    // The second order scheme is exact for a quadratic function.
        const auto
    d = central_finite_difference ([] (double x) { return x * x; }, 0.5, 1.0);
    return d == 2.0 ? 0 : 1;
}
