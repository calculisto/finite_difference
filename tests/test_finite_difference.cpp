#include <doctest/doctest.h>
#include <calculisto/finite_difference/finite_difference.hpp>
    using calculisto::finite_difference::central_finite_difference;
#include <cmath>
#include <cstddef>

namespace
{
        // x^N, computed by repeated multiplication: exact on the dyadic
        // values used below.
        template <int N>
        constexpr auto
    monomial = [] (double x)
    {
            double
        r = 1.0;
        for (auto i = 0; i < N; ++i)
        {
            r *= x;
        }
        return r;
    };

        // The derivative of order K of x^N.
        template <int N, int K>
        double
    monomial_derivative (double x)
    {
            double
        r = 1.0;
        for (auto i = 0; i < K; ++i)
        {
            r *= N - i;
        }
        for (auto i = 0; i < N - K; ++i)
        {
            r *= x;
        }
        return r;
    }

        // A central scheme of order K and accuracy P is exact for the
        // polynomials of degree less than K + P. With x = 1 and h = 1/2, all
        // the values are exactly representable, so the result is exact.
        template <std::size_t K, std::size_t P, int N>
        void
    check_exact ()
    {
            volatile double
        x = 1.0;
            volatile double
        h = 0.5;
        CAPTURE (N);
        CHECK (central_finite_difference <0, K, P> (monomial <N>, h, x)
            == monomial_derivative <N, int (K)> (x));
    }

        template <std::size_t K, std::size_t P, int... N>
        void
    check_exact_up_to (std::integer_sequence <int, N...>)
    {
        (check_exact <K, P, N> (), ...);
    }
} // namespace

TEST_CASE ("exact for the polynomials of degree less than order + accuracy")
{
    check_exact_up_to <1, 2> (std::make_integer_sequence <int, 3> {});
    check_exact_up_to <1, 4> (std::make_integer_sequence <int, 5> {});
    check_exact_up_to <2, 2> (std::make_integer_sequence <int, 4> {});
    check_exact_up_to <2, 4> (std::make_integer_sequence <int, 6> {});
}

TEST_CASE ("leading error term, for the polynomial of degree order + accuracy")
{
    // The error is c h^P f^(K+P) (x), with c = 1/6, -1/30, 1/12 and -1/90.
    // With x = 1 and h = 1/2, it is exact.
        volatile double
    x = 1.0;
        volatile double
    h = 0.5;
    CHECK (central_finite_difference <0, 1, 2> (monomial <3>, h, x)
        == 3.0 + 0.25);
    CHECK (central_finite_difference <0, 1, 4> (monomial <5>, h, x)
        == 5.0 - 0.25);
    CHECK (central_finite_difference <0, 2, 2> (monomial <4>, h, x)
        == 12.0 + 0.5);
    CHECK (central_finite_difference <0, 2, 4> (monomial <6>, h, x)
        == 30.0 - 0.5);
}

TEST_CASE ("default order and accuracy: first derivative, second order")
{
        volatile double
    x = 1.0;
        volatile double
    h = 0.5;
    CHECK (central_finite_difference (monomial <3>, h, x) == 3.0 + 0.25);
}

TEST_CASE ("derivative with respect to the argument of the given index")
{
    // f (x, y) = x^3 y^5 + x y, exact on dyadic values.
        const auto
    f = [] (double x, double y)
    {
        return monomial <3> (x) * monomial <5> (y) + x * y;
    };
        volatile double
    x = 1.0;
        volatile double
    y = 0.5;
        volatile double
    h = 0.25;
    // Degree 3 in x: the fourth order schemes are exact.
    CHECK (central_finite_difference <0, 1, 4> (f, h, x, y)
        == 3.0 * x * x * monomial <5> (y) + y);
    CHECK (central_finite_difference <0, 2, 4> (f, h, x, y)
        == 6.0 * x * monomial <5> (y));
    // Degree 5 in y: the sixth derivative vanishes, so the fourth order
    // scheme is exact for the second derivative...
    CHECK (central_finite_difference <1, 2, 4> (f, h, x, y)
        == x * x * x * 20.0 * monomial <3> (y));
    // ... but not for the first derivative: error -h^4/30 * 120 x^3.
    CHECK (central_finite_difference <1, 1, 4> (f, h, x, y)
        == x * x * x * 5.0 * monomial <4> (y) + x - 4.0 * h * h * h * h * x * x * x);
}

TEST_CASE ("convergence, compared with std::exp")
{
    // Halving the step divides the error by about 2^P.
        volatile double
    x = 1.0;
        const auto
    f = [] (double t) { return std::exp (t); };
        const double
    exact = std::exp (x);
        const auto
    ratio = [&] (auto scheme)
    {
        return (scheme (0.1) - exact) / (scheme (0.05) - exact);
    };
    CHECK (ratio ([&] (double h) { return central_finite_difference <0, 1, 2> (f, h, x); })
        == doctest::Approx { 4.0 }.epsilon (0.01));
    CHECK (ratio ([&] (double h) { return central_finite_difference <0, 1, 4> (f, h, x); })
        == doctest::Approx { 16.0 }.epsilon (0.01));
    CHECK (ratio ([&] (double h) { return central_finite_difference <0, 2, 2> (f, h, x); })
        == doctest::Approx { 4.0 }.epsilon (0.01));
    CHECK (ratio ([&] (double h) { return central_finite_difference <0, 2, 4> (f, h, x); })
        == doctest::Approx { 16.0 }.epsilon (0.01));
}
