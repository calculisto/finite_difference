#include <doctest/doctest.h>
#include <calculisto/finite_difference/finite_difference.hpp>
    using calculisto::finite_difference::central_finite_difference;
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>

namespace
{
        // Whether central_finite_difference <Index> (a...) is a valid call.
        template <std::size_t Index, class... A>
        concept
    callable = requires (A&&... a)
    {
        central_finite_difference <Index> (std::forward <A> (a)...);
    };

        // Whether central_finite_difference <0, Order, Accuracy> (a...) is a
        // valid call.
        template <std::size_t Order, std::size_t Accuracy, class... A>
        concept
    callable_scheme = requires (A&&... a)
    {
        central_finite_difference <0, Order, Accuracy> (
              std::forward <A> (a)...
        );
    };

        constexpr auto
    identity = [] (auto x) { return x; };

        constexpr auto
    identity_2 = [] (auto x, auto) { return x; };

        using
    identity_t = decltype (identity);

        using
    identity_2_t = decltype (identity_2);
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

        // The ratio of the errors of the scheme of order K and accuracy P on
        // exp at 1, with the steps 0.1 and 0.05.
        template <std::size_t K, std::size_t P>
        double
    convergence_ratio ()
    {
            volatile double
        x = 1.0;
            const auto
        f = [] (double t) { return std::exp (t); };
            const double
        exact = std::exp (x);
            const auto
        error = [&] (double h)
        {
            return central_finite_difference <0, K, P> (f, h, x) - exact;
        };
        return error (0.1) / error (0.05);
    }
} // namespace

TEST_CASE ("exact for the polynomials of degree less than order + accuracy")
{
    check_exact_up_to <1, 2> (std::make_integer_sequence <int, 3> {});
    check_exact_up_to <1, 4> (std::make_integer_sequence <int, 5> {});
    check_exact_up_to <2, 2> (std::make_integer_sequence <int, 4> {});
    check_exact_up_to <2, 4> (std::make_integer_sequence <int, 6> {});
}

TEST_CASE ("leading error term for the polynomial of degree order + accuracy")
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

TEST_CASE ("default order and accuracy: first derivative to second order")
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
        == x * x * x * 5.0 * monomial <4> (y) + x
            - 4.0 * h * h * h * h * x * x * x);
}

TEST_CASE ("convergence compared with std::exp")
{
    // Halving the step divides the error by about 2^P.
    CHECK (convergence_ratio <1, 2> ()
        == doctest::Approx { 4.0 }.epsilon (0.01));
    CHECK (convergence_ratio <1, 4> ()
        == doctest::Approx { 16.0 }.epsilon (0.01));
    CHECK (convergence_ratio <2, 2> ()
        == doctest::Approx { 4.0 }.epsilon (0.01));
    CHECK (convergence_ratio <2, 4> ()
        == doctest::Approx { 16.0 }.epsilon (0.01));
}

TEST_CASE ("the function is called as an lvalue: it is called several times")
{
    // The rvalue overload stands for one that consumes the object (e.g.
    // moves from a member): calling it more than once would be a use after
    // move.
        struct
    function_t
    {
            double
        operator () (double x) const&
        {
            return x * x;
        }
            double
        operator () (double) &&
        {
            return 0.0;
        }
    };
        volatile double
    x = 1.0;
        volatile double
    h = 0.5;
    CHECK (central_finite_difference <0, 1, 2> (function_t {}, h, x) == 2.0);
    CHECK (central_finite_difference <0, 1, 4> (function_t {}, h, x) == 2.0);
    CHECK (central_finite_difference <0, 2, 2> (function_t {}, h, x) == 2.0);
    CHECK (central_finite_difference <0, 2, 4> (function_t {}, h, x) == 2.0);
}

TEST_CASE ("the shifted argument keeps its type: no narrowing")
{
    // A double step would be truncated by an int argument, and rounded by a
    // float one (1000.f + 1e-6 == 1000.f).
    static_assert (not callable <0, identity_t, double, int>);
    static_assert (not callable <0, identity_t, double, float>);
    static_assert (not callable <1, identity_2_t, double, double, int>);
    static_assert (callable <0, identity_t, double, double>);
    static_assert (callable <0, identity_t, float, float>);
    static_assert (callable <0, identity_t, int, int>);
    static_assert (callable <1, identity_2_t, double, int, double>);
    // A float step on a double argument is exact.
    static_assert (callable <0, identity_t, float, double>);
        volatile double
    x = 1.0;
        volatile float
    h = 0.5f;
    CHECK (central_finite_difference (monomial <2>, h, x) == 2.0);
}

TEST_CASE ("only the implemented schemes are callable")
{
    static_assert (callable_scheme <1, 2, identity_t, double, double>);
    static_assert (callable_scheme <1, 4, identity_t, double, double>);
    static_assert (callable_scheme <2, 2, identity_t, double, double>);
    static_assert (callable_scheme <2, 4, identity_t, double, double>);
    static_assert (not callable_scheme <0, 2, identity_t, double, double>);
    static_assert (not callable_scheme <3, 2, identity_t, double, double>);
    static_assert (not callable_scheme <1, 0, identity_t, double, double>);
    static_assert (not callable_scheme <1, 3, identity_t, double, double>);
    static_assert (not callable_scheme <1, 6, identity_t, double, double>);
}

TEST_CASE ("the index must designate an argument")
{
    static_assert (not callable <1, identity_t, double, double>);
    static_assert (not callable <2, identity_2_t, double, double, double>);
    static_assert (not callable <0, identity_t, double>);
}

TEST_CASE ("the function must be callable with the arguments")
{
        constexpr auto
    takes_string = [] (std::string const&) { return 0.0; };
        constexpr auto
    takes_two = [] (double, double) { return 0.0; };
        struct
    rvalue_only
    {
            double
        operator () (double x) &&
        {
            return x;
        }
    };
    static_assert (not callable <0, decltype (takes_string), double, double>);
    static_assert (not callable <0, decltype (takes_two), double, double>);
    // The function is called several times: it must be callable as an
    // lvalue.
    static_assert (not callable <0, rvalue_only, double, double>);
}

TEST_CASE ("the values of the function support the arithmetic of the scheme")
{
        constexpr auto
    returns_void = [] (double) {};
        constexpr auto
    returns_string = [] (double) { return std::string {}; };
    static_assert (not callable <0, decltype (returns_void), double, double>);
    static_assert (not callable <0, decltype (returns_string), double, double>);
}
