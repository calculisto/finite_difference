#include <doctest/doctest.h>
#include <calculisto/finite_difference/finite_difference.hpp>
    using namespace calculisto::finite_difference;
    using doctest::Approx;

#include <cmath>

TEST_CASE("finite difference")
{
    SUBCASE("Unary function")
    {
            const auto
        f = [](auto x)
        {
            return x * x * x;
        };
            const auto
        d1f = [](auto x)
        {
            return 3 * x * x;
        };
            const auto
        d2f = [](auto x)
        {
            return 6 * x;
        };
            const auto
        x = 7.36;
            const auto
        expected_dy_1 = d1f (x);
            const auto
        dy_1_2 = central_finite_difference           (f, 1e-6, x);
        CHECK(dy_1_2 == Approx { expected_dy_1 });
            const auto
        dy_1_4 = central_finite_difference <0, 1, 4> (f, 1e-6, x);
        CHECK(dy_1_4 == Approx { expected_dy_1 });
            const auto
        expected_dy_2 = d2f (x);
            const auto
        dy_2_2 = central_finite_difference <0, 2> (f, 1e-2, x);
        CHECK(dy_2_2 == Approx { expected_dy_2 });
            const auto
        dy_2_4 = central_finite_difference <0, 2, 4> (f, 1e-2, x);
        CHECK(dy_2_4 == Approx { expected_dy_2 });
    }
    SUBCASE("Binary function")
    {
            const auto
        f = [](auto x, auto y)
        {
            return x * x * x + sin (y) + x * y;
        };
            const auto
        d1xf = [](auto x, auto y)
        {
            return 3 * x * x + y;
        };
            const auto
        d2xf = [](auto x, auto)
        {
            return 6 * x;
        };
            const auto
        d1yf = [](auto x, auto y)
        {
            return cos (y) + x;
        };
            const auto
        d2yf = [](auto, auto y)
        {
            return -sin (y);
        };
            const auto
        x = 7.36;
            const auto
        y = -1.658;
            const auto
        expected_dz_1x = d1xf (x, y);
            const auto
        dz_1x_2 = central_finite_difference           (f, 1e-6, x, y);
        CHECK(dz_1x_2 == Approx { expected_dz_1x });
            const auto
        dz_1x_4 = central_finite_difference <0, 1, 4> (f, 1e-6, x, y);
        CHECK(dz_1x_4 == Approx { expected_dz_1x });
            const auto
        expected_dz_2x = d2xf (x, y);
            const auto
        dz_2x_2 = central_finite_difference <0, 2> (f, 1e-2, x, y);
        CHECK(dz_2x_2 == Approx { expected_dz_2x });
            const auto
        dz_2x_4 = central_finite_difference <0, 2, 4> (f, 1e-2, x, y);
        CHECK(dz_2x_4 == Approx { expected_dz_2x });
            const auto
        expected_dz_1y = d1yf (x, y);
            const auto
        dz_1y_2 = central_finite_difference <1>       (f, 1e-6, x, y);
        CHECK(dz_1y_2 == Approx { expected_dz_1y });
            const auto
        dz_1y_4 = central_finite_difference <1, 1, 4> (f, 1e-6, x, y);
        CHECK(dz_1y_4 == Approx { expected_dz_1y });
            const auto
        expected_dz_2y = d2yf (x, y);
            const auto
        dz_2y_2 = central_finite_difference <1, 2> (f, 1e-2, x, y);
        CHECK(dz_2y_2 == Approx { expected_dz_2y });
            const auto
        dz_2y_4 = central_finite_difference <1, 2, 4> (f, 1e-2, x, y);
        CHECK(dz_2y_4 == Approx { expected_dz_2y });
    }
}

