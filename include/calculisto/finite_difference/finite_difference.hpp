#pragma once
#include <concepts>
#include <cstddef>
#include <iterator>
#include <tuple>
#include <type_traits>
#include <utility>

    namespace
calculisto::finite_difference
{
    namespace
detail
{
        // The arguments, with the one at Index shifted by delta.
        template <std::size_t Index, class T, class... Args>
        constexpr auto
    shift_by (T const& delta, Args const&... args)
    {
            std::tuple <Args...>
        r { args... };
        std::get <Index> (r) = std::get <Index> (r) + delta;
        return r;
    }

        // The central scheme for the derivative of order Order with the given
        // accuracy:
        //   f^(Order) (x) ~ sum (coefficients[i] * f (x + offsets[i] * h))
        //                   / denominator / h^Order
        template <std::size_t Order, std::size_t Accuracy>
        struct
    stencil;

        template <>
        struct
    stencil <1, 2>
    {
            static constexpr int
        offsets[] { -1, 1 };
            static constexpr int
        coefficients[] { -1, 1 };
            static constexpr int
        denominator = 2;
    };

        template <>
        struct
    stencil <1, 4>
    {
            static constexpr int
        offsets[] { -2, -1, 1, 2 };
            static constexpr int
        coefficients[] { 1, -8, 8, -1 };
            static constexpr int
        denominator = 12;
    };

        template <>
        struct
    stencil <2, 2>
    {
            static constexpr int
        offsets[] { -1, 0, 1 };
            static constexpr int
        coefficients[] { 1, -2, 1 };
            static constexpr int
        denominator = 1;
    };

        template <>
        struct
    stencil <2, 4>
    {
            static constexpr int
        offsets[] { -2, -1, 0, 1, 2 };
            static constexpr int
        coefficients[] { -1, 16, -30, 16, -1 };
            static constexpr int
        denominator = 12;
    };

        // The term i of the sum: the argument at Index is not shifted for a
        // zero offset, so that it is passed unchanged (e.g. -0.0 stays -0.0).
        template <std::size_t Index, class Stencil, std::size_t I>
        constexpr auto
    term (auto& F, auto const& delta, auto const&... args)
    {
            constexpr auto
        offset = Stencil::offsets[I];
            constexpr auto
        coefficient = Stencil::coefficients[I];
        if constexpr (offset == 0)
        {
            return coefficient * F (args...);
        }
        else
        {
            return coefficient * std::apply (F, shift_by <Index> (offset * delta, args...));
        }
    }

        template <std::size_t Index, std::size_t Order, std::size_t Accuracy>
        constexpr auto
    evaluate (auto& F, auto const& delta, auto const&... args)
    {
            using
        S = stencil <Order, Accuracy>;
            constexpr auto
        size = std::size (S::offsets);
            const auto
        sum = [&] <std::size_t... I> (std::index_sequence <I...>)
        {
            return (... + term <Index, S, I> (F, delta, args...));
        } (std::make_index_sequence <size> {});
        if constexpr (Order == 1)
        {
            return sum / S::denominator / delta;
        }
        else
        {
            return sum / S::denominator / delta / delta;
        }
    }

        // The argument X can be shifted by multiples of the step D without
        // changing its type: an int argument would truncate a double step, a
        // float one would round it.
        template <class X, class D>
        concept
    shiftable_by = requires (X const& x, int k, D const& delta)
    {
        { x + k * delta } -> std::same_as <X>;
    };

        template <std::size_t Index, class... Args>
        using
    argument_t = std::remove_cvref_t <std::tuple_element_t <Index, std::tuple <Args...>>>;
} // namespace detail

    template <std::size_t I> struct S;

    template <
          std::size_t Index = 0
        , std::size_t Order = 1
        , std::size_t Accuracy = 2
    >
    constexpr const auto
central_finite_difference (auto&& F, auto&& delta, auto&&... args)
    requires (Index < sizeof... (args))
        and detail::shiftable_by <
              detail::argument_t <Index, decltype (args)...>
            , std::remove_cvref_t <decltype (delta)>
        >
{
    static_assert (Order > 0, "The order of derivation must be at least one");
    static_assert (Accuracy > 1, "The accuracy of derivation must be at least two");
    static_assert (Accuracy % 2 == 0, "The accuracy of derivation must be even");
    static_assert (Index < sizeof... (args), "Invalid index");

    return detail::evaluate <Index, Order, Accuracy> (F, delta, args...);
}
} // namespace calculisto::finite_difference
