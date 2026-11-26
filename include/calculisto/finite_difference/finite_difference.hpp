#pragma once
#include <tuple>

    namespace
calculisto::finite_difference
{
    namespace
detail
{
        template <std::size_t Index, class T, class... Args>
        constexpr auto
    shift_by (T&& delta, Args&&... args)
    {
            auto
        r = std::tuple { std::forward <Args> (args)... };
        std::get <Index> (r) += delta;
        return r;
    }
} // namespace detail

    template <std::size_t I> struct S;

    template <
          std::size_t Index = 0
        , std::size_t Order = 1
        , std::size_t Accuracy = 2
    >
    constexpr const auto
central_finite_difference (auto&& F, auto&& delta, auto&&... args)
{
    static_assert (Order > 0, "The order of derivation must be at least one");
    static_assert (Accuracy > 1, "The accuracy of derivation must be at least two");
    static_assert (Accuracy % 2 == 0, "The accuracy of derivation must be even");
    static_assert (Index < sizeof... (args), "Invalid index");

        using
    detail::shift_by;
    if constexpr (Order == 1)
    {
        if constexpr (Accuracy == 2)
        {
                const auto
            xm1 = shift_by <Index> (-1 * delta, args...);
                auto
            xp1 = shift_by <Index> (+1 * delta, args...);
                const auto
            r = (
                - std::apply (F, xm1)
                + std::apply (F, xp1)
            ) / 2 / delta;
            return r;
        }
        else if constexpr (Accuracy == 4)
        {
                const auto
            xm1 = shift_by <Index> (-1 * delta, args...);
                const auto
            xm2 = shift_by <Index> (-2 * delta, args...);
                auto
            xp1 = shift_by <Index> (+1 * delta, args...);
                auto
            xp2 = shift_by <Index> (+2 * delta, args...);
                const auto
            r = (
                  1 * std::apply (F, xm2)
                - 8 * std::apply (F, xm1)
                + 8 * std::apply (F, xp1)
                - 1 * std::apply (F, xp2)
            ) / 12 / delta;
            return r;
        }
        else 
        {
            static_assert (false, "Accuracy higher than 4 not implemented");
        }
    }
    else if constexpr (Order == 2)
    {
        if constexpr (Accuracy == 2)
        {
                const auto
            xm1 = shift_by <Index> (-1 * delta, args...);
                const auto
            x   = std::tuple { args... };
                auto
            xp1 = shift_by <Index> (+1 * delta, args...);
                const auto
            r = (
                  1 * std::apply (F, xm1)
                - 2 * std::apply (F, x)
                + 1 * std::apply (F, xp1)
            ) / delta / delta;
            return r;
        }
        else if constexpr (Accuracy == 4)
        {
                const auto
            xm1 = shift_by <Index> (-1 * delta, args...);
                const auto
            xm2 = shift_by <Index> (-2 * delta, args...);
                const auto
            x   = std::tuple { args... };
                auto
            xp1 = shift_by <Index> (+1 * delta, args...);
                auto
            xp2 = shift_by <Index> (+2 * delta, args...);
                const auto
            r = (
                - 1  * std::apply (F, xm2)
                + 16 * std::apply (F, xm1)
                - 30 * std::apply (F, x)
                + 16 * std::apply (F, xp1)
                - 1  * std::apply (F, xp2)
            ) / 12 / delta / delta;
            return r;
        }
        else 
        {
            static_assert (false, "Accuracy higher than 4 not implemented");
        }

    }
    else
    {
        static_assert (false, "Order higher than 2 not implemented");
    }
}
} // namespace calculisto::finite_difference
