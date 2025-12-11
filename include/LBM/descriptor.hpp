#include <concepts>
#include <cstdint>
#include <array>

#ifndef UTILS_HPP
#define UTILS_HPP

template <int D, int Q>
struct LatticeDescriptor {

    static constexpr int d = D;
    static constexpr int q = Q;

};


template <typename T>
concept isDescriptor = requires{

    T::d;
    T::q;

} && std::derived_from<T, LatticeDescriptor<T::d, T::q>>;


template<std::floating_point float_type>
struct D2Q9 : public LatticeDescriptor<2, 9> {

    static constexpr std::array<std::array<int, 2>, 9> c = {{
        {{ 0, 0}},
        {{ 1, 0}}, {{ 0, 1}}, {{-1, 0}}, {{ 0,-1}},
        {{ 1, 1}}, {{-1, 1}}, {{-1,-1}}, {{ 1,-1}}
    }};
    static constexpr std::array<float_type, 9> w = {
        4.0 / 9.0,
        1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
    };
    static constexpr std::array<int, 9> opposite = {0, 3, 4, 1, 2, 7, 8, 5, 6};
};


template<std::floating_point float_type>
struct D3Q19 : public LatticeDescriptor<3, 19> {

    static constexpr std::array<std::array<int, 3>, 19> c = {{
        {{ 0, 0, 0}},
        {{ 1, 0, 0}}, {{ -1, 0, 0}}, {{ 0, 1, 0}}, {{ 0, -1, 0}}, {{ 0, 0, 1}}, {{ 0, 0, -1}},
        {{ 1, 1, 0}}, {{ 1, -1, 0}}, {{ -1, 1, 0}}, {{ -1, -1, 0}},
        {{ 1, 0, 1}}, {{ 1, 0, -1}}, {{ -1, 0, 1}}, {{ -1, 0, -1}},
        {{ 0, 1, 1}}, {{ 0, 1, -1}}, {{ 0, -1, 1}}, {{ 0, -1, -1}}
    }};
    static constexpr std::array<float_type, 19> w = {
        1.0 / 3.0,
        1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
    };
    static constexpr std::array<int, 19> opposite = {0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15};
};

template <int D>
struct Direction;

template <>
struct Direction<3> {
    enum class Value : unsigned char {
        NODIR   = 0,
        LEFT  = 1,
        RIGHT = 2,
        BOTTOM = 3,
        TOP   = 4,
        FRONT  = 5,
        BACK = 6
    };
    using enum Value;
};

template <>
struct Direction<2> {
    enum class Value : unsigned char {
        NODIR = 0,
        LEFT  = 1,
        RIGHT = 2,
        BOTTOM= 3,
        TOP   = 4
    };
    using enum Value;
};

#endif
