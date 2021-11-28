#pragma once

#include <complex>

namespace dspbb {

template <class T>
constexpr T pi_v = T(3.1415926535897932384626433);

template <class T>
constexpr std::complex<T> i_v(T(0), T(1));

}