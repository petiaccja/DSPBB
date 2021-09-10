#pragma once

#include <type_traits>


namespace dspbb {

struct MethodTagWindowed {};
struct MethodTagLeastSquares {};

constexpr MethodTagWindowed WINDOWED;
constexpr MethodTagLeastSquares LEAST_SQUARES;


//------------------------------------------------------------------------------
// Lowpass desc
//------------------------------------------------------------------------------
template <class Method, class... Args>
struct LowpassDesc;

template <class T, class WindowFunc>
struct LowpassDesc<MethodTagWindowed, T, WindowFunc> {
	LowpassDesc& Cutoff(T value) {
		cutoff = value;
		return *this;
	}
	T cutoff = T(0.5);
	
};


template <class T>
struct LowpassDesc<MethodTagLeastSquares, T> {
	LowpassDesc& Cutoff(T pass, T stop) {
		cutoffPass = pass;
		cutoffStop = stop;
		return *this;
	}
	LowpassDesc& Weight(T pass, T transition, T stop) {
		weightPass = pass;
		weightTransition = transition;
		weightStop = stop;
		return *this;
	}
	LowpassDesc& Smooth(bool enable) {
		smooth = enable;
	}
	
	T cutoffPass = T(0.45);
	T cutoffStop = T(0.55);
	T weightPass = T(1.0);
	T weightTransition = T(0.0);
	T weightStop = T(1.0);
	bool smooth = false;
};



} // namespace dspbb
