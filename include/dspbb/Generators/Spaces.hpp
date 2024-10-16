#pragma once

#include "../Math/Functions.hpp"
#include "../Signal/Signal.hpp"
#include "../Signal/Traits.hpp"
#include "../Utility/TypeTraits.hpp"

#include <type_traits>


namespace dspbb {

/// <summary> Generate a signal that linearly interpolates between two values. </summary>
/// <param name="output"> The generated signal is written here. </param>
/// <param name="start"> The value at the start of the signal. </param>
/// <param name="end"> The value at the end of the signal. </param>
/// <param name="inclusive"> True if the signal is a closed interval [start, end], false if [start, end). </param>
template <mutable_signal_or_view_r SignalOut, class Real = remove_complex_t<typename std::decay_t<SignalOut>::value_type>>
void LinSpace(SignalOut&& output,
			  Real start,
			  Real end,
			  bool inclusive = true) {
	const auto count = output.size();
	for (size_t i = 0; i < count; ++i) {
		output[i] = Real(i);
	}
	const Real scale = (end - start) / Real(std::max(intptr_t(1), intptr_t(count) - intptr_t(inclusive)));
	const Real offset = start;
	output *= scale;
	output += offset;
}


/// <summary> Generate a signal that linearly interpolates between two values. </summary>
/// <param name="start"> The value at the start of the signal. </param>
/// <param name="end"> The value at the end of the signal. </param>
/// <param name="count"> The number of samples in the output. </param>
/// <param name="inclusive"> True if the signal is a closed interval [start, end], false if [start, end). </param>
/// <returns> The linear signal. </returns>
template <class T, eSignalDomain Domain>
auto LinSpace(remove_complex_t<T> start,
			  remove_complex_t<T> end,
			  size_t count,
			  bool inclusive = true) {
	BasicSignal<T, Domain> s(count);
	LinSpace(s, start, end, inclusive);
	return s;
}


/// <summary> Generate a signal that logarithmically interpolates between two values. </summary>
/// <param name="output"> The generated signal is written here. </param>
/// <param name="start"> The value at the start of the signal. </param>
/// <param name="end"> The value at the end of the signal. </param>
/// <param name="base"> The base of the logarithm. </param>
/// <param name="inclusive"> True if the signal is a closed interval [start, end], false if [start, end). </param>
template <mutable_signal_or_view_r SignalOut, class Real = remove_complex_t<typename std::decay_t<SignalOut>::value_type>>
void LogSpace(SignalOut&& output,
			  Real start,
			  Real end,
			  Real base = Real(10),
			  bool inclusive = true) {
	LinSpace(output, start, end, inclusive);
	output *= std::log(base);
	Exp(output, output);
}


/// <summary> Generate a signal that logarithmically interpolates between two values. </summary>
/// <param name="start"> The value at the start of the signal. </param>
/// <param name="end"> The value at the end of the signal. </param>
/// <param name="count"> The number of samples in the output. </param>
/// <param name="base"> The base of the logarithm. </param>
/// <param name="inclusive"> True if the signal is a closed interval [start, end], false if [start, end). </param>
/// <returns> The logarithmic signal. </returns>
template <class T, eSignalDomain Domain>
auto LogSpace(remove_complex_t<T> start,
			  remove_complex_t<T> end,
			  size_t count,
			  remove_complex_t<T> base = remove_complex_t<T>(10),
			  bool inclusive = true) {
	BasicSignal<T, Domain> s(count);
	LogSpace(s, start, end, base, inclusive);
	return s;
}

} // namespace dspbb