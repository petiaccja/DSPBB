#pragma once

#include "../Math/Functions.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalTraits.hpp"
#include "../Utility/TypeTraits.hpp"

#include <numeric>


namespace dspbb {

template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
auto LinSpace(SignalR&& output,
			  remove_complex_t<typename signal_traits<std::decay_t<SignalR>>::type> start,
			  remove_complex_t<typename signal_traits<std::decay_t<SignalR>>::type> end,
			  bool inclusive = true) {
	using R = remove_complex_t<typename signal_traits<std::decay_t<SignalR>>::type>;
	const auto count = output.Size();
	for (size_t i = 0; i < count; ++i) {
		output[i] = R(i);
	}
	const R scale = (end - start) / R(std::max(intptr_t(1), intptr_t(count) - intptr_t(inclusive)));
	const R offset = start;
	output *= scale;
	output += offset;
}

template <class T, eSignalDomain Domain>
auto LinSpace(remove_complex_t<T> start, remove_complex_t<T> end, size_t count, bool inclusive = true) {
	BasicSignal<T, Domain> s(count);
	LinSpace(s, start, end, inclusive);
	return s;
}

template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0, class R = remove_complex_t<typename signal_traits<std::decay_t<SignalR>>::type>>
auto LogSpace(SignalR&& output,
			  R start,
			  R end,
			  R base = R(10),
			  bool inclusive = true) {
	LinSpace(output, start, end, inclusive);
	output *= std::log(base);
	Exp(output, output);
}

template <class T, eSignalDomain Domain>
auto LogSpace(remove_complex_t<T> start, remove_complex_t<T> end, size_t count, remove_complex_t<T> base = remove_complex_t<T>(10), bool inclusive = true) {
	BasicSignal<T, Domain> s(count);
	LogSpace(s, start, end, base, inclusive);
	return s;
}


} // namespace dspbb