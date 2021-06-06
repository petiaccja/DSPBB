#pragma once

#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalTraits.hpp"

#include <numeric>


namespace dspbb {

template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
auto LinSpace(SignalR&& output,
			  remove_complex_t<typename signal_traits<std::decay_t<SignalR>>::type> start,
			  remove_complex_t<typename signal_traits<std::decay_t<SignalR>>::type> end,
			  bool inclusive = true) {
	using R = remove_complex_t<typename signal_traits<std::decay_t<SignalR>>::type>;
	const auto count = output.Size();
	std::iota(output.begin(), output.end(), size_t(0));
	const R scale = (end - start) / R(std::max(intptr_t(1), intptr_t(count) - intptr_t(inclusive)));
	const R offset = start;
	output *= scale;
	output += offset;
}

template <class T, eSignalDomain Domain>
auto LinSpace(remove_complex_t<T> start, remove_complex_t<T> end, size_t count, bool inclusive = true) {
	Signal<T, Domain> s(count);
	LinSpace(s, start, end, inclusive);
	return s;
}


} // namespace dspbb