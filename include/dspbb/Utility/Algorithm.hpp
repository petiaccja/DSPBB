#pragma once
#include "dspbb/Primitives/Signal.hpp"
#include "dspbb/Primitives/SignalView.hpp"

namespace dspbb {

//------------------------------------------------------------------------------
// Apply function to signal.
//------------------------------------------------------------------------------

template <class T, eSignalDomain Domain, class Func, class... Args, std::enable_if_t<!std::is_const_v<T> && std::is_same_v<std::decay_t<std::invoke_result_t<Func, T, Args...>>, T>, int> = 0>
void Apply(SignalView<T, Domain> inout, Func func, const Args&... args) {
	std::for_each(inout.begin(), inout.end(), [&func, &args...](T& item) { item = func(item, args...); });
}

template <class T, dspbb::eSignalDomain Domain, class Func, class... Args, std::enable_if_t<std::is_same_v<std::decay_t<std::invoke_result_t<Func, T, Args...>>, T>, int> = 0>
Signal<T, Domain> Apply(Signal<T, Domain>&& inout, Func func, const Args&... args) {
	Apply(AsView(inout), func, args...);
	return std::move(inout);
}

template <class T, eSignalDomain Domain, class Func, class... Args, std::enable_if_t<std::is_same_v<std::decay_t<std::invoke_result_t<Func, T, Args...>>, T>, int> = 0>
Signal<T, Domain> Apply(SignalView<const T, Domain> view, Func func, const Args&... args) {
	return Apply(Signal<T, Domain>{ view.begin(), view.end() }, func, args...);
}

template <class T, eSignalDomain Domain, class Func, class... Args, std::enable_if_t<std::is_same_v<std::decay_t<std::invoke_result_t<Func, T, Args...>>, T>, int> = 0>
Signal<T, Domain> Apply(const Signal<T, Domain>& signal, Func func, const Args&... args) {
	return Apply(Signal<T, Domain>{ signal }, func, args...);
}

template <class T, eSignalDomain Domain, class Func, class... Args, std::enable_if_t<!std::is_same_v<std::decay_t<std::invoke_result_t<Func, T, Args...>>, T>, int> = 0>
auto Apply(SignalView<const T, Domain> view, Func func, const Args&... args) {
	using R = std::decay_t<std::invoke_result_t<Func, T, Args...>>;
	Signal<R, Domain> r(view.Size());
	for (size_t i = 0; i < view.Size(); ++i) {
		r[i] = func(view[i], args...);
	}
	return r;
}

template <class T, eSignalDomain Domain, class Func, class... Args, std::enable_if_t<!std::is_same_v<std::decay_t<std::invoke_result_t<Func, T, Args...>>, T>, int> = 0>
auto Apply(const Signal<T, Domain>& signal, Func func, const Args&... args) {
	return Apply(AsConstView(signal), func, args...);
}

} // namespace dspbb