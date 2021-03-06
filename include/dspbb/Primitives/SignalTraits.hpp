#pragma once

#include "Signal.hpp"
#include "SignalView.hpp"

namespace dspbb {

enum class eSignalDomain;
template <class T, eSignalDomain Domain>
class Signal;
template <class T, eSignalDomain Domain>
class SignalView;

}

namespace dspbb {


template <class T>
struct is_signal : std::false_type {};

template <class T, eSignalDomain Domain>
struct is_signal<Signal<T, Domain>> : std::true_type {};

template <class T>
constexpr bool is_signal_v = is_signal<T>::value;

template <class T>
struct is_signal_view : std::false_type {};

template <class T, eSignalDomain Domain>
struct is_signal_view<SignalView<T, Domain>> : std::true_type {};

template <class T>
constexpr bool is_signal_view_v = is_signal_view<T>::value;

template <class T>
struct is_signal_like {
	static constexpr bool value = is_signal<T>::value || is_signal_view<T>::value;
};

template <class T>
constexpr bool is_signal_like_v = is_signal_like<T>::value;

template <class SignalT>
struct signal_traits;

template <class T, eSignalDomain Domain>
struct signal_traits<Signal<T, Domain>> {
	using type = T;
	static constexpr auto domain = Domain;
};

template <class T, eSignalDomain Domain>
struct signal_traits<SignalView<T, Domain>> {
	using type = T;
	static constexpr auto domain = Domain;
};


template <class...>
struct conjunction_compat : std::true_type {};
template <class B1>
struct conjunction_compat<B1> : B1 {};
template <class B1, class... Bn>
struct conjunction_compat<B1, Bn...>
	: std::conditional_t<bool(B1::value), conjunction_compat<Bn...>, B1> {};


template <class... Signals>
struct is_same_domain {
	static constexpr bool compare() { return true; }
	template <class H1>
	static constexpr bool compare() { return true; }
	template <class H1, class H2, class... Tail>
	static constexpr bool compare() {
		return signal_traits<H1>::domain == signal_traits<H2>::domain && compare<H2, Tail...>();
	}
	template <class... Signals_, std::enable_if_t<conjunction_compat<is_signal_like<Signals_>...>::value, int> = 0>
	static constexpr bool test(int) {
		return compare<Signals_...>();
	}
	template <class... Signals_>
	static constexpr bool test(...) {
		return false;
	}
	static constexpr bool value = test<std::decay_t<Signals>...>(0);
};

template <class... Signals>
constexpr bool is_same_domain_v = is_same_domain<Signals...>::value;

template <class Signal>
struct is_mutable_signal {
	template <class Signal_, std::enable_if_t<is_signal_v<std::decay_t<Signal_>>, int> = 0>
	static constexpr bool test(int) {
		return !std::is_const<std::remove_reference_t<Signal_>>::value;
	}
	template <class SignalView_, std::enable_if_t<is_signal_view_v<std::decay_t<SignalView_>>, int> = 0>
	static constexpr bool test(int) {
		return !std::is_const<typename signal_traits<std::decay_t<SignalView_>>::type>::value;
	}
	template <class Signal_>
	static constexpr bool test(...) {
		return false;
	}
	static constexpr bool value = test<Signal>(0);
};

template <class Signal>
constexpr bool is_mutable_signal_v = is_mutable_signal<Signal>::value;



} // namespace dspbb