#pragma once

#include "Definitions.hpp"

#include <type_traits>


namespace dspbb {

//------------------------------------------------------------------------------
// Determine if object is a signal.
//------------------------------------------------------------------------------

template <class>
struct is_signal : std::false_type {};

template <class T, eSignalDomain Domain>
struct is_signal<BasicSignal<T, Domain>> : std::true_type {};

template <class T>
inline constexpr bool is_signal_v = is_signal<T>::value;


template <class>
struct is_signal_view : std::false_type {};

template <class T, eSignalDomain Domain>
struct is_signal_view<BasicSignalView<T, Domain>> : std::true_type {};

template <class T>
constexpr bool is_signal_view_v = is_signal_view<T>::value;


template <class T>
struct is_signal_like {
	static constexpr bool value = is_signal<T>::value || is_signal_view<T>::value;
};

template <class T>
inline constexpr bool is_signal_like_v = is_signal_like<T>::value;

template <class T>
concept signal_like = is_signal_like_v<T>;


//------------------------------------------------------------------------------
// Signal traits.
//------------------------------------------------------------------------------

template <class>
struct scalar_type {};

template <class T, eSignalDomain Domain>
struct scalar_type<BasicSignal<T, Domain>> {
	using type = T;
};

template <class T, eSignalDomain Domain>
struct scalar_type<BasicSignalView<T, Domain>> {
	using type = T;
};

template <class T>
using scalar_type_t = typename scalar_type<T>::type;


template <class>
struct domain {};

template <class T, eSignalDomain Domain>
struct domain<BasicSignal<T, Domain>> {
	static constexpr auto value = Domain;
};

template <class T, eSignalDomain Domain>
struct domain<BasicSignalView<T, Domain>> {
	static constexpr auto value = Domain;
};

template <class T>
inline constexpr eSignalDomain domain_v = domain<T>::value;


template <class>
struct is_mutable {};

template <class T, eSignalDomain Domain>
struct is_mutable<BasicSignal<T, Domain>> {
	static constexpr bool value = true;
};

template <class T, eSignalDomain Domain>
struct is_mutable<BasicSignalView<T, Domain>> {
	static constexpr bool value = !std::is_const_v<T>;
};

template <class T>
inline constexpr eSignalDomain is_mutable_v = is_mutable<T>::value;


//------------------------------------------------------------------------------
// Old stuff.
//------------------------------------------------------------------------------


template <class... Signals>
struct is_same_domain {
	static constexpr bool compare() { return true; }
	template <class H1>
	static constexpr bool compare() { return true; }
	template <class H1, class H2, class... Tail>
	static constexpr bool compare() {
		return domain_v<H1> == domain_v<H2> && compare<H2, Tail...>();
	}
	template <class... Signals_, std::enable_if_t<std::conjunction_v<is_signal_like<Signals_>...>, int> = 0>
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
		return !std::is_const_v<std::remove_reference_t<Signal_>>;
	}
	template <class SignalView_, std::enable_if_t<is_signal_view_v<std::decay_t<SignalView_>>, int> = 0>
	static constexpr bool test(int) {
		return !std::is_const_v<scalar_type_t<std::decay_t<SignalView_>>>;
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