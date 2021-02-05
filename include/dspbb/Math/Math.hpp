#pragma once

#include "../Primitives/Signal.hpp"

#include <complex>


namespace dspbb {

template <class T, eSignalDomain Domain>
Signal<T, Domain> Abs(SignalView<const T, Domain> signal) {
	Signal<T, Domain> absval = signal;
	for (auto& v : absval) {
		v = std::abs(v);
	}
	return absval;
}


template <class T, eSignalDomain Domain>
Signal<T, Domain> Abs(SignalView<const std::complex<T>, Domain> signal) {
	Signal<T, Domain> absval(signal.Length());
	for (size_t i = 0; i < signal.Length(); ++i) {
		absval[i] = std::abs(signal[i]);
	}
	return absval;
}


template <class T, eSignalDomain Domain>
Signal<T, Domain> Log(SignalView<const T, Domain> signal) {
	Signal<T, Domain> ret = signal;
	for (auto& v : ret) {
		v = std::log(v);
	}
	return ret;
}

template <class T, eSignalDomain Domain>
const Signal<T, Domain> Real(SignalView<const T, Domain> signal) {
	return { signal.begin(), signal.end() };
}


template <class T, eSignalDomain Domain>
Signal<T, Domain> Real(SignalView<const std::complex<T>, Domain> signal) {
	return { signal.Real(), signal.Real() + signal.Size() };
}


template <class T, eSignalDomain Domain>
Signal<T, Domain> Imag(SignalView<const std::complex<T>, Domain> signal) {
	return { signal.Imag(), signal.Imag() + signal.Size() };
}



// Wrappers

template <class T, eSignalDomain Domain>
auto Abs(const Signal<std::complex<T>, Domain>& signal) { return Abs(AsConstSpan(signal)); }

template <class T, eSignalDomain Domain>
auto Log(const Signal<std::complex<T>, Domain>& signal) { return Log(AsConstSpan(signal)); }

template <class T, eSignalDomain Domain>
auto Real(const Signal<T, Domain>& signal) { return Real(AsConstSpan(signal)); }

template <class T, eSignalDomain Domain>
auto Imag(const Signal<std::complex<T>, Domain>& signal) { return Imag(AsConstSpan(signal)); }



} // namespace dspbb