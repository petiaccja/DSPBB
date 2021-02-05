#pragma once

#include "../Primitives/Signal.hpp"

#include <complex>


namespace enl {

template <class T, eSignalDomain Domain>
Signal<T, Domain> abs(const Signal<T, Domain>& signal) {
	Signal<T, Domain> absval = signal;
	for (auto& v : absval) {
		v = std::abs(v);
	}
	return absval;
}


template <class T, eSignalDomain Domain>
Signal<T, Domain> abs(const Signal<std::complex<T>, Domain>& signal) {
	Signal<T, Domain> absval(signal.Length());
	for (size_t i = 0; i < signal.Length(); ++i) {
		absval[i] = std::abs(signal[i]);
	}
	return absval;
}


template <class T, eSignalDomain Domain>
Signal<T, Domain> log(const Signal<T, Domain>& signal) {
	Signal<T, Domain> ret = signal;
	for (auto& v : ret) {
		v = std::log(v);
	}
	return ret;
}

template <class T, eSignalDomain Domain>
const Signal<T, Domain>& real(const Signal<T, Domain>& signal) {
	return signal;
}


template <class T, eSignalDomain Domain>
Signal<T, Domain> real(const Signal<std::complex<T>, Domain>& signal) {
	return { signal.Real(), signal.Real() + signal.Size() };
}


template <class T, eSignalDomain Domain>
Signal<T, Domain> imag(const Signal<std::complex<T>, Domain>& signal) {
	return { signal.Imag(), signal.Imag() + signal.Size() };
}








} // namespace enl