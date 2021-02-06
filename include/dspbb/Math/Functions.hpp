#pragma once

#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"

#include <complex>


namespace dspbb {

template <class T, eSignalDomain Domain>
Signal<T, Domain> Abs(SignalView<const T, Domain> signal) {
	Signal<T, Domain> absval{ signal.begin(), signal.end() };
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
	Signal<T, Domain> ret{ signal.begin(), signal.end() };
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
	Signal<T, Domain> real(signal.Size());
	for (size_t i = 0; i < signal.Size(); ++i) {
		real[i] = signal[i].real();
	}
	return real;
}

template <class T, eSignalDomain Domain>
Signal<T, Domain> Imag(SignalView<const std::complex<T>, Domain> signal) {
	Signal<T, Domain> imag(signal.Size());
	for (size_t i = 0; i < signal.Size(); ++i) {
		imag[i] = signal[i].imag();
	}
	return imag;
}



// Wrappers

template <class T, eSignalDomain Domain>
auto Abs(const Signal<T, Domain>& signal) { return Abs(AsConstView(signal)); }

template <class T, eSignalDomain Domain>
auto Log(const Signal<T, Domain>& signal) { return Log(AsConstView(signal)); }

template <class T, eSignalDomain Domain>
auto Real(const Signal<T, Domain>& signal) { return Real(AsConstView(signal)); }

template <class T, eSignalDomain Domain>
auto Imag(const Signal<std::complex<T>, Domain>& signal) { return Imag(AsConstView(signal)); }



} // namespace dspbb