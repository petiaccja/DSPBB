#pragma once

#include <dspbb/Math/DotProduct.hpp>
#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/Primitives/SignalView.hpp>

#include <cmath>
#include <numeric>

namespace dspbb {

template <class T, eSignalDomain Domain>
T Sum(SignalView<T, Domain> signal) {
	return std::accumulate(signal.begin(), signal.end(), T(0));
}


template <class T, eSignalDomain Domain>
T Mean(SignalView<T, Domain> signal) {
	return !signal.Empty() ? Sum(signal) / T(signal.Size()) : T(0);
}


template <class T, eSignalDomain Domain>
T SumSquare(SignalView<T, Domain> signal) {
	return DotProduct(signal, signal, signal.Size());
}


template <class T, eSignalDomain Domain>
T MeanSquare(SignalView<T, Domain> signal) {
	return !signal.Empty() ? SumSquare(signal) / T(signal.Size()) : T(0);
}


template <class T, eSignalDomain Domain>
T RootMeanSquare(SignalView<T, Domain> signal) {
	return std::sqrt(MeanSquare(signal));
}


template <class T, eSignalDomain Domain>
T StandardDeviation(SignalView<T, Domain> signal) {
	const T mean = Mean(signal);
	const Signal<std::remove_cv_t<T>, Domain> diff = signal - mean;
	T sum = SumSquare(diff);
	return signal.Size() >= 1 ? std::sqrt(sum / T(signal.Size())) : T(0);
}


template <class T, eSignalDomain Domain>
T Norm(SignalView<T, Domain> signal) {
	return std::sqrt(SumSquare(signal));
}


template <class T, eSignalDomain Domain>
T Max(SignalView<T, Domain> signal) {
	assert(!signal.Empty());
	return *std::max_element(signal.begin(), signal.end());
}


template <class T, eSignalDomain Domain>
T Min(SignalView<T, Domain> signal) {
	assert(!signal.Empty());
	return *std::min_element(signal.begin(), signal.end());
}


// Wrappers.

template <class T, eSignalDomain Domain>
T Sum(const Signal<T, Domain>& signal) { return Sum(AsConstView(signal)); }

template <class T, eSignalDomain Domain>
T Mean(const Signal<T, Domain>& signal) { return Mean(AsConstView(signal)); }

template <class T, eSignalDomain Domain>
T SumSquare(const Signal<T, Domain>& signal) { return SumSquare(AsConstView(signal)); }

template <class T, eSignalDomain Domain>
T MeanSquare(const Signal<T, Domain>& signal) { return MeanSquare(AsConstView(signal)); }

template <class T, eSignalDomain Domain>
T RootMeanSquare(const Signal<T, Domain>& signal) { return RootMeanSquare(AsConstView(signal)); }

template <class T, eSignalDomain Domain>
T StandardDeviation(const Signal<T, Domain>& signal) { return StandardDeviation(AsConstView(signal)); }

template <class T, eSignalDomain Domain>
T Norm(const Signal<T, Domain>& signal) { return Norm(AsConstView(signal)); }

template <class T, eSignalDomain Domain>
T Max(const Signal<T, Domain>& signal) { return Max(AsConstView(signal)); }

template <class T, eSignalDomain Domain>
T Min(const Signal<T, Domain>& signal) { return Min(AsConstView(signal)); }

} // namespace dspbb