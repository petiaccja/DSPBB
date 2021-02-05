#include "dspbb/Math/DotProduct.hpp"
#include "dspbb/Primitives/Signal.hpp"
#include "dspbb/Primitives/SignalView.hpp"

#include <cmath>
#include <numeric>

namespace dspbb {

template <class T, eSignalDomain Domain>
T Sum(SignalView<const T, Domain> signal) {
	return std::accumulate(signal.begin(), signal.end(), T(0));
}


template <class T, eSignalDomain Domain>
T Mean(SignalView<const T, Domain> signal) {
	return !signal.Empty() ? Sum(signal) / T(signal.Size()) : T(0);
}


template <class T, eSignalDomain Domain>
T SumSquare(SignalView<const T, Domain> signal) {
	return DotProduct(signal, signal, signal.Size());
}


template <class T, eSignalDomain Domain>
T RootMeanSquare(SignalView<const T, Domain> signal) {
	return !signal.Empty() ? std::sqrt(Sum(signal) / T(signal.Size())) : T(0);
}


template <class T, eSignalDomain Domain>
T StandardDeviation(SignalView<const T, Domain> signal) {
	const T mean = Mean(signal);
	const Signal<T, Domain> diff = Signal<T, Domain>{ signal.begin(), signal.end() } - mean;
	T sum = SumSquare(diff);
	return signal.Size() >= 2 ? std::sqrt(sum / T(signal.Size() - 1)) : T(0);
}


template <class T, eSignalDomain Domain>
T Norm(SignalView<const T, Domain> signal) {
	return std::sqrt(SumSquare(signal));
}


// Wrappers.

template <class T, eSignalDomain Domain>
T Sum(const Signal<T, Domain>& signal) { return Sum(AsConstSpan(signal)); }

template <class T, eSignalDomain Domain>
T Mean(const Signal<T, Domain>& signal) { return Mean(AsConstSpan(signal)); }

template <class T, eSignalDomain Domain>
T SumSquare(const Signal<T, Domain>& signal) { return SumSquare(AsConstSpan(signal)); }

template <class T, eSignalDomain Domain>
T RootMeanSquare(const Signal<T, Domain>& signal) { return RootMeanSquare(AsConstSpan(signal)); }

template <class T, eSignalDomain Domain>
T StandardDeviation(const Signal<T, Domain>& signal) { return StandardDeviation(AsConstSpan(signal)); }

template <class T, eSignalDomain Domain>
T Norm(const Signal<T, Domain>& signal) { return Norm(AsConstSpan(signal)); }

} // namespace dspbb