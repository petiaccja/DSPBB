#pragma once

#include <cmath>
#include <dspbb/Math/DotProduct.hpp>
#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/Primitives/SignalView.hpp>
#include <dspbb/Vectorization/MathFunctions.hpp>
#include <dspbb/Vectorization/Kernels.hpp>
#include <numeric>

namespace dspbb {

template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto CentralMoment(const SignalT& signal, size_t k) {
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	const auto add = [](const auto& a, const auto& b) { return a + b; };
	const auto m2 = [](const auto& a, auto mean) { const auto d = a - mean; return d*d; };
	const auto m3 = [](const auto& a, auto mean) { const auto d = a - mean; return d*d*d; };
	const auto m4 = [](const auto& a, auto mean) { const auto d = a - mean; return d*d*d*d; };
	const auto mg = [](const auto& a, auto mean, auto k) { const auto d = a - mean; return d * math_functions::pow(math_functions::abs(d), decltype(d)(k-1)); };

	const T den = T(signal.Size());
	const T mean = Reduce(signal.Data(), signal.Size(), T(0), add) / den;

	switch (k) {
		case 0: return T(0);
		case 1: return T(0);
		case 2: return MapReduceVectorized(signal.Data(), signal.Size(), T(0), add, [mean, m2](const auto& a) { return m2(a, mean); }) / den;
		case 3: return MapReduceVectorized(signal.Data(), signal.Size(), T(0), add, [mean, m3](const auto& a) { return m3(a, mean); }) / den;
		case 4: return MapReduceVectorized(signal.Data(), signal.Size(), T(0), add, [mean, m4](const auto& a) { return m4(a, mean); }) / den;
		default: return MapReduceVectorized(signal.Data(), signal.Size(), T(0), add, [mean, mg, k](const auto& a) { return mg(a, mean, k); }) / den;
	}
}

template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto StandardizedMoment(const SignalT& signal, size_t k) {
	const auto variance = CentralMoment(signal, 2);
	return CentralMoment(signal, k) / std::pow(variance, decltype(variance)(k) / decltype(variance)(2));
}

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