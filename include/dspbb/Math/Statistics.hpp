#pragma once

#include <cmath>
#include <dspbb/Math/DotProduct.hpp>
#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/Primitives/SignalView.hpp>
#include <dspbb/Vectorization/Kernels.hpp>
#include <dspbb/Vectorization/MathFunctions.hpp>
#include <numeric>

namespace dspbb {


//------------------------------------------------------------------------------
// General stat functions
//------------------------------------------------------------------------------
template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto Sum(const SignalT& signal) {
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	return ReduceVectorized(signal.Data(), signal.Size(), T(0), [](const auto& a, const auto& b) { return a + b; });
}


template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto Mean(const SignalT& signal) {
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	return !signal.Empty() ? Sum(signal) / T(signal.Size()) : T(0);
}


template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto SumSquare(const SignalT& signal) {
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	return MapReduceVectorized(
		signal.Data(),
		signal.Size(),
		T(0),
		[](const auto& a, const auto& b) { return a + b; },
		[](const auto& a) { return a * a; });
}


template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto MeanSquare(const SignalT& signal) {
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	return !signal.Empty() ? SumSquare(signal) / T(signal.Size()) : T(0);
}


template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto RootMeanSquare(const SignalT& signal) {
	return std::sqrt(MeanSquare(signal));
}

template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto Norm(const SignalT& signal) {
	return std::sqrt(SumSquare(signal));
}


template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto Max(const SignalT& signal) {
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	assert(!signal.Empty());
	return ReduceVectorized(signal.Data(), signal.Size(), signal[0], [](const auto& a, const auto& b) { return math_functions::max(a, b); });
}


template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto Min(const SignalT& signal) {
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	assert(!signal.Empty());
	return ReduceVectorized(signal.Data(), signal.Size(), signal[0], [](const auto& a, const auto& b) { return math_functions::min(a, b); });
}


//------------------------------------------------------------------------------
// Moments
//------------------------------------------------------------------------------
template <class SignalT, class U, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto CentralMoment(const SignalT& signal, size_t k, U mean) {
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	const auto add = [](const auto& a, const auto& b) { return a + b; };
	const auto m2 = [](const auto& a, auto mean) { const auto d = a - mean; return d*d; };
	const auto m3 = [](const auto& a, auto mean) { const auto d = a - mean; return d*d*d; };
	const auto m4 = [](const auto& a, auto mean) { const auto d = a - mean; return d*d*d*d; };
	const auto mg = [](const auto& a, auto mean, auto k) { const auto d = a - mean; return d * math_functions::pow(math_functions::abs(d), decltype(d)(k-1)); };

	const T den = T(signal.Size());

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
auto CentralMoment(const SignalT& signal, size_t k) {
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	const T mean = ReduceVectorized(signal.Data(), signal.Size(), T(0), add) / den;
	return CentralMoment(signal, k, mean);
}

template <class SignalT, class U, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto StandardizedMoment(const SignalT& signal, size_t k, U mean) {
	const auto variance = CentralMoment(signal, 2, mean);
	return CentralMoment(signal, k, mean) / std::pow(variance, decltype(variance)(k) / decltype(variance)(2));
}

template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto StandardizedMoment(const SignalT& signal, size_t k) {
	return CentralMoment(signal, k, Mean(signal));
}


//------------------------------------------------------------------------------
// Moments with special name
//------------------------------------------------------------------------------

template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto StandardDeviation(const SignalT& signal) {
	return std::sqrt(CentralMoment(signal, 2));
}

template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto Variance(const SignalT& signal) {
	return CentralMoment(signal, 2);
}

template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto Skewness(const SignalT& signal) {
	return StandardizedMoment(signal, 3);
}

template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto Kurtosis(const SignalT& signal) {
	return StandardizedMoment(signal, 4);
}

template <class SignalT, class U, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto StandardDeviation(const SignalT& signal, U mean) {
	return std::sqrt(CentralMoment(signal, 2, mean));
}

template <class SignalT, class U, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto Variance(const SignalT& signal, U mean) {
	return CentralMoment(signal, 2, mean);
}

template <class SignalT, class U, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto Skewness(const SignalT& signal, U mean) {
	return StandardizedMoment(signal, 3, mean);
}

template <class SignalT, class U, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto Kurtosis(const SignalT& signal, U mean) {
	return StandardizedMoment(signal, 4, mean);
}

template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto CorrectedStandardDeviation(const SignalT& signal) {
	using T = typename SignalT::value_type;
	return std::sqrt(CentralMoment(signal, 2) * T(signal.Size()) / T(signal.Size() - 1));
}

template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto CorrectedVariance(const SignalT& signal) {
	return CentralMoment(signal, 2) * T(signal.Size()) / T(signal.Size() - 1);
}

template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto CorrectedSkewness(const SignalT& signal) {
	//return StandardizedMoment(signal, 3);
	throw std::logic_error("not implemented");
}

template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto CorrectedKurtosis(const SignalT& signal) {
	//return StandardizedMoment(signal, 4);
	throw std::logic_error("not implemented");
}



//------------------------------------------------------------------------------
// Covariance & correlation
//------------------------------------------------------------------------------

template <class SignalT, class SignalU, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>> && is_signal_like_v<std::decay_t<SignalU>>, int> = 0>
auto Covariance(const SignalT& a, const SignalU& b) {
	assert(a.Size() == b.Size());
	const auto size = a.Size();
	using R = decltype(std::declval<typename SignalT::value_type>() * std::declval<typename SignalU::value_type>());
	const auto ma = Mean(a);
	const auto mb = Mean(b);
	return InnerProduct(
			   a.Data(),
			   b.Data(),
			   size,
			   R(0),
			   [ma, mb](const auto& a, const auto& b) { return (a - ma) * (b - mb); },
			   [](const auto& acc, const auto& x) { return acc + x; })
		   / R(size);
}


template <class SignalT, class SignalU, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>> && is_signal_like_v<std::decay_t<SignalU>>, int> = 0>
auto Correlation(const SignalT& a, const SignalU& b) {
	assert(a.Size() == b.Size());
	return Covariance(a, b) / (StandardDeviation(a) * StandardDeviation(b));
}


} // namespace dspbb