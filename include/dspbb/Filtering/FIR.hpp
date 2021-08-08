#pragma once

#include "../Math/Functions.hpp"
#include "../Math/Statistics.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "../Utility/Numbers.hpp"
#include "FFT.hpp"
#include "WindowFunctions.hpp"

#include "FIR2.hpp"

namespace dspbb {

template <class T>
struct FirGeneralQuality {
	T cosineSimilarity;
	T maxDifference;
	T maxRelDifference;
	T meanDifference;
	T meanRelDifference;
	T meanSquareDifference;
	T meanSquareRelDifference;
};


/// <summary>
/// Measures how well the <paramref name="filter"/> approximates the desired frequency response.
/// </summary>
/// <param name="filter"> The impulse response of the filter. </param>
/// <param name="desiredResponse"> The desired frequency response of the filter to which the real filter is compared against.
///		Bins range from 0Hz to sampleRate/2, and the number of bins must be at least half of the number of filter taps. </param>
/// <returns> A number between 0 (extremely bad match) and 1 (perfect match). </returns>
template <class T>
FirGeneralQuality<T> FirQuality(const TimeSignal<T>& filter, const Spectrum<T>& desiredResponse) {
	const size_t numBins = 2 * desiredResponse.Length() - 1 - desiredResponse.Length() % 2;
	assert(numBins >= filter.Size());
	TimeSignal<T> extendedFilter = filter;
	extendedFilter.Resize(numBins, T(0));
	const auto actualResponse = Abs(FourierTransform(extendedFilter, false));
	assert(actualResponse.Size() == desiredResponse.Size());
	const auto difference = actualResponse - desiredResponse;
	const auto relDifference = difference / desiredResponse;

	auto magActual = DotProduct(AsConstView(actualResponse), AsConstView(actualResponse), desiredResponse.Size());
	auto magDesired = DotProduct(AsConstView(desiredResponse), AsConstView(desiredResponse), desiredResponse.Size());
	auto similarity = DotProduct(AsConstView(desiredResponse), AsConstView(actualResponse), desiredResponse.Size());

	FirGeneralQuality<T> quality;
	quality.cosineSimilarity = similarity / std::max(magActual, magDesired);
	quality.maxDifference = Max(Abs(difference));
	quality.maxRelDifference = Max(Abs(relDifference));
	quality.meanDifference = Mean(Abs(difference));
	quality.meanRelDifference = Mean(Abs(relDifference));
	quality.meanSquareDifference = RootMeanSquare(difference);
	quality.meanSquareRelDifference = RootMeanSquare(relDifference);

	return quality;
}


} // namespace dspbb