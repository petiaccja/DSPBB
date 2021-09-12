#include <array>
#include <catch2/catch.hpp>
#include <dspbb/Filtering/Convolution.hpp>
#include <dspbb/Filtering/FIR.hpp>
#include <dspbb/Filtering/Interpolation.hpp>
#include <dspbb/Generators/Waveforms.hpp>
#include <dspbb/Math/Statistics.hpp>

using namespace dspbb;



//------------------------------------------------------------------------------
// Helpers
//------------------------------------------------------------------------------

TimeSignal<float> GenTestSignal(size_t sampleRate, float frequency, float length = 1.0f) {
	const size_t numSamples = std::max(size_t(1), size_t(double(sampleRate) * double(length)));
	return SineWave<float, TIME_DOMAIN>(numSamples, sampleRate, frequency);
}

template <class SignalT>
bool IsSymmetric(const SignalT& signal) {
	auto beg = signal.begin();
	auto end = signal.rbegin();
	while (beg <= end.base()) {
		if (*beg != Approx(*end).margin(1e-7f)) {
			return false;
		}
		++beg;
		++end;
	}
	return true;
}

template <class SignalT>
bool IsAntiSymmetric(const SignalT& signal) {
	auto beg = signal.begin();
	auto end = signal.rbegin();
	while (beg <= end.base()) {
		if (*beg != Approx(-*end).margin(1e-7f)) {
			return false;
		}
		++beg;
		++end;
	}
	return true;
}

template <class SignalT>
auto MeasureResponse(float frequency, const SignalT& filter) {
	const float period = 1.0f / frequency;
	const float length = std::max(float(filter.Size()), 25.f * period);
	auto testSignal = GenTestSignal(2, frequency, length);
	testSignal *= BlackmanWindow<float, TIME_DOMAIN>(testSignal.Size());
	const auto filteredSignal = Convolution(testSignal, filter, convolution::full);
	const auto rmsTest = std::sqrt(SumSquare(testSignal));
	const auto rmsFiltered = std::sqrt(SumSquare(filteredSignal));
	return rmsFiltered / rmsTest;
}

template <class SignalT>
void RequireResponse(const SignalT& impulse, std::vector<std::pair<float, float>> desired, double margin = 0.03) {
	for (const auto& [frequency, desiredResponse] : desired) {
		const auto actualResponse = MeasureResponse(frequency, impulse);
		REQUIRE(actualResponse == Approx(desiredResponse).margin(margin));
	}
}

constexpr auto TestArbitraryResponse = [](float x) {
	return 2.0f * x - 1.5f * x * x - 0.5f * std::pow(x - 1.f, 3.f);
};


//------------------------------------------------------------------------------
// Window method
//------------------------------------------------------------------------------

TEST_CASE("Windowed low-pass", "[FIR]") {
	constexpr size_t numTaps = 255;
	static constexpr float cutoff = 0.3f;

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Lowpass(WINDOWED).Cutoff(cutoff).Window(windows::blackman));
	REQUIRE(impulse.Size() == numTaps);
	REQUIRE(IsSymmetric(impulse));
	REQUIRE(Sum(impulse) == Approx(1));

	RequireResponse(impulse, { { cutoff - 0.04f, 1.0f }, { cutoff + 0.04f, 0.0f } });
}

TEST_CASE("Windowed high-pass", "[FIR]") {
	constexpr size_t numTaps = 255;
	static constexpr float cutoff = 0.3f;

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Highpass(WINDOWED).Cutoff(cutoff).Window(windows::blackman));
	REQUIRE(impulse.Size() == numTaps);
	REQUIRE(IsSymmetric(impulse));
	REQUIRE(Sum(impulse) == Approx(0).margin(0.01));

	RequireResponse(impulse, { { cutoff - 0.04f, 0.0f }, { cutoff + 0.04f, 1.0f } });
}

TEST_CASE("Windowed band-pass", "[FIR]") {
	constexpr size_t numTaps = 255;
	static constexpr float bandLow = 0.3f;
	static constexpr float bandHigh = 0.6f;

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Bandpass(WINDOWED).Band(bandLow, bandHigh).Window(windows::blackman));
	REQUIRE(impulse.Size() == numTaps);
	REQUIRE(IsSymmetric(impulse));
	REQUIRE(Sum(impulse) == Approx(0).margin(0.01));

	RequireResponse(impulse, {
								 { bandLow - 0.05f, 0.0f },
								 { bandLow + 0.05f, 1.0f },
								 { bandHigh - 0.05f, 1.0f },
								 { bandHigh + 0.05f, 0.0f },
							 });
}

TEST_CASE("Windowed band-stop", "[FIR]") {
	constexpr size_t numTaps = 255;
	static constexpr float bandLow = 0.3f;
	static constexpr float bandHigh = 0.6f;

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Bandstop(WINDOWED).Band(bandLow, bandHigh).Window(windows::blackman));
	REQUIRE(impulse.Size() == numTaps);
	REQUIRE(IsSymmetric(impulse));
	REQUIRE(Sum(impulse) == Approx(1.0).margin(0.01));

	RequireResponse(impulse, {
								 { bandLow - 0.05f, 1.0f },
								 { bandLow + 0.05f, 0.0f },
								 { bandHigh - 0.05f, 0.0f },
								 { bandHigh + 0.05f, 1.0f },
							 });
}

TEST_CASE("Windowed arbitrary", "[FIR]") {
	constexpr size_t numTaps = 255;

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Arbitrary(WINDOWED).Response(TestArbitraryResponse).Window(windows::blackman));
	REQUIRE(impulse.Size() == numTaps);
	REQUIRE(IsSymmetric(impulse));

	RequireResponse(impulse, {
								 { 0.12f, TestArbitraryResponse(0.12f) },
								 { 0.12f, TestArbitraryResponse(0.12f) },
								 { 0.32f, TestArbitraryResponse(0.32f) },
								 { 0.67f, TestArbitraryResponse(0.67f) },
								 { 0.88f, TestArbitraryResponse(0.88f) },
							 });
}

TEST_CASE("Windowed hilbert magnitude", "[FIR]") {
	const auto odd = FirFilter<float, TIME_DOMAIN>(377, Hilbert(WINDOWED).Window(windows::blackman));
	const auto even = FirFilter<float, TIME_DOMAIN>(376, Hilbert(WINDOWED).Window(windows::blackman));

	std::vector<std::pair<float, float>> required = {
		{ 0.1f, 1.0f },
		{ 0.5f, 1.0f },
		{ 0.9f, 1.0f },
	};
	RequireResponse(odd, required);
	RequireResponse(even, required);
}

TEST_CASE("Windowed methods equal", "[FIR]") {
	constexpr size_t numTaps = 127;
	constexpr float cutoff = 0.3f;
	constexpr float bandHigh = 0.2f;
	constexpr float bandLow = 0.6f;

	const auto lp1 = FirFilter<float, TIME_DOMAIN>(numTaps, Lowpass(WINDOWED).Cutoff(cutoff).Window(windows::blackman));
	const auto lp2 = FirFilter<float, TIME_DOMAIN>(numTaps, Lowpass(WINDOWED).Cutoff(cutoff).Window(windows::blackman.operator()<float, TIME_DOMAIN>(numTaps)));

	const auto hp1 = FirFilter<float, TIME_DOMAIN>(numTaps, Highpass(WINDOWED).Cutoff(cutoff).Window(windows::blackman));
	const auto hp2 = FirFilter<float, TIME_DOMAIN>(numTaps, Highpass(WINDOWED).Cutoff(cutoff).Window(windows::blackman.operator()<float, TIME_DOMAIN>(numTaps)));

	const auto bp1 = FirFilter<float, TIME_DOMAIN>(numTaps, Bandpass(WINDOWED).Band(bandLow, bandHigh).Window(windows::blackman));
	const auto bp2 = FirFilter<float, TIME_DOMAIN>(numTaps, Bandpass(WINDOWED).Band(bandLow, bandHigh).Window(windows::blackman.operator()<float, TIME_DOMAIN>(numTaps)));

	const auto bs1 = FirFilter<float, TIME_DOMAIN>(numTaps, Bandstop(WINDOWED).Band(bandLow, bandHigh).Window(windows::blackman));
	const auto bs2 = FirFilter<float, TIME_DOMAIN>(numTaps, Bandstop(WINDOWED).Band(bandLow, bandHigh).Window(windows::blackman.operator()<float, TIME_DOMAIN>(numTaps)));

	REQUIRE(Max(Abs(lp1 - lp2)) < 1e-4f);
	REQUIRE(Max(Abs(hp1 - hp2)) < 1e-4f);
	REQUIRE(Max(Abs(bp1 - bp2)) < 1e-4f);
	REQUIRE(Max(Abs(bs1 - bs2)) < 1e-4f);
}


//------------------------------------------------------------------------------
// Least squares method
//------------------------------------------------------------------------------

TEST_CASE("Least squares low-pass", "[FIR]") {
	constexpr size_t numTaps = 255;
	static constexpr float cutoffBegin = 0.28f;
	static constexpr float cutoffEnd = 0.32f;

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Lowpass(LEAST_SQUARES).Cutoff(cutoffBegin, cutoffEnd));
	REQUIRE(impulse.Size() == numTaps);
	REQUIRE(IsSymmetric(impulse));
	REQUIRE(Sum(impulse) == Approx(1).margin(0.01));

	RequireResponse(impulse, { { cutoffBegin - 0.01f, 1.0f }, { cutoffEnd + 0.01f, 0.0f } });
}

TEST_CASE("Least squares high-pass", "[FIR]") {
	constexpr size_t numTaps = 255;
	static constexpr float cutoffBegin = 0.28f;
	static constexpr float cutoffEnd = 0.32f;

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Highpass(LEAST_SQUARES).Cutoff(cutoffBegin, cutoffEnd));
	REQUIRE(impulse.Size() == numTaps);
	REQUIRE(IsSymmetric(impulse));
	REQUIRE(Sum(impulse) == Approx(0).margin(0.01));

	RequireResponse(impulse, { { cutoffBegin - 0.04f, 0.0f }, { cutoffEnd + 0.04f, 1.0f } });
}


TEST_CASE("Least squares band-pass", "[FIR]") {
	constexpr size_t numTaps = 255;
	static constexpr float bandLowBegin = 0.28f;
	static constexpr float bandLowEnd = 0.32f;
	static constexpr float bandHighBegin = 0.58f;
	static constexpr float bandHighEnd = 0.65f;

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Bandpass(LEAST_SQUARES).Band(bandLowBegin, bandLowEnd, bandHighBegin, bandHighEnd));
	REQUIRE(impulse.Size() == numTaps);
	REQUIRE(IsSymmetric(impulse));
	REQUIRE(Sum(impulse) == Approx(0).margin(0.01));

	RequireResponse(impulse, {
								 { bandLowBegin - 0.01f, 0.0f },
								 { bandLowEnd + 0.01f, 1.0f },
								 { bandHighBegin - 0.01f, 1.0f },
								 { bandHighEnd + 0.01f, 0.0f },
							 });
}

TEST_CASE("Least squares band-stop", "[FIR]") {
	constexpr size_t numTaps = 255;
	static constexpr float bandLowBegin = 0.28f;
	static constexpr float bandLowEnd = 0.32f;
	static constexpr float bandHighBegin = 0.58f;
	static constexpr float bandHighEnd = 0.65f;

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Bandstop(LEAST_SQUARES).Band(bandLowBegin, bandLowEnd, bandHighBegin, bandHighEnd));
	REQUIRE(impulse.Size() == numTaps);
	REQUIRE(IsSymmetric(impulse));
	REQUIRE(Sum(impulse) == Approx(1.0).margin(0.01));

	RequireResponse(impulse, {
								 { bandLowBegin - 0.01f, 1.0f },
								 { bandLowEnd + 0.01f, 0.0f },
								 { bandHighBegin - 0.01f, 0.0f },
								 { bandHighEnd + 0.01f, 1.0f },
							 });
}

TEST_CASE("Least squares arbitrary", "[FIR]") {
	constexpr size_t numTaps = 255;

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Arbitrary(LEAST_SQUARES).Response(TestArbitraryResponse));
	REQUIRE(impulse.Size() == numTaps);
	REQUIRE(IsSymmetric(impulse));

	RequireResponse(impulse, {
								 { 0.12f, TestArbitraryResponse(0.12f) },
								 { 0.12f, TestArbitraryResponse(0.12f) },
								 { 0.32f, TestArbitraryResponse(0.32f) },
								 { 0.67f, TestArbitraryResponse(0.67f) },
								 { 0.88f, TestArbitraryResponse(0.88f) },
							 });
}

TEST_CASE("Least squares hilbert magnitude", "[FIR]") {
	const float transition = 0.03f;
	const auto odd = FirFilter<float, TIME_DOMAIN>(155, Hilbert(LEAST_SQUARES).TransitionWidth(transition));
	const auto even = FirFilter<float, TIME_DOMAIN>(154, Hilbert(LEAST_SQUARES).TransitionWidth(transition));

	std::vector<std::pair<float, float>> requiredOdd = {
		{ 0.031f, 1.0f },
		{ 0.5f, 1.0f },
		{ 0.969f, 1.0f },
	};
	std::vector<std::pair<float, float>> requiredEven = {
		{ 0.062f, 1.0f },
		{ 0.5f, 1.0f },
		{ 0.999f, 1.0f },
	};
	RequireResponse(odd, requiredOdd, 0.01);
	RequireResponse(even, requiredEven, 0.01);
	REQUIRE(MeasureResponse(0.020f, odd) < 0.95f);
	REQUIRE(MeasureResponse(0.980f, odd) < 0.95f);
	REQUIRE(MeasureResponse(0.020f, even) < 0.95f);
}

TEST_CASE("Least squares weights", "[FIR]") {
	const auto response = [](float f) {
		if (f < 0.5f) {
			return 1.0f;
		}
		return 0.0f;
	};
	const auto weightL = [](float f) {
		if (f < 0.45f) {
			return 3.0f;
		}
		if (f < 0.55f) {
			return 0.0f;
		}
		return 1.0f;
	};
	const auto weightH = [](float f) {
		if (f < 0.45f) {
			return 1.0f;
		}
		if (f < 0.55f) {
			return 0.0f;
		}
		return 3.0f;
	};

	auto filterL = FirFilter<float, TIME_DOMAIN>(27, Arbitrary(LEAST_SQUARES).Response(response).Weight(weightL));
	auto filterH = FirFilter<float, TIME_DOMAIN>(27, Arbitrary(LEAST_SQUARES).Response(response).Weight(weightH));
	filterL.Resize(1024, 0.0f);
	filterH.Resize(1024, 0.0f);
	const auto responseL = Abs(FourierTransform(filterL, false));
	const auto responseH = Abs(FourierTransform(filterH, false));
	const float stdLL = StandardDeviation(AsView(responseL).SubSignal(0, 230));
	const float stdLH = StandardDeviation(AsView(responseL).SubSignal(280));
	const float stdHL = StandardDeviation(AsView(responseH).SubSignal(0, 230));
	const float stdHH = StandardDeviation(AsView(responseH).SubSignal(280));
	REQUIRE(stdLL < stdHL * 0.6f);
	REQUIRE(stdLH * 0.6f > stdHH);
}

//------------------------------------------------------------------------------
// Hilbert band transform special checks
//------------------------------------------------------------------------------

TEST_CASE("Hilbert odd form", "[FIR]") {
	const auto filter = FirFilter<float, TIME_DOMAIN>(247, Hilbert(WINDOWED));
	REQUIRE(filter.Size() == 247);
	REQUIRE(IsAntiSymmetric(filter));
	const auto nonZeroSamples = Decimate(filter, 2);
	const auto zeroSamples = Decimate(AsView(filter).SubSignal(1), 2);
	REQUIRE(Max(zeroSamples) == 0.0f);
	REQUIRE(Min(Abs(nonZeroSamples)) > 0.0f);
	const auto firstHalf = AsView(nonZeroSamples).SubSignal(0, nonZeroSamples.Size() / 2);
	const auto secondHalf = AsView(nonZeroSamples).SubSignal(nonZeroSamples.Size() / 2);
	REQUIRE(Max(firstHalf) < 0.0f);
	REQUIRE(Min(secondHalf) > 0.0f);
}

TEST_CASE("Hilbert even form", "[FIR]") {
	const auto filter = FirFilter<float, TIME_DOMAIN>(246, Hilbert(WINDOWED));
	REQUIRE(filter.Size() == 246);
	REQUIRE(IsAntiSymmetric(filter));
	REQUIRE(Min(Abs(filter)) > 0.0f);
	const auto firstHalf = AsView(filter).SubSignal(0, filter.Size() / 2);
	const auto secondHalf = AsView(filter).SubSignal(filter.Size() / 2);
	REQUIRE(Max(firstHalf) < 0.0f);
	REQUIRE(Min(secondHalf) > 0.0f);
}

TEST_CASE("Hilbert odd small form", "[FIR]") {
	const auto filter = FirFilter<float, TIME_DOMAIN>(19, Hilbert(WINDOWED));
	REQUIRE(filter.Size() == 19);
	REQUIRE(IsAntiSymmetric(filter));
	const auto nonZeroSamples = Decimate(filter, 2);
	const auto zeroSamples = Decimate(AsView(filter).SubSignal(1), 2);
	REQUIRE(Max(zeroSamples) == 0.0f);
	REQUIRE(Min(Abs(nonZeroSamples)) > 0.0f);
	const auto firstHalf = AsView(nonZeroSamples).SubSignal(0, nonZeroSamples.Size() / 2);
	const auto secondHalf = AsView(nonZeroSamples).SubSignal(nonZeroSamples.Size() / 2);
	REQUIRE(Max(firstHalf) < 0.0f);
	REQUIRE(Min(secondHalf) > 0.0f);
}

TEST_CASE("Hilbert even small form", "[FIR]") {
	const auto filter = FirFilter<float, TIME_DOMAIN>(10, Hilbert(WINDOWED));
	REQUIRE(filter.Size() == 10);
	REQUIRE(IsAntiSymmetric(filter));
	REQUIRE(Min(Abs(filter)) > 0.0f);
	const auto firstHalf = AsView(filter).SubSignal(0, filter.Size() / 2);
	const auto secondHalf = AsView(filter).SubSignal(filter.Size() / 2);
	REQUIRE(Max(firstHalf) < 0.0f);
	REQUIRE(Min(secondHalf) > 0.0f);
}


TEST_CASE("Hilbert odd phase shift", "[FIR]") {
	constexpr size_t testSignalSize = 4096;
	const auto filter = FirFilter<float, TIME_DOMAIN>(377, Hilbert(WINDOWED));
	const auto testSignal = SineWave<float, TIME_DOMAIN>(testSignalSize, testSignalSize, 60.0) * GaussianWindow<float, TIME_DOMAIN>(testSignalSize, 0.25);
	const auto imaginarySignal = Convolution(filter, testSignal, convolution::central);
	const auto realSignal = AsConstView(testSignal).SubSignal(filter.Size() / 2, imaginarySignal.Size());
	REQUIRE(std::abs(DotProduct(realSignal, imaginarySignal) / testSignalSize) < 0.000001f);
	REQUIRE(Mean(realSignal) == Approx(Mean(imaginarySignal)).margin(0.001f));
}

TEST_CASE("Hilbert even phase shift", "[FIR]") {
	constexpr size_t testSignalSize = 4096;
	const auto filter = FirFilter<float, TIME_DOMAIN>(376, Hilbert(WINDOWED));
	const auto testSignal = SineWave<float, TIME_DOMAIN>(testSignalSize, testSignalSize, 60.0) * GaussianWindow<float, TIME_DOMAIN>(testSignalSize, 0.25);
	const auto imaginarySignal = Convolution(filter, testSignal, convolution::central);
	const auto realSignal = AsConstView(testSignal).SubSignal(filter.Size() / 2, imaginarySignal.Size());
	REQUIRE(std::abs(DotProduct(realSignal, imaginarySignal) / testSignalSize) < 0.01f);
	REQUIRE(Mean(realSignal) == Approx(Mean(imaginarySignal)).margin(0.001f));
}