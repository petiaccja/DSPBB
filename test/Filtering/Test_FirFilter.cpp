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
void RequireResponse(SignalT&& impulse, std::vector<std::pair<float, float>> desired, double margin = 0.03) {
	for (auto& point : desired) {
		const auto response = MeasureResponse(point.first, impulse);
		REQUIRE(response == Approx(point.second).margin(margin));
	}
}

constexpr auto TestArbitraryResponse = [](float x) {
	return 2.0f * x - 1.5f * x * x - 0.5f * std::pow(x - 1.f, 3.f);
};


//------------------------------------------------------------------------------
// Window method
//------------------------------------------------------------------------------

TEST_CASE("Windowed low-pass", "[FirFilter]") {
	constexpr size_t numTaps = 255;
	static constexpr float cutoff = 0.3f;

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Lowpass(WINDOWED).Cutoff(cutoff).Window(windows::blackman));
	REQUIRE(impulse.Size() == numTaps);
	REQUIRE(IsSymmetric(impulse));
	REQUIRE(Sum(impulse) == Approx(1));

	RequireResponse(impulse, { { cutoff - 0.04f, 1.0f }, { cutoff + 0.04f, 0.0f } });
}

TEST_CASE("Windowed high-pass", "[FirFilter]") {
	constexpr size_t numTaps = 255;
	static constexpr float cutoff = 0.3f;

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Highpass(WINDOWED).Cutoff(cutoff).Window(windows::blackman));
	REQUIRE(impulse.Size() == numTaps);
	REQUIRE(IsSymmetric(impulse));
	REQUIRE(Sum(impulse) == Approx(0).margin(0.01));

	RequireResponse(impulse, { { cutoff - 0.04f, 0.0f }, { cutoff + 0.04f, 1.0f } });
}

TEST_CASE("Windowed band-pass", "[FirFilter]") {
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

TEST_CASE("Windowed band-stop", "[FirFilter]") {
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

TEST_CASE("Windowed arbitrary", "[FirFilter]") {
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

TEST_CASE("Windowed methods equal", "[FirFilter]") {
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


/*
TEST_CASE("Windowed arbitrary filter", "[FirFilter]") {
	constexpr size_t numTaps = 255;
	constexpr std::array<float, 4> amplitudes = { 1.0f, 0.2f, 0.6f, 1.2f };
	constexpr std::array<float, 4> frequencies = { 0.0625f, 0.1875f, 0.375f, 0.75f };

	const auto impulse1 = FirFilter<float, TIME_DOMAIN>(numTaps, Arbitrary(WINDOWED).Response(response).Window(windows::hamming));
	const auto impulse2 = FirFilter<float, TIME_DOMAIN>(numTaps, Arbitrary(WINDOWED).Response(response).Window(windows::hamming.operator()<float, TIME_DOMAIN>(numTaps)));
	REQUIRE(IsSymmetric(impulse1));
	REQUIRE(impulse1.Size() == numTaps);
	REQUIRE(impulse2.Size() == numTaps);
	REQUIRE(Max(Abs(impulse1 - impulse2)) < 1e-4f);

	for (size_t i = 0; i < amplitudes.size(); ++i) {
		const auto response = MeasureResponse(sampleRate, frequencies[i] * sampleRate / 2.0f, impulse1);
		REQUIRE(response == Approx(amplitudes[i]).margin(0.05f));
	}
}

TEST_CASE("Highpass", "[FirFilter]") {
	constexpr size_t numTaps = 255;
	static constexpr float cutoff = 3800.f;
	const auto normalizedCutoff = NormalizedFrequency(cutoff, sampleRate);

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Highpass(normalizedCutoff), Windowed(windows::hamming));
	REQUIRE(IsSymmetric(impulse));
	REQUIRE(Sum(impulse) < 1e-4f);
	REQUIRE(impulse.Size() == numTaps);

	const float passResponse = MeasureResponse(sampleRate, cutoff * 1.15f, impulse);
	const float stopResponse = MeasureResponse(sampleRate, cutoff * 0.85f, impulse);

	REQUIRE(passResponse > 0.95f);
	REQUIRE(passResponse < 1.05f);
	REQUIRE(stopResponse < 0.05f);
}

TEST_CASE("Bandpass", "[FirFilter]") {
	constexpr size_t numTaps = 255;
	static constexpr float bandLow = 3800.f;
	static constexpr float bandHigh = 14500.f;
	const auto normalizedLow = NormalizedFrequency(bandLow, sampleRate);
	const auto normalizedHigh = NormalizedFrequency(bandHigh, sampleRate);

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Bandpass(normalizedLow, normalizedHigh), Windowed(windows::hamming));
	REQUIRE(IsSymmetric(impulse));
	REQUIRE(Sum(impulse) < 1e-3f);
	REQUIRE(impulse.Size() == numTaps);

	const float lowStopResponse = MeasureResponse(sampleRate, bandLow * 0.9f, impulse);
	const float lowPassResponse = MeasureResponse(sampleRate, bandLow * 1.1f, impulse);
	const float highPassResponse = MeasureResponse(sampleRate, bandHigh * 0.9f, impulse);
	const float highStopResponse = MeasureResponse(sampleRate, bandHigh * 1.1f, impulse);

	REQUIRE(highPassResponse > 0.95f);
	REQUIRE(highPassResponse < 1.05f);
	REQUIRE(highStopResponse < 0.05f);
	REQUIRE(lowPassResponse > 0.95f);
	REQUIRE(lowPassResponse < 1.05f);
	REQUIRE(lowStopResponse < 0.05f);
}


TEST_CASE("Bandstop", "[FirFilter]") {
	constexpr size_t numTaps = 255;
	static constexpr float bandLow = 3800.f;
	static constexpr float bandHigh = 14500.f;
	const auto normalizedLow = NormalizedFrequency(bandLow, sampleRate);
	const auto normalizedHigh = NormalizedFrequency(bandHigh, sampleRate);

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Bandstop(normalizedLow, normalizedHigh), Windowed(windows::hamming));
	REQUIRE(IsSymmetric(impulse));
	REQUIRE(Sum(impulse) == Approx(1).epsilon(0.005f));
	REQUIRE(impulse.Size() == numTaps);

	const float lowPassResponse = MeasureResponse(sampleRate, bandLow * 0.9f, impulse);
	const float lowStopResponse = MeasureResponse(sampleRate, bandLow * 1.1f, impulse);
	const float highStopResponse = MeasureResponse(sampleRate, bandHigh * 0.9f, impulse);
	const float highPassResponse = MeasureResponse(sampleRate, bandHigh * 1.1f, impulse);

	REQUIRE(highPassResponse > 0.95f);
	REQUIRE(highPassResponse < 1.05f);
	REQUIRE(highStopResponse < 0.05f);
	REQUIRE(lowPassResponse > 0.95f);
	REQUIRE(lowPassResponse < 1.05f);
	REQUIRE(lowStopResponse < 0.05f);
}

TEST_CASE("Hilbert odd form", "[Hilbert]") {
	const auto filter = FirFilter<float, TIME_DOMAIN>(247, Hilbert(), Windowed(windows::hamming));
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

TEST_CASE("Hilbert even form", "[Hilbert]") {
	const auto filter = FirFilter<float, TIME_DOMAIN>(246, Hilbert(), Windowed(windows::hamming));
	REQUIRE(filter.Size() == 246);
	REQUIRE(IsAntiSymmetric(filter));
	REQUIRE(Min(Abs(filter)) > 0.0f);
	const auto firstHalf = AsView(filter).SubSignal(0, filter.Size() / 2);
	const auto secondHalf = AsView(filter).SubSignal(filter.Size() / 2);
	REQUIRE(Max(firstHalf) < 0.0f);
	REQUIRE(Min(secondHalf) > 0.0f);
}

TEST_CASE("Hilbert odd small form", "[Hilbert]") {
	const auto filter = FirFilter<float, TIME_DOMAIN>(19, Hilbert(), Windowed(windows::hamming));
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

TEST_CASE("Hilbert even small form", "[Hilbert]") {
	const auto filter = FirFilter<float, TIME_DOMAIN>(10, Hilbert(), Windowed(windows::hamming));
	REQUIRE(filter.Size() == 10);
	REQUIRE(IsAntiSymmetric(filter));
	REQUIRE(Min(Abs(filter)) > 0.0f);
	const auto firstHalf = AsView(filter).SubSignal(0, filter.Size() / 2);
	const auto secondHalf = AsView(filter).SubSignal(filter.Size() / 2);
	REQUIRE(Max(firstHalf) < 0.0f);
	REQUIRE(Min(secondHalf) > 0.0f);
}


TEST_CASE("Hilbert odd response", "[Hilbert]") {
	constexpr size_t testSignalSize = 4096;
	const auto filter = FirFilter<float, TIME_DOMAIN>(377, Hilbert(), Windowed(windows::hamming));
	const auto testSignal = SineWave<float, TIME_DOMAIN>(testSignalSize, testSignalSize, 60.0) * GaussianWindow<float, TIME_DOMAIN>(testSignalSize, 0.25);
	const auto imaginarySignal = Convolution(filter, testSignal, convolution::central);
	const auto realSignal = AsConstView(testSignal).SubSignal(filter.Size() / 2, imaginarySignal.Size());
	REQUIRE(std::abs(DotProduct(realSignal, imaginarySignal) / testSignalSize) < 0.000001f);
	REQUIRE(Mean(realSignal) == Approx(Mean(imaginarySignal)).margin(0.001f));
}

TEST_CASE("Hilbert even response", "[Hilbert]") {
	constexpr size_t testSignalSize = 4096;
	const auto filter = FirFilter<float, TIME_DOMAIN>(376, Hilbert(), Windowed(windows::hamming));
	const auto testSignal = SineWave<float, TIME_DOMAIN>(testSignalSize, testSignalSize, 60.0) * GaussianWindow<float, TIME_DOMAIN>(testSignalSize, 0.25);
	const auto imaginarySignal = Convolution(filter, testSignal, convolution::central);
	const auto realSignal = AsConstView(testSignal).SubSignal(filter.Size() / 2, imaginarySignal.Size());
	REQUIRE(std::abs(DotProduct(realSignal, imaginarySignal) / testSignalSize) < 0.01f);
	REQUIRE(Mean(realSignal) == Approx(Mean(imaginarySignal)).margin(0.001f));
}
*/