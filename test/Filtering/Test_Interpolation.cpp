#include "../TestUtils.hpp"

#include <dspbb/Filtering/Interpolation.hpp>
#include <dspbb/Math/Convolution.hpp>

#include <catch2/catch.hpp>
#include <cmath>

using namespace dspbb;

auto MakeRamp(size_t size) {
	Signal<float> signal;
	for (size_t i = 0; i < size; ++i) {
		signal.PushBack(float(i));
	}
	return signal;
}


template <class SignalT>
SignalT InterpolateRefImpl(const SignalT& signal, const SignalT& filter, size_t rate, size_t offset, size_t length) {
	return Convolution(Expand(signal, rate), filter, offset, length) * rate;
}


TEST_CASE("Decimate", "[Interpolation]") {
	const Signal<float> s = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
	const Signal<float> d = Decimate(s, 3);
	REQUIRE(d.Size() == 4);
	REQUIRE(d[0] == 1);
	REQUIRE(d[1] == 4);
	REQUIRE(d[2] == 7);
	REQUIRE(d[3] == 10);
}


TEST_CASE("Expand", "[Interpolation]") {
	const Signal<float> s = { 1, 2, 3 };
	const Signal<float> e = Expand(s, 3);
	const Signal<float> exp = { 1, 0, 0, 2, 0, 0, 3, 0, 0 };

	REQUIRE(e.Size() == 9);
	REQUIRE(Max(Abs(e - exp)) == Approx(0.0f));
}

TEST_CASE("Interpolation full", "[Interpolation]") {
	constexpr int interpRate = 5;
	constexpr int signalSize = 1024;

	for (const int filterSize : { 31, 33, 2047 }) {
		const auto signal = RandomSignal<float, TIME_DOMAIN>(signalSize);
		const auto filter = FirFilter<float, TIME_DOMAIN>(filterSize, Lowpass(WINDOWED).Cutoff(1.0f / interpRate));
		const auto polyphase = PolyphaseDecompose(filter, interpRate);

		const auto length = ConvolutionLength(signal.Size() * interpRate, filter.Size(), CONV_FULL);
		const auto reference = InterpolateRefImpl(signal, filter, interpRate, 0, length);
		const auto answer = Interpolate(signal, polyphase, 0, length);

		INFO("filterSize=" << filterSize);
		REQUIRE(reference.Size() == answer.Size());
		REQUIRE(Max(Abs(reference - answer)) < 1e-6f);
	}
}

TEST_CASE("Interpolation central", "[Interpolation]") {
	constexpr int interpRate = 5;
	constexpr int signalSize = 1024;

	for (const int filterSize : { 31, 33, 2047 }) {
		const auto signal = RandomSignal<float, TIME_DOMAIN>(signalSize);
		const auto filter = FirFilter<float, TIME_DOMAIN>(filterSize, Lowpass(WINDOWED).Cutoff(1.0f / interpRate));
		const auto polyphase = PolyphaseDecompose(filter, interpRate);

		const auto length = ConvolutionLength(signal.Size() * interpRate, filter.Size(), CONV_CENTRAL);
		const auto reference = InterpolateRefImpl(signal, filter, interpRate, filterSize - 1, length);
		const auto answer = Interpolate(signal, polyphase, filterSize - 1, length);

		INFO("filterSize=" << filterSize);
		REQUIRE(reference.Size() == answer.Size());
		REQUIRE(Max(Abs(reference - answer)) < 1e-6f);
	}
}


TEST_CASE("Resampling length full", "[Interpolation]") {
	SECTION("Upsample exact") {
		constexpr std::pair<uint64_t, uint64_t> sampleRates = { 2, 3 };
		constexpr size_t signalSize = 2000;
		constexpr size_t filterSize = 1001;
		constexpr size_t numPhases = 5;

		constexpr auto size = ResamplingLength(signalSize, filterSize, numPhases, sampleRates, CONV_FULL);
		const double real = double(size.first) / double(size.second);
		REQUIRE(real == Approx(16500.0 / 5).margin(0.01));
	}
	SECTION("Upsample inexact") {
		constexpr std::pair<uint64_t, uint64_t> sampleRates = { 3, 5 };
		constexpr size_t signalSize = 2000;
		constexpr size_t filterSize = 1001;
		constexpr size_t numPhases = 5;

		constexpr auto size = ResamplingLength(signalSize, filterSize, numPhases, sampleRates, CONV_FULL);
		const double real = double(size.first) / double(size.second);
		REQUIRE(real == Approx(18333.333 / 5).margin(0.01));
	}
	SECTION("Downsample exact") {
		constexpr std::pair<uint64_t, uint64_t> sampleRates = { 11000, 3500 };
		constexpr size_t signalSize = 2000;
		constexpr size_t filterSize = 1001;
		volatile size_t numPhases = 5;

		const auto size = ResamplingLength(signalSize, filterSize, numPhases, sampleRates, CONV_FULL);
		const double real = double(size.first) / double(size.second);
		REQUIRE(real == Approx(3500.0 / 5).margin(0.01));
	}
	SECTION("Downsample inexact") {
		constexpr std::pair<uint64_t, uint64_t> sampleRates = { 22000, 7001 };
		constexpr size_t signalSize = 2000;
		constexpr size_t filterSize = 1001;
		constexpr size_t numPhases = 5;

		constexpr auto size = ResamplingLength(signalSize, filterSize, numPhases, sampleRates, CONV_FULL);
		const double real = double(size.first) / double(size.second);
		REQUIRE(real == Approx(3500.5 / 5).margin(0.01));
	}
}

TEST_CASE("Resampling length central", "[Interpolation]") {
	SECTION("Upsample exact") {
		constexpr std::pair<uint64_t, uint64_t> sampleRates = { 9000, 14000 };
		constexpr size_t signalSize = 2000;
		constexpr size_t filterSize = 1001;
		constexpr size_t numPhases = 5;

		constexpr auto size = ResamplingLength(signalSize, filterSize, numPhases, sampleRates, CONV_CENTRAL);
		const double real = double(size.first) / double(size.second);
		REQUIRE(real == Approx(14000.0 / 5).margin(0.01));
	}
	SECTION("Upsample inexact") {
		constexpr std::pair<uint64_t, uint64_t> sampleRates = { 27000, 14000 };
		constexpr size_t signalSize = 2000;
		constexpr size_t filterSize = 1001;
		constexpr size_t numPhases = 5;

		constexpr auto size = ResamplingLength(signalSize, filterSize, numPhases, sampleRates, CONV_CENTRAL);
		const double real = double(size.first) / double(size.second);
		REQUIRE(real == Approx(4666.667 / 5).margin(0.01));
	}
	SECTION("Downsample exact") {
		constexpr std::pair<uint64_t, uint64_t> sampleRates = { 9000, 3500 };
		constexpr size_t signalSize = 2000;
		constexpr size_t filterSize = 1001;
		volatile size_t numPhases = 5;

		const auto size = ResamplingLength(signalSize, filterSize, numPhases, sampleRates, CONV_CENTRAL);
		const double real = double(size.first) / double(size.second);
		REQUIRE(real == Approx(3500.0 / 5).margin(0.01));
	}
	SECTION("Downsample inexact") {
		constexpr std::pair<uint64_t, uint64_t> sampleRates = { 18000, 7001 };
		constexpr size_t signalSize = 2000;
		constexpr size_t filterSize = 1001;
		constexpr size_t numPhases = 5;

		constexpr auto size = ResamplingLength(signalSize, filterSize, numPhases, sampleRates, CONV_CENTRAL);
		const double real = double(size.first) / double(size.second);
		REQUIRE(real == Approx(3500.5 / 5).margin(0.01));
	}
}

TEST_CASE("Resampling change sample rate", "[Interpolation]") {
	constexpr int inputRate = 7;
	constexpr int outputRate = 17;

	constexpr std::pair originalSample = { 28ull, 42ull };

	SECTION("Regular") {
		constexpr auto newSample = resample::ChangeSampleRate({ inputRate, outputRate }, originalSample, false);

		const double outputIndexReal = double(originalSample.first) / double(originalSample.second);
		const double inputIndexRealExpected = outputIndexReal / double(inputRate) * double(outputRate);
		const double inputIndexReal = double(newSample.first) / double(newSample.second);

		REQUIRE(inputIndexReal == Approx(inputIndexRealExpected));
	}
	SECTION("Simplify") {
		constexpr auto newSample = resample::ChangeSampleRate({ inputRate, outputRate }, originalSample, true);

		const double outputIndexReal = double(originalSample.first) / double(originalSample.second);
		const double inputIndexRealExpected = outputIndexReal / double(inputRate) * double(outputRate);
		const double inputIndexReal = double(newSample.first) / double(newSample.second);

		REQUIRE(inputIndexReal == Approx(inputIndexRealExpected));
	}
}

TEST_CASE("Resampling input index 2 samples", "[Interpolation]") {
	SECTION("Zero weight") {
		const auto [firstSample, secondSample] = resample::InputIndex2Sample({ 43, 7 }, 7);
		REQUIRE(firstSample.inputIndex == 6);
		REQUIRE(firstSample.phaseIndex == 1);
		REQUIRE(firstSample.weight == 7);

		REQUIRE(secondSample.inputIndex == 6);
		REQUIRE(secondSample.phaseIndex == 2);
		REQUIRE(secondSample.weight == 0);
	}
	SECTION("Split weight") {
		const auto [firstSample, secondSample] = resample::InputIndex2Sample({ 87, 14 }, 5);
		REQUIRE(firstSample.inputIndex == 6);
		REQUIRE(firstSample.phaseIndex == 1);
		REQUIRE(firstSample.weight == 13);

		REQUIRE(secondSample.inputIndex == 6);
		REQUIRE(secondSample.phaseIndex == 2);
		REQUIRE(secondSample.weight == 1);
	}
	SECTION("Rollover") {
		const auto [firstSample, secondSample] = resample::InputIndex2Sample({ 27, 14 }, 5);
		REQUIRE(firstSample.inputIndex == 1);
		REQUIRE(firstSample.phaseIndex == 4);
		REQUIRE(firstSample.weight == 5);

		REQUIRE(secondSample.inputIndex == 2);
		REQUIRE(secondSample.phaseIndex == 0);
		REQUIRE(secondSample.weight == 9);
	}
}

TEST_CASE("Resampling dot product sample", "[Interpolation]") {
	const Signal<int> signal = { 1, 2, 3, 6, 5, 7 };
	const Signal<int> filter = { -1, 3, -2 }; // Convolution: -2, 3, -1
	REQUIRE(-2 == resample::DotProductSample(signal, filter, 0));
	REQUIRE(-1 == resample::DotProductSample(signal, filter, 2));
	REQUIRE(-5 == resample::DotProductSample(signal, filter, 5));
	REQUIRE(-7 == resample::DotProductSample(signal, filter, 7));
}

TEST_CASE("Resampling filter cutoff", "[Interpolation]") {
	REQUIRE(ResamplingFilterCutoff({ 4, 6 }, 5) == Approx(0.2));
	REQUIRE(ResamplingFilterCutoff({ 6, 4 }, 5) == Approx(0.1333333333));
	REQUIRE(ResamplingFilterCutoff({ 4, 71 }, 12) == Approx(0.0833333333));
	REQUIRE(ResamplingFilterCutoff({ 40, 6 }, 12) == Approx(0.0125));
}


template <class SignalT, class SignalU>
auto ResampledSimilarity(std::pair<uint64_t, uint64_t> sampleRates, SignalT original, SignalU resampled) {
	const size_t rescale = std::max(original.Size() / sampleRates.first, resampled.Size() / sampleRates.second) + 1;
	original.Resize(rescale * sampleRates.first);
	resampled.Resize(rescale * sampleRates.second);

	const auto fftSignal = Abs(Fft(original, FFT_HALF));
	const auto fftResampled = Abs(Fft(resampled, FFT_HALF));

	const size_t fftCompareSize = std::min(fftSignal.Size(), fftResampled.Size());
	const auto fftSignalCompare = AsView(fftSignal).SubSignal(0, fftCompareSize);
	const auto fftResampledCompare = AsView(fftResampled).SubSignal(0, fftCompareSize);

	const auto similarity = DotProduct(fftSignalCompare, fftResampledCompare) / Norm(fftSignalCompare) / Norm(fftResampledCompare);

	return similarity;
}


TEST_CASE("Resampling spectrum invariance - upsample mild", "[Interpolation]") {
	constexpr int inputRate = 7;
	constexpr int outputRate = 11;
	constexpr int supersamplingRate = 16;
	constexpr int signalSize = 1024;
	constexpr auto filterCutoff = ResamplingFilterCutoff({ inputRate, outputRate }, supersamplingRate);

	for (const int filterSize : { 513, 2047 }) {
		const auto signal = RandomSignal<float, TIME_DOMAIN>(signalSize);
		const auto filter = FirFilter<float, TIME_DOMAIN>(filterSize, Lowpass(WINDOWED).Cutoff(filterCutoff));
		const auto polyphase = PolyphaseDecompose(filter, supersamplingRate);

		const auto length = ResamplingLength(signalSize, filterSize, supersamplingRate, { inputRate, outputRate }, CONV_FULL);
		const auto resampled = Resample(signal, polyphase, { inputRate, outputRate }, { 0, 1 }, length.first / length.second);
		const auto similarity = ResampledSimilarity({ inputRate, outputRate }, signal, resampled);

		INFO("filterSize=" << filterSize)
		REQUIRE(similarity > 0.98f);
	}
}

TEST_CASE("Resampling spectrum invariance - upsample strong", "[Interpolation]") {
	constexpr int inputRate = 9;
	constexpr int outputRate = 210;
	constexpr int supersamplingRate = 32;
	constexpr int signalSize = 2048;
	constexpr auto filterCutoff = ResamplingFilterCutoff({ inputRate, outputRate }, supersamplingRate);

	for (const int filterSize : { 1023, 4047 }) {
		const auto signal = RandomSignal<float, TIME_DOMAIN>(signalSize);
		const auto filter = FirFilter<float, TIME_DOMAIN>(filterSize, Lowpass(WINDOWED).Cutoff(filterCutoff));
		const auto polyphase = PolyphaseDecompose(filter, supersamplingRate);

		const auto length = ResamplingLength(signalSize, filterSize, supersamplingRate, { inputRate, outputRate }, CONV_FULL);
		const auto resampled = Resample(signal, polyphase, { inputRate, outputRate }, { 0, 1 }, length.first / length.second);
		const auto similarity = ResampledSimilarity({ inputRate, outputRate }, signal, resampled);

		INFO("filterSize=" << filterSize)
		REQUIRE(similarity > 0.98f);
	}
}

TEST_CASE("Resampling spectrum invariance - downsample mild", "[Interpolation]") {
	constexpr int inputRate = 11;
	constexpr int outputRate = 7;
	constexpr int supersamplingRate = 16;
	constexpr int signalSize = 16384;
	constexpr auto filterCutoff = ResamplingFilterCutoff({ inputRate, outputRate }, supersamplingRate);

	for (const int filterSize : { 4095, 20001 }) {
		const auto signal = RandomSignal<float, TIME_DOMAIN>(signalSize);
		const auto filter = FirFilter<float, TIME_DOMAIN>(filterSize, Lowpass(WINDOWED).Cutoff(filterCutoff));
		const auto polyphase = PolyphaseDecompose(filter, supersamplingRate);

		const auto length = ResamplingLength(signalSize, filterSize, supersamplingRate, { inputRate, outputRate }, CONV_FULL);
		const auto resampled = Resample(signal, polyphase, { inputRate, outputRate }, { 0, 1 }, length.first / length.second);
		const auto similarity = ResampledSimilarity({ inputRate, outputRate }, signal, resampled);

		INFO("filterSize=" << filterSize)
		REQUIRE(similarity > 0.98f);
	}
}

TEST_CASE("Resampling spectrum invariance - downsample strong", "[Interpolation]") {
	constexpr int inputRate = 210;
	constexpr int outputRate = 9;
	constexpr int supersamplingRate = 16;
	constexpr int signalSize = 16384;
	constexpr auto filterCutoff = ResamplingFilterCutoff({ inputRate, outputRate }, supersamplingRate);

	for (const int filterSize : { 4095, 20001 }) {
		const auto signal = RandomSignal<float, TIME_DOMAIN>(signalSize);
		const auto filter = FirFilter<float, TIME_DOMAIN>(filterSize, Lowpass(WINDOWED).Cutoff(filterCutoff));
		const auto polyphase = PolyphaseDecompose(filter, supersamplingRate);

		const auto length = ResamplingLength(signalSize, filterSize, supersamplingRate, { inputRate, outputRate }, CONV_FULL);
		const auto resampled = Resample(signal, polyphase, { inputRate, outputRate }, { 0, 1 }, length.first / length.second);
		const auto similarity = ResampledSimilarity({ inputRate, outputRate }, signal, resampled);

		INFO("filterSize=" << filterSize)
		REQUIRE(similarity > 0.98f);
	}
}


template <class SignalT>
double FindCrossing(const SignalT& signal, double value) {
	const auto it = std::adjacent_find(signal.begin(), signal.end(), [&value](auto left, auto right) {
		return left <= value && value < right;
	});
	if (it != signal.end()) {
		const size_t firstIndex = it - signal.begin();
		const auto difference = it[1] - it[0];
		const auto t = (value - it[0]) / difference;
		return double(firstIndex) + double(t);
	}
	return -1.0;
}


TEST_CASE("Resampling delay - upsample mild", "[Interpolation]") {
	// Resample a ramp function.
	// The exact crossing (i.e. f(x) = 10, x = ?) can be easily found by linear interpolation.
	// The exact crossing can be used to correlate delays on the input and output signals.

	constexpr int inputRate = 7;
	constexpr int outputRate = 11;
	constexpr int supersamplingRate = 16;
	constexpr int signalSize = 1024;
	constexpr auto filterCutoff = ResamplingFilterCutoff({ inputRate, outputRate }, supersamplingRate);

	for (const int filterSize : { 513, 2047 }) {
		auto signal = Signal<float>(signalSize);
		std::iota(signal.begin(), signal.end(), 0.0f);
		const auto filter = FirFilter<float, TIME_DOMAIN>(filterSize, Lowpass(WINDOWED).Cutoff(filterCutoff));
		const auto polyphase = PolyphaseNormalized(PolyphaseDecompose(filter, supersamplingRate));

		const auto length = ResamplingLength(signalSize, filterSize, supersamplingRate, { inputRate, outputRate }, CONV_FULL);
		const auto resampled = Resample(signal, polyphase, { inputRate, outputRate }, { 0, 1 }, length.first / length.second);

		const double crossingSignal = FindCrossing(signal, 500.0);
		const double crossingResampled = FindCrossing(resampled, 500.0);
		const auto resamplingDelay = ResamplingDelay(filterSize, supersamplingRate, { inputRate, outputRate });
		const double crossingExpected = double(resamplingDelay.first) / resamplingDelay.second + crossingSignal * outputRate / inputRate;

		INFO("filterSize=" << filterSize)
		REQUIRE(crossingExpected == Approx(crossingResampled));
	}
}


TEST_CASE("Interplation continuation calculation", "[Interpolation]") {
	// Interpolate a ramp and see if it's continuous.
	REQUIRE(false);
}


TEST_CASE("Resampling continuation calculation", "[Interpolation]") {
	constexpr size_t numPhases = 6;
	constexpr size_t filterSize = 31;
	constexpr std::pair sampleRates = { 4ULL, 7ULL };

	SECTION("Initial point") {
		constexpr std::pair nextOutputSample = { 0ULL, 1ULL };
		const auto [inputIndex, startPoint] = resample::Continuation(nextOutputSample, filterSize, numPhases, sampleRates);

		REQUIRE(inputIndex == 0);
		REQUIRE(double(startPoint.first) / double(startPoint.second) == Approx(0));
	}
	SECTION("One off") {
		constexpr std::pair nextOutputSample = { 7ULL, 7ULL };
		const auto [inputIndex, startPoint] = resample::Continuation(nextOutputSample, filterSize, numPhases, sampleRates);

		REQUIRE(inputIndex == 0);
		REQUIRE(double(startPoint.first) / double(startPoint.second) == Approx(1));
	}
	SECTION("Middle point") {
		constexpr std::pair nextOutputSample = { 6ULL * 7ULL, 4ULL };
		const auto [inputIndex, startPoint] = resample::Continuation(nextOutputSample, filterSize, numPhases, sampleRates);

		REQUIRE(inputIndex == 1);
		const double expectedTotalOffset = double(nextOutputSample.first) / double(nextOutputSample.second);
		const double actualTotalOffset = double(inputIndex) * sampleRates.second / sampleRates.first + double(startPoint.first) / double(startPoint.second);
		REQUIRE(expectedTotalOffset == Approx(actualTotalOffset));
	}
}