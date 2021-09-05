#include <dspbb/Filtering/Interpolation.hpp>
#include <dspbb/Generators/Waveforms.hpp>
#include <dspbb/Math/Statistics.hpp>

#include <array>
#include <catch2/catch.hpp>
#include <dspbb/Filtering/Convolution.hpp>
#include <dspbb/Filtering/FIR.hpp>

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
	while (beg != end.base()) {
		if (*beg != *end) {
			return false;
		}
	}
	return true;
}


//------------------------------------------------------------------------------
// Tests
//------------------------------------------------------------------------------


static constexpr size_t sampleRate = 44100;


TEST_CASE("Windowed Lowpass", "[FirFilter]") {
	constexpr size_t numTaps = 255;
	static constexpr float cutoff = 3800.f;
	const auto normalizedCutoff = NormalizedFrequency(cutoff, sampleRate);


	const auto impulse1 = FirFilter<float, TIME_DOMAIN>(numTaps, Lowpass(normalizedCutoff), Windowed(windows::hamming));
	const auto impulse2 = FirFilter<float, TIME_DOMAIN>(numTaps, Lowpass(normalizedCutoff), Windowed(windows::hamming.operator()<float, TIME_DOMAIN>(numTaps)));
	REQUIRE(Sum(impulse1) == Approx(1));
	REQUIRE(impulse1.Size() == numTaps);
	REQUIRE(impulse2.Size() == numTaps);
	REQUIRE(Max(Abs(impulse1 - impulse2)) < 1e-4f);

	// Generate two signals just above and just below the cutoff and see their attenuation.
	const auto passSignal = GenTestSignal(sampleRate, cutoff * 0.85f);
	const auto rejectSignal = GenTestSignal(sampleRate, cutoff * 1.15f);

	const auto filteredPassSignal = Convolution(passSignal, impulse1, convolution::full);
	const auto filteredRejectSignal = Convolution(rejectSignal, impulse1, convolution::full);

	const float energyPass = SumSquare(passSignal);
	const float energyReject = SumSquare(rejectSignal);
	const float energyFilteredPass = SumSquare(filteredPassSignal);
	const float energyFilteredReject = SumSquare(filteredRejectSignal);

	REQUIRE(energyFilteredPass / energyPass > 0.95f);
	REQUIRE(energyFilteredPass / energyPass < 1.05f);
	REQUIRE(energyFilteredReject / energyReject < 0.05f);
}


TEST_CASE("Windowed arbitrary filter", "[FirFilter]") {
	constexpr size_t numTaps = 255;
	constexpr std::array<float, 4> amplitudes = { 1.0f, 0.2f, 0.6f, 1.2f };
	constexpr std::array<float, 4> frequencies = { 0.0625f, 0.1875f, 0.375f, 0.75f };

	const auto response = [amplitudes](float freq) {
		if (freq < 0.125f) {
			return amplitudes[0];
		}
		if (freq < 0.25f) {
			return amplitudes[1];
		}
		if (freq < 0.5f) {
			return amplitudes[2];
		}
		return amplitudes[3];
	};

	const auto impulse1 = FirFilter<float, TIME_DOMAIN>(numTaps, Arbitrary(response), Windowed(windows::hamming));
	const auto impulse2 = FirFilter<float, TIME_DOMAIN>(numTaps, Arbitrary(response), Windowed(windows::hamming.operator()<float, TIME_DOMAIN>(numTaps)));
	REQUIRE(impulse1.Size() == numTaps);
	REQUIRE(impulse2.Size() == numTaps);
	REQUIRE(Max(Abs(impulse1 - impulse2)) < 1e-4f);

	for (size_t i = 0; i < amplitudes.size(); ++i) {
		const auto signal = GenTestSignal(sampleRate, frequencies[i] * sampleRate / 2.0f);
		const auto filtered = Convolution(signal, impulse1, convolution::full);
		const float energy = std::sqrt(SumSquare(signal));
		const float filteredEnergy = std::sqrt(SumSquare(filtered));
		REQUIRE(filteredEnergy / energy == Approx(amplitudes[i]).margin(0.05f));
	}
}


TEST_CASE("Highpass", "[FirFilter]") {
	constexpr size_t numTaps = 255;
	static constexpr float cutoff = 3800.f;
	const auto normalizedCutoff = NormalizedFrequency(cutoff, sampleRate);

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Highpass(normalizedCutoff), Windowed(windows::hamming));
	REQUIRE(Sum(impulse) < 1e-4f);
	REQUIRE(impulse.Size() == numTaps);

	// Generate two signals just above and just below the cutoff and see their attenuation.
	const auto passSignal = GenTestSignal(sampleRate, cutoff * 1.15f);
	const auto rejectSignal = GenTestSignal(sampleRate, cutoff * 0.85f);

	const auto filteredPassSignal = Convolution(passSignal, impulse, convolution::full);
	const auto filteredRejectSignal = Convolution(rejectSignal, impulse, convolution::full);

	const float energyPass = SumSquare(passSignal);
	const float energyReject = SumSquare(rejectSignal);
	const float energyFilteredPass = SumSquare(filteredPassSignal);
	const float energyFilteredReject = SumSquare(filteredRejectSignal);

	REQUIRE(energyFilteredPass / energyPass > 0.95f);
	REQUIRE(energyFilteredPass / energyPass < 1.05f);
	REQUIRE(energyFilteredReject / energyReject < 0.05f);
}


TEST_CASE("Bandpass", "[FirFilter]") {
	constexpr size_t numTaps = 255;
	static constexpr float bandLow = 3800.f;
	static constexpr float bandHigh = 14500.f;
	const auto normalizedLow = NormalizedFrequency(bandLow, sampleRate);
	const auto normalizedHigh = NormalizedFrequency(bandHigh, sampleRate);

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Bandpass(normalizedLow, normalizedHigh), Windowed(windows::hamming));
	REQUIRE(Sum(impulse) < 1e-3f);
	REQUIRE(impulse.Size() == numTaps);

	auto extended = impulse;
	extended.Resize(44100, 0.0f);
	const auto spectrum = Abs(FourierTransform(extended));

	const auto passSignal1 = GenTestSignal(sampleRate, bandLow * 1.1f);
	const auto passSignal2 = GenTestSignal(sampleRate, bandHigh * 0.9f);
	const auto rejectSignal1 = GenTestSignal(sampleRate, bandLow * 0.9f);
	const auto rejectSignal2 = GenTestSignal(sampleRate, bandHigh * 1.1f);

	const auto filteredPassSignal1 = Convolution(passSignal1, impulse, convolution::full);
	const auto filteredPassSignal2 = Convolution(passSignal2, impulse, convolution::full);
	const auto filteredRejectSignal1 = Convolution(rejectSignal1, impulse, convolution::full);
	const auto filteredRejectSignal2 = Convolution(rejectSignal2, impulse, convolution::full);

	const float energyPass1 = SumSquare(passSignal1);
	const float energyPass2 = SumSquare(passSignal2);
	const float energyReject1 = SumSquare(rejectSignal1);
	const float energyReject2 = SumSquare(rejectSignal2);
	const float energyFilteredPass1 = SumSquare(filteredPassSignal1);
	const float energyFilteredPass2 = SumSquare(filteredPassSignal2);
	const float energyFilteredReject1 = SumSquare(filteredRejectSignal1);
	const float energyFilteredReject2 = SumSquare(filteredRejectSignal2);

	REQUIRE(energyFilteredPass1 / energyPass1 > 0.95f);
	REQUIRE(energyFilteredPass1 / energyPass1 < 1.05f);
	REQUIRE(energyFilteredReject1 / energyReject1 < 0.05f);
	REQUIRE(energyFilteredPass2 / energyPass2 > 0.95f);
	REQUIRE(energyFilteredPass2 / energyPass2 < 1.05f);
	REQUIRE(energyFilteredReject2 / energyReject2 < 0.05f);
}


TEST_CASE("Bandstop", "[FirFilter]") {
	constexpr size_t numTaps = 255;
	static constexpr float bandLow = 3800.f;
	static constexpr float bandHigh = 14500.f;
	const auto normalizedLow = NormalizedFrequency(bandLow, sampleRate);
	const auto normalizedHigh = NormalizedFrequency(bandHigh, sampleRate);

	const auto impulse = FirFilter<float, TIME_DOMAIN>(numTaps, Bandstop(normalizedLow, normalizedHigh), Windowed(windows::hamming));
	REQUIRE(Sum(impulse) == Approx(1).epsilon(0.005f));
	REQUIRE(impulse.Size() == numTaps);

	auto extended = impulse;
	extended.Resize(44100, 0.0f);
	const auto spectrum = Abs(FourierTransform(extended));

	const auto rejectSignal1 = GenTestSignal(sampleRate, bandLow * 1.1f);
	const auto rejectSignal2 = GenTestSignal(sampleRate, bandHigh * 0.9f);
	const auto passSignal1 = GenTestSignal(sampleRate, bandLow * 0.9f);
	const auto passSignal2 = GenTestSignal(sampleRate, bandHigh * 1.1f);

	const auto filteredPassSignal1 = Convolution(passSignal1, impulse, convolution::full);
	const auto filteredPassSignal2 = Convolution(passSignal2, impulse, convolution::full);
	const auto filteredRejectSignal1 = Convolution(rejectSignal1, impulse, convolution::full);
	const auto filteredRejectSignal2 = Convolution(rejectSignal2, impulse, convolution::full);

	const float energyPass1 = SumSquare(passSignal1);
	const float energyPass2 = SumSquare(passSignal2);
	const float energyReject1 = SumSquare(rejectSignal1);
	const float energyReject2 = SumSquare(rejectSignal2);
	const float energyFilteredPass1 = SumSquare(filteredPassSignal1);
	const float energyFilteredPass2 = SumSquare(filteredPassSignal2);
	const float energyFilteredReject1 = SumSquare(filteredRejectSignal1);
	const float energyFilteredReject2 = SumSquare(filteredRejectSignal2);

	REQUIRE(energyFilteredPass1 / energyPass1 > 0.95f);
	REQUIRE(energyFilteredPass1 / energyPass1 < 1.05f);
	REQUIRE(energyFilteredReject1 / energyReject1 < 0.05f);
	REQUIRE(energyFilteredPass2 / energyPass2 > 0.95f);
	REQUIRE(energyFilteredPass2 / energyPass2 < 1.05f);
	REQUIRE(energyFilteredReject2 / energyReject2 < 0.05f);
}



TEST_CASE("Hilbert odd form", "[Hilbert]") {
	const auto filter = FirFilter<float, TIME_DOMAIN>(247, Hilbert(), Windowed(windows::hamming));
	REQUIRE(filter.Size() == 247);
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
	REQUIRE(Min(Abs(filter)) > 0.0f);
	const auto firstHalf = AsView(filter).SubSignal(0, filter.Size() / 2);
	const auto secondHalf = AsView(filter).SubSignal(filter.Size() / 2);
	REQUIRE(Max(firstHalf) < 0.0f);
	REQUIRE(Min(secondHalf) > 0.0f);
}

TEST_CASE("Hilbert odd small form", "[Hilbert]") {
	const auto filter = FirFilter<float, TIME_DOMAIN>(19, Hilbert(), Windowed(windows::hamming));
	REQUIRE(filter.Size() == 19);
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
