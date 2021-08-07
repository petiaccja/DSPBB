#include "dspbb/Math/Statistics.hpp"

#include <dspbb/Filtering/Convolution.hpp>
#include <dspbb/Filtering/FIR.hpp>

#include <catch2/catch.hpp>
#include <array>

using namespace dspbb;




//------------------------------------------------------------------------------
// Helpers
//------------------------------------------------------------------------------

TimeSignal<float> GenTestSignal(size_t sampleRate, float frequency, float length = 1.0f) {
	const size_t numSamples = std::max(size_t(1), size_t(sampleRate * double(length)));
	const float phasePerSample = 2 * pi_v<float> * frequency / sampleRate;
	TimeSignal<float> signal(numSamples);

	float phase = 0.0f;
	for (auto& sample : signal) {
		sample = std::sin(phase);
		phase += phasePerSample;
	}
	return signal;
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


TEST_CASE("Lowpass windowed odd", "[FirFilter]") {
	constexpr size_t numTaps = 255;
	static constexpr float cutoff = 3800.f;
	
	const auto impulse1 = FirLowpassWin<float, TIME_DOMAIN>(NormalizedFrequency(cutoff, sampleRate), numTaps, windows::hamming);
	const auto impulse2 = FirLowpassWin<float, TIME_DOMAIN>(NormalizedFrequency(cutoff, sampleRate), windows::hamming.operator()<float, TIME_DOMAIN>(numTaps));
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


TEST_CASE("Arbitrary filter", "[FirFilter]") {
	constexpr size_t numTaps = 255;
	constexpr std::array<float, 4> amplitudes = { 1.0f, 0.2f, 0.6f, 1.2f };
	constexpr float nyquistLimit = float(sampleRate / 2);
	constexpr std::array<float, 4> frequencies = { 0.0625f, 0.1875f, 0.375f, 0.75f };
	
	Spectrum<float> frequencyResponse(8192, amplitudes[0]);
	std::fill(frequencyResponse.begin() + 1024, frequencyResponse.begin() + 2048, amplitudes[1]);
	std::fill(frequencyResponse.begin() + 2048, frequencyResponse.begin() + 4096, amplitudes[2]);
	std::fill(frequencyResponse.begin() + 4096, frequencyResponse.end(), amplitudes[3]);

	const auto impulse1 = FirArbitraryWin<float, TIME_DOMAIN>(frequencyResponse, numTaps, windows::hamming);
	const auto impulse2 = FirArbitraryWin<float, TIME_DOMAIN>(frequencyResponse, windows::hamming.operator()<float, TIME_DOMAIN>(numTaps));
	REQUIRE(impulse1.Size() == numTaps);
	REQUIRE(impulse2.Size() == numTaps);
	REQUIRE(Max(Abs(impulse1 - impulse2)) < 1e-4f);

	const auto quality = FirQuality(impulse1, frequencyResponse);
	REQUIRE(quality.cosineSimilarity > 0.95f);
	
	for (size_t i = 0; i < amplitudes.size(); ++i) {
		const auto signal = GenTestSignal(sampleRate, frequencies[i]*nyquistLimit);
		const auto filtered = Convolution(signal, impulse1, convolution::full);
		const float energy = std::sqrt(SumSquare(signal));
		const float filteredEnergy = std::sqrt(SumSquare(filtered));
		REQUIRE(std::abs(amplitudes[i]*energy - filteredEnergy)/energy < 0.05f);
	}
}