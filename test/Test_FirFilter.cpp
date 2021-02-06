#include "dspbb/Math/Statistics.hpp"

#include <dspbb/Filtering/Convolution.hpp>
#include <dspbb/Filtering/FIR.hpp>

#include <Catch2/catch.hpp>
#include <array>

using namespace dspbb;


static constexpr size_t sampleRate = 44100;
static constexpr float cutoff = 3800.f;


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


TEST_CASE("Low pass filter", "[AudioFramework:FirFilter]") {
	constexpr size_t numTaps = 255;
	const auto impulse = FirLowPassWindowed<float>(sampleRate, cutoff, numTaps, &HammingWindow<float, TIME_DOMAIN>);
	REQUIRE(impulse.Size() == numTaps);

	// Generate two signals just above and just below the cutoff and see their attenuation.
	const auto passSignal = GenTestSignal(sampleRate, cutoff * 0.85f);
	const auto rejectSignal = GenTestSignal(sampleRate, cutoff * 1.15f);

	const auto filteredPassSignal = ConvolutionFast(passSignal, impulse, convolution::full);
	const auto filteredRejectSignal = ConvolutionFast(rejectSignal, impulse, convolution::full);

	const float energyPass = SumSquare(passSignal);
	const float energyReject = SumSquare(rejectSignal);
	const float energyFilteredPass = SumSquare(filteredPassSignal);
	const float energyFilteredReject = SumSquare(filteredRejectSignal);

	REQUIRE(energyFilteredPass / energyPass > 0.85f);
	REQUIRE(energyFilteredReject / energyReject < 0.15f);
}


TEST_CASE("Arbitrary filter", "[AudioFramework:FirFilter]") {
	constexpr size_t numTaps = 255;
	constexpr std::array amplitudes = { 1.0f, 0.2f, 0.6f, 1.2f };
	constexpr float nyquistLimit = float(sampleRate / 2);
	constexpr std::array frequencies = { 0.0625f, 0.1875f, 0.375f, 0.75f };
	
	Spectrum<float> frequencyResponse(8192, amplitudes[0]);
	std::fill(frequencyResponse.begin() + 1024, frequencyResponse.begin() + 2048, amplitudes[1]);
	std::fill(frequencyResponse.begin() + 2048, frequencyResponse.begin() + 4096, amplitudes[2]);
	std::fill(frequencyResponse.begin() + 4096, frequencyResponse.end(), amplitudes[3]);

	const auto impulse = FirGeneralWindowed(frequencyResponse, numTaps);
	REQUIRE(impulse.Size() == numTaps);

	const float quality = FirAccuracy(impulse, frequencyResponse);
	REQUIRE(quality > 0.9f);
	
	for (size_t i = 0; i < amplitudes.size(); ++i) {
		const auto signal = GenTestSignal(sampleRate, frequencies[i]*nyquistLimit);
		const auto filtered = ConvolutionFast(signal, impulse, convolution::full);
		const float energy = std::sqrt(SumSquare(signal));
		const float filteredEnergy = std::sqrt(SumSquare(filtered));
		REQUIRE(std::abs(amplitudes[i]*energy - filteredEnergy)/energy < 0.1);
	}

}