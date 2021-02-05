#include "dspbb/Generators/Sine.hpp"
#include "dspbb/Math/DotProduct.hpp"
#include "dspbb/Math/Statistics.hpp"


#include <Catch2/catch.hpp>
#include <algorithm>
#include <dspbb/DSP/FFT.hpp>
#include <dspbb/Math/Math.hpp>
#include <dspbb/Primitives/Signal.hpp>

using namespace dspbb;

constexpr uint64_t sampleRate = 16000;
constexpr float frequency = 2000;
constexpr size_t fftSize = 1024;



TEST_CASE("FFT - Real spectral peak", "[AudioFramework:FFT]") {
	const auto signal = SineWave<float, TIME_DOMAIN>(fftSize, sampleRate, frequency);

	Spectrum<std::complex<float>> complexSpectrum = FourierTransform(signal);
	Spectrum<float> powerSpectrum = abs(complexSpectrum);

	REQUIRE(complexSpectrum.Length() == fftSize);

	auto maxPos = std::max_element(powerSpectrum.begin(), powerSpectrum.begin() + powerSpectrum.Size() / 2) - powerSpectrum.begin();
	auto maxPosExpected = FourierFrequency2Bin(frequency, fftSize, sampleRate);
	REQUIRE(std::abs(intptr_t(maxPos) - intptr_t(maxPosExpected)) <= 1);
}


TEST_CASE("FFT - Complex spectral peak", "[AudioFramework:FFT]") {
	const auto signal = SineWave<std::complex<float>, TIME_DOMAIN>(fftSize, sampleRate, frequency);

	Spectrum<std::complex<float>> complexSpectrum = FourierTransform(signal);
	Spectrum<float> powerSpectrum = abs(complexSpectrum);

	REQUIRE(complexSpectrum.Length() == fftSize);

	auto maxPos = std::max_element(powerSpectrum.begin(), powerSpectrum.begin() + powerSpectrum.Size() / 2) - powerSpectrum.begin();
	auto maxPosExpected = FourierFrequency2Bin(frequency, fftSize, sampleRate);
	REQUIRE(std::abs(intptr_t(maxPos) - intptr_t(maxPosExpected)) <= 1);
}


TEST_CASE("IFFT - Real identity", "[AudioFramework:FFT]") {
	const auto signal = SineWave<float, TIME_DOMAIN>(fftSize, sampleRate, frequency);
	Spectrum<std::complex<float>> spectrum = FourierTransform(signal);
	const auto repro = FourierTransform<float>(spectrum);

	const float norm = Norm(signal);
	const float rnorm = Norm(repro);
	const float similarity = DotProduct(signal, repro, signal.Length()) / norm / rnorm;
	REQUIRE(similarity == Approx(1));
	REQUIRE(norm == Approx(rnorm));
}


TEST_CASE("IFFT - Complex identity", "[AudioFramework:FFT]") {
	const auto signal = SineWave<std::complex<float>, TIME_DOMAIN>(fftSize, sampleRate, frequency);
	Spectrum<std::complex<float>> spectrum = FourierTransform(signal);
	const auto repro = FourierTransform<std::complex<float>>(spectrum);

	const float norm = std::abs(Norm(signal));
	const float rnorm = std::abs(Norm(repro));
	const float similarity = std::abs(DotProduct(signal, repro, signal.Length())) / norm / rnorm;
	REQUIRE(similarity == Approx(1));
	REQUIRE(norm == Approx(rnorm));
}


TEST_CASE("Parseval's relation", "[AudioFramework:FFT]") {
	const auto signal = SineWave<float, TIME_DOMAIN>(fftSize, sampleRate, frequency);
	Spectrum<std::complex<float>> spectrum = FourierTransform(signal);

	const float signalSum = SumSquare(signal);
	const float spectrumSum = SumSquare(abs(spectrum));
	
	REQUIRE(signalSum == Approx(spectrumSum / fftSize));
}