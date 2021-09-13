#include <algorithm>
#include <catch2/catch.hpp>
#include <dspbb/Math/FFT.hpp>
#include <dspbb/Generators/Waveforms.hpp>
#include <dspbb/Math/DotProduct.hpp>
#include <dspbb/Math/Functions.hpp>
#include <dspbb/Math/Statistics.hpp>
#include <dspbb/Primitives/Signal.hpp>

using namespace dspbb;

constexpr uint64_t sampleRate = 16000;
constexpr float frequency = 2000;
constexpr size_t fftSize = 1024;



TEST_CASE("FFT - Real spectral peak", "[FFT]") {
	const auto signal = SineWave<float, TIME_DOMAIN>(fftSize, sampleRate, frequency);

	Spectrum<std::complex<float>> complexSpectrum = FourierTransform(signal, true);
	Spectrum<float> powerSpectrum = Abs(complexSpectrum);

	REQUIRE(complexSpectrum.Length() == fftSize);

	auto maxPos = std::max_element(powerSpectrum.begin(), powerSpectrum.begin() + powerSpectrum.Size() / 2) - powerSpectrum.begin();
	auto maxPosExpected = FourierFrequency2Bin(frequency, fftSize, sampleRate);
	REQUIRE(std::abs(intptr_t(maxPos) - intptr_t(maxPosExpected)) <= 1);
}


TEST_CASE("FFT - Complex spectral peak", "[FFT]") {
	const auto signal = SineWave<std::complex<float>, TIME_DOMAIN>(fftSize, sampleRate, frequency);

	Spectrum<std::complex<float>> complexSpectrum = FourierTransform(signal);
	Spectrum<float> powerSpectrum = Abs(complexSpectrum);

	REQUIRE(complexSpectrum.Length() == fftSize);

	auto maxPos = std::max_element(powerSpectrum.begin(), powerSpectrum.begin() + powerSpectrum.Size() / 2) - powerSpectrum.begin();
	auto maxPosExpected = FourierFrequency2Bin(frequency, fftSize, sampleRate);
	REQUIRE(std::abs(intptr_t(maxPos) - intptr_t(maxPosExpected)) <= 1);
}


TEST_CASE("IFFT - Real identity", "[FFT]") {
	const auto signal = SineWave<float, TIME_DOMAIN>(fftSize, sampleRate, frequency);
	Spectrum<std::complex<float>> spectrum = FourierTransform(signal, false);
	const auto repro = InverseFourierTransformR(spectrum, signal.Size());

	const float norm = Norm(signal);
	const float rnorm = Norm(repro);
	const float similarity = DotProduct(signal, repro) / norm / rnorm;
	REQUIRE(similarity == Approx(1));
	REQUIRE(norm == Approx(rnorm));
}


TEST_CASE("IFFT - Complex identity", "[FFT]") {
	const auto signal = SineWave<std::complex<float>, TIME_DOMAIN>(fftSize, sampleRate, frequency);
	Spectrum<std::complex<float>> spectrum = FourierTransform(signal);
	const auto repro = InverseFourierTransformC(spectrum);

	const float norm = std::abs(Norm(signal));
	const float rnorm = std::abs(Norm(repro));
	const float similarity = std::abs(DotProduct(signal, repro)) / norm / rnorm;
	REQUIRE(similarity == Approx(1));
	REQUIRE(norm == Approx(rnorm));
}


TEST_CASE("Parseval's relation", "[FFT]") {
	const auto signal = SineWave<float, TIME_DOMAIN>(fftSize, sampleRate, frequency);
	Spectrum<std::complex<float>> spectrum = FourierTransform(signal, true);

	const float signalSum = SumSquare(signal);
	const float spectrumSum = SumSquare(Abs(spectrum));

	REQUIRE(signalSum == Approx(spectrumSum / fftSize));
}


TEST_CASE("FFT - Full real even", "[FFT]") {
	TimeSignal<float> even(64, 0.0f);
	even[30] = 1.0f;
	const Spectrum<std::complex<float>> evenHalf = FourierTransform(even, false);
	const Spectrum<std::complex<float>> evenFull = FourierTransform(even, true);
	REQUIRE(evenHalf.Length() == 33);
	REQUIRE(evenFull.Length() == 64);
	REQUIRE(std::all_of(evenHalf.begin(), evenHalf.end(), [](const auto& v) { return std::abs(v) == Approx(1.0f); }));
	REQUIRE(std::all_of(evenFull.begin(), evenFull.end(), [](const auto& v) { return std::abs(v) == Approx(1.0f); }));
}

TEST_CASE("FFT - Full real odd", "[FFT]") {
	TimeSignal<float> odd(63, 0.0f);
	odd[30] = 1.0f;
	Spectrum<std::complex<float>> evenHalf = FourierTransform(odd, false);
	Spectrum<std::complex<float>> evenFull = FourierTransform(odd, true);
	REQUIRE(evenHalf.Length() == 32);
	REQUIRE(evenFull.Length() == 63);
	REQUIRE(std::all_of(evenHalf.begin(), evenHalf.end(), [](const auto& v) { return std::abs(v) == Approx(1.0f); }));
	REQUIRE(std::all_of(evenFull.begin(), evenFull.end(), [](const auto& v) { return std::abs(v) == Approx(1.0f); }));
}


TEST_CASE("FFT - Full real identity", "[FFT]") {
	const std::array<size_t, 7> sizes = { 63, 64, 65, 66, 67, 68, 69 };

	for (auto s : sizes) {
		TimeSignal<float> signal(s, 0.0f);
		signal[signal.Size() / 2] = 1.0f;
		const Spectrum<std::complex<float>> spectrum = FourierTransform(signal, true);
		const TimeSignal<float> repro = InverseFourierTransformR(spectrum, 0);
		REQUIRE(signal.Size() == repro.Size());
		REQUIRE(Max(Abs(signal - repro)) < 0.001f);
	}
}

TEST_CASE("FFT - Half real identity", "[FFT]") {
	const std::array<size_t, 7> sizes = { 63, 64, 65, 66, 67, 68, 69 };

	for (auto s : sizes) {
		TimeSignal<float> signal(s, 0.0f);
		signal[signal.Size() / 2] = 1.0f;
		const Spectrum<std::complex<float>> spectrum = FourierTransform(signal, false);
		const TimeSignal<float> repro = InverseFourierTransformR(spectrum, signal.Size());
		REQUIRE(signal.Size() == repro.Size());
		REQUIRE(Max(Abs(signal - repro)) < 0.001f);
	}
}
