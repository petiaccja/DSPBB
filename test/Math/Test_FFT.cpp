#include "../TestUtils.hpp"

#include <dspbb/Generators/Waveforms.hpp>
#include <dspbb/Math/DotProduct.hpp>
#include <dspbb/Math/FFT.hpp>
#include <dspbb/Math/Functions.hpp>
#include <dspbb/Math/Statistics.hpp>
#include <dspbb/Primitives/Signal.hpp>

#include <algorithm>
#include <catch2/catch.hpp>

using namespace dspbb;

constexpr uint64_t sampleRate = 16000;
constexpr float frequency = 2000;
constexpr size_t fftSize = 1024;



TEST_CASE("FFT - Real spectral peak", "[FFT]") {
	const auto signal = SineWave<float, TIME_DOMAIN>(fftSize, sampleRate, frequency);

	Spectrum<std::complex<float>> complexSpectrum = Fft(signal, FFT_FULL);
	Spectrum<float> powerSpectrum = Abs(complexSpectrum);

	REQUIRE(complexSpectrum.size() == fftSize);

	auto maxPos = std::max_element(powerSpectrum.begin(), powerSpectrum.begin() + powerSpectrum.size() / 2) - powerSpectrum.begin();
	auto maxPosExpected = FourierFrequency2Bin(frequency, fftSize, sampleRate);
	REQUIRE(std::abs(intptr_t(maxPos) - intptr_t(maxPosExpected)) <= 1);
}


TEST_CASE("FFT - Complex spectral peak", "[FFT]") {
	const auto signal = SineWave<std::complex<float>, TIME_DOMAIN>(fftSize, sampleRate, frequency);

	Spectrum<std::complex<float>> complexSpectrum = Fft(signal);
	Spectrum<float> powerSpectrum = Abs(complexSpectrum);

	REQUIRE(complexSpectrum.size() == fftSize);

	auto maxPos = std::max_element(powerSpectrum.begin(), powerSpectrum.begin() + powerSpectrum.size() / 2) - powerSpectrum.begin();
	auto maxPosExpected = FourierFrequency2Bin(frequency, fftSize, sampleRate);
	REQUIRE(std::abs(intptr_t(maxPos) - intptr_t(maxPosExpected)) <= 1);
}


TEST_CASE("IFFT - Real identity", "[FFT]") {
	const auto signal = SineWave<float, TIME_DOMAIN>(fftSize, sampleRate, frequency);
	Spectrum<std::complex<float>> spectrum = Fft(signal, FFT_HALF);
	const auto repro = Ifft(spectrum, FFT_HALF, signal.size() % 2 == 0);

	REQUIRE(Max(Abs(signal - repro)) < 1e-4f);
}


TEST_CASE("IFFT - Complex identity", "[FFT]") {
	const auto signal = RandomSignal<std::complex<float>, TIME_DOMAIN>(fftSize);
	Spectrum<std::complex<float>> spectrum = Fft(signal);
	const auto repro = Ifft(spectrum);

	REQUIRE(Max(Abs(signal - repro)) < 1e-4f);
}


TEST_CASE("Parseval's relation", "[FFT]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(fftSize);
	Spectrum<std::complex<float>> spectrum = Fft(signal, FFT_FULL);

	const float signalSum = SumSquare(signal);
	const float spectrumSum = SumSquare(Abs(spectrum));

	REQUIRE(signalSum == Approx(spectrumSum / fftSize));
}


TEST_CASE("FFT - Full real even", "[FFT]") {
	Signal<float> even(64, 0.0f);
	even[30] = 1.0f;
	const Spectrum<std::complex<float>> evenHalf = Fft(even, FFT_HALF);
	const Spectrum<std::complex<float>> evenFull = Fft(even, FFT_FULL);
	REQUIRE(evenHalf.size() == 33);
	REQUIRE(evenFull.size() == 64);
	REQUIRE(std::all_of(evenHalf.begin(), evenHalf.end(), [](const auto& v) { return std::abs(v) == Approx(1.0f); }));
	REQUIRE(std::all_of(evenFull.begin(), evenFull.end(), [](const auto& v) { return std::abs(v) == Approx(1.0f); }));
	const auto normal = Spectrum<std::complex<float>>(evenFull.begin() + 1, evenFull.begin() + 32);
	const auto conjugate = Spectrum<std::complex<float>>(evenFull.rbegin(), evenFull.rbegin() + 31);
	REQUIRE(Max(Abs(normal - Conj(conjugate))) < 1e-4f);
}

TEST_CASE("FFT - Full real odd", "[FFT]") {
	Signal<float> odd(63, 0.0f);
	odd[30] = 1.0f;
	Spectrum<std::complex<float>> evenHalf = Fft(odd, FFT_HALF);
	Spectrum<std::complex<float>> evenFull = Fft(odd, FFT_FULL);
	REQUIRE(evenHalf.size() == 32);
	REQUIRE(evenFull.size() == 63);
	REQUIRE(std::all_of(evenHalf.begin(), evenHalf.end(), [](const auto& v) { return std::abs(v) == Approx(1.0f); }));
	REQUIRE(std::all_of(evenFull.begin(), evenFull.end(), [](const auto& v) { return std::abs(v) == Approx(1.0f); }));
	const auto normal = Spectrum<std::complex<float>>(evenFull.begin() + 1, evenFull.begin() + 32);
	const auto conjugate = Spectrum<std::complex<float>>(evenFull.rbegin(), evenFull.rbegin() + 31);
	REQUIRE(Max(Abs(normal - Conj(conjugate))) < 1e-4f);
}


TEST_CASE("FFT - Full real identity", "[FFT]") {
	const std::array<size_t, 7> sizes = { 63, 64, 65, 66, 67, 68, 69 };

	for (auto s : sizes) {
		const auto signal = RandomSignal<float, TIME_DOMAIN>(s);
		const Spectrum<std::complex<float>> spectrum = Fft(signal, FFT_FULL);
		const Signal<float> repro = Ifft(spectrum, FFT_FULL);
		REQUIRE(signal.size() == repro.size());
		REQUIRE(Max(Abs(signal - repro)) < 0.001f);
	}
}

TEST_CASE("FFT - Half real identity", "[FFT]") {
	const std::array<size_t, 7> sizes = { 63, 64, 65, 66, 67, 68, 69 };

	for (auto s : sizes) {
		const auto signal = RandomSignal<float, TIME_DOMAIN>(s);
		const Spectrum<std::complex<float>> spectrum = Fft(signal, FFT_HALF);
		const Signal<float> repro = Ifft(spectrum, FFT_HALF, signal.size() % 2 == 0);
		REQUIRE(signal.size() == repro.size());
		REQUIRE(Max(Abs(signal - repro)) < 0.001f);
	}
}

TEST_CASE("FFT shift even", "[FFT]") {
	const Spectrum<float> s = { 0, 1, 2, 3, 4, 5 };
	const Spectrum<float> e = { 3, 4, 5, 0, 1, 2 };
	const auto r = FftShift(s);
	REQUIRE(Max(r - e) == 0);
}

TEST_CASE("FFT shift odd", "[FFT]") {
	const Spectrum<float> s = { 0, 1, 2, 3, 4, 5, 6 };
	const Spectrum<float> e = { 4, 5, 6, 0, 1, 2, 3 };
	const auto r = FftShift(s);
	REQUIRE(Max(r - e) == 0);
}

TEST_CASE("FFT shift 1", "[FFT]") {
	const Spectrum<float> s = { 0 };
	const auto r = FftShift(s);
	REQUIRE(r[0] == 0);
}

TEST_CASE("FFT shift empty", "[FFT]") {
	const Spectrum<float> s = {};
	const auto r = FftShift(s);
	REQUIRE(r.size() == 0);
}

TEST_CASE("FFT inverse shift even", "[FFT]") {
	const Spectrum<float> s = { 3, 4, 5, 0, 1, 2 };
	const Spectrum<float> e = { 0, 1, 2, 3, 4, 5 };
	const auto r = IfftShift(s);
	REQUIRE(Max(r - e) == 0);
}

TEST_CASE("FFT inverse shift odd", "[FFT]") {
	const Spectrum<float> s = { 4, 5, 6, 0, 1, 2, 3 };
	const Spectrum<float> e = { 0, 1, 2, 3, 4, 5, 6 };
	const auto r = IfftShift(s);
	REQUIRE(Max(r - e) == 0);
}