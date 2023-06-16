#include <dspbb/Filtering/MeasureFilter.hpp>
#include <dspbb/Filtering/Windowing.hpp>
#include <dspbb/Math/Functions.hpp>
#include <dspbb/Math/Statistics.hpp>
#include <dspbb/Primitives/Signal.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>


using namespace dspbb;
using Catch::Approx;


TEST_CASE("Coherent gain", "[WindowFunctions]") {
	BasicSignal<float, TIME_DOMAIN> window(32, 0.5f);
	REQUIRE(CoherentGain(window) == Approx(0.5f));
}

TEST_CASE("Energy gain", "[WindowFunctions]") {
	BasicSignal<float, TIME_DOMAIN> window(32, 0.5f);
	REQUIRE(EnergyGain(window) == Approx(0.25f));
}

template <class T, eSignalDomain Domain>
static bool IsSymmetric(const BasicSignal<T, Domain>& window) {
	auto fw = window.begin();
	auto rev = window.rbegin();
	for (; fw <= rev.base(); ++fw, ++rev) {
		auto diff = std::abs(*fw - *rev);
		if (diff > 0.001f) {
			return false;
		}
	}
	return true;
}

template <class T, eSignalDomain Domain>
static bool IsPeakCentered(const BasicSignal<T, Domain>& window) {
	return std::abs(Max(Abs(window)) - std::abs(window[window.size() / 2])) < 0.01f;
}

TEST_CASE("Hamming window", "[WindowFunctions]") {
	auto window = windows::hamming.operator()<float>(256);

	REQUIRE(window.size() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(CoherentGain(window) == Approx(0.54f).margin(0.01f));
}


TEST_CASE("Hamming window complex", "[WindowFunctions]") {
	Signal<std::complex<float>> window(256);
	windows::hamming(window);

	REQUIRE(window.size() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(std::abs(CoherentGain(window)) == Approx(0.54f).margin(0.01f));
	REQUIRE(Sum(Abs(Imag(window))) == Approx(0.0f));
}


TEST_CASE("Flat top window", "[WindowFunctions]") {
	auto window = windows::flattop.operator()<float>(256);

	REQUIRE(window.size() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(CoherentGain(window) == Approx(0.22f).margin(0.01f));
}


TEST_CASE("Flat top complex", "[WindowFunctions]") {
	Signal<std::complex<float>> window(256);
	windows::flattop(window);

	REQUIRE(window.size() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(std::abs(CoherentGain(window)) == Approx(0.22f).margin(0.01f));
	REQUIRE(Sum(Abs(Imag(window))) == Approx(0.0f));
}


TEST_CASE("Rectangular window", "[WindowFunctions]") {
	auto window = windows::rectangular.operator()<float>(256);

	REQUIRE(window.size() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(CoherentGain(window) == Approx(1.0f).margin(0.01f));
}


TEST_CASE("Rectangular complex", "[WindowFunctions]") {
	Signal<std::complex<float>> window(256);
	windows::rectangular(window);

	REQUIRE(window.size() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(std::abs(CoherentGain(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(Sum(Abs(Imag(window))) == Approx(0.0f));
}


TEST_CASE("Triangular window", "[WindowFunctions]") {
	auto window = windows::triangular.operator()<float>(256);

	REQUIRE(window.size() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(CoherentGain(window) == Approx(0.50f).margin(0.01f));
}


TEST_CASE("Triangular complex", "[WindowFunctions]") {
	Signal<std::complex<float>> window(256);
	windows::triangular(window);

	REQUIRE(window.size() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(std::abs(CoherentGain(window)) == Approx(0.50f).margin(0.01f));
	REQUIRE(Sum(Abs(Imag(window))) == Approx(0.0f));
}

TEST_CASE("Blackman window", "[WindowFunctions]") {
	auto window = windows::blackman.operator()<float>(256);

	REQUIRE(window.size() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(CoherentGain(window) == Approx(0.42f).margin(0.01f));
}


TEST_CASE("Blackman complex", "[WindowFunctions]") {
	Signal<std::complex<float>> window(256);
	windows::blackman(window);

	REQUIRE(window.size() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(std::abs(CoherentGain(window)) == Approx(0.42f).margin(0.01f));
	REQUIRE(Sum(Abs(Imag(window))) == Approx(0.0f));
}

TEST_CASE("Blackman-Harris window", "[WindowFunctions]") {
	auto window = windows::blackmanHarris.operator()<float>(256);

	REQUIRE(window.size() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(CoherentGain(window) == Approx(0.36f).margin(0.01f));
}


TEST_CASE("Blackman-Harris complex", "[WindowFunctions]") {
	Signal<std::complex<float>> window(256);
	windows::blackmanHarris(window);

	REQUIRE(window.size() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(std::abs(CoherentGain(window)) == Approx(0.36f).margin(0.01f));
	REQUIRE(Sum(Abs(Imag(window))) == Approx(0.0f));
}

TEST_CASE("Gaussian window", "[WindowFunctions]") {
	auto window = windows::gaussian.sigma(0.3f).operator()<float>(256);

	REQUIRE(window.size() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(CoherentGain(window) == Approx(0.37f).margin(0.01f));
}


TEST_CASE("Gaussian complex", "[WindowFunctions]") {
	Signal<std::complex<float>> window(256);
	windows::gaussian.sigma(0.3f)(window);

	REQUIRE(window.size() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(std::abs(CoherentGain(window)) == Approx(0.37f).margin(0.01f));
	REQUIRE(Sum(Abs(Imag(window))) == Approx(0.0f));
}

TEST_CASE("Kaiser window", "[WindowFunctions]") {
	auto window = windows::kaiser.alpha(1.0f).operator()<float>(256);

	REQUIRE(window.size() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(CoherentGain(window) == Approx(0.67f).margin(0.01f));
}

TEST_CASE("Kaiser complex", "[WindowFunctions]") {
	Signal<std::complex<float>> window(256);
	windows::kaiser.alpha(0.5f)(window);

	REQUIRE(window.size() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(std::abs(CoherentGain(window)) == Approx(0.85f).margin(0.01f));
	REQUIRE(Sum(Abs(Imag(window))) == Approx(0.0f));
}

TEST_CASE("Lanczos window", "[WindowFunctions]") {
	auto window = windows::lanczos.operator()<float>(255);

	REQUIRE(window.size() == 255);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(CoherentGain(window) == Approx(0.59f).margin(0.01f));
}

TEST_CASE("Lanczos complex", "[WindowFunctions]") {
	Signal<std::complex<float>> window(255);
	windows::lanczos(window);

	REQUIRE(window.size() == 255);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(std::abs(CoherentGain(window)) == Approx(0.59f).margin(0.01f));
	REQUIRE(Sum(Abs(Imag(window))) == Approx(0.0f));
}

TEST_CASE("Dolph-Chebyshev window", "[WindowFunctions]") {
	auto window = windows::dolphChebyshev.attenuation(0.01f).operator()<float>(255);

	REQUIRE(window.size() == 255);
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));

	const auto [amplitude, phase] = FrequencyResponse(window, 2048);
	const auto atten = Max(AsView(amplitude).subsignal(200)) / amplitude[0];
	REQUIRE(atten == Approx(0.01f).epsilon(1e-3f));
}

TEST_CASE("Dolph-Chebyshev complex", "[WindowFunctions]") {
	Signal<std::complex<float>> window(256);
	windows::dolphChebyshev.attenuation(0.001f)(window);
	auto dbg = Real(window);

	REQUIRE(window.size() == 256);
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(Sum(Abs(Imag(window))) == Approx(0.0f));

	const auto [amplitude, phase] = FrequencyResponse(Real(window), 2048);
	const auto atten = Max(AsView(amplitude).subsignal(200)) / amplitude[0];
	REQUIRE(atten == Approx(0.001f).epsilon(1e-3f));
}
