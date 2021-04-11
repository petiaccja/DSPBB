#include <catch2/catch.hpp>
#include <algorithm>
#include <dspbb/Filtering/WindowFunctions.hpp>
#include <dspbb/Math/Functions.hpp>
#include <dspbb/Math/Statistics.hpp>
#include <dspbb/Primitives/Signal.hpp>

using namespace dspbb;


TEST_CASE("Coherent gain", "[WindowFunctions]") {
	Signal<float, TIME_DOMAIN> window(32, 0.5f);
	REQUIRE(CoherentGain(window) == Approx(0.5f));
}

TEST_CASE("Energy gain", "[WindowFunctions]") {
	Signal<float, TIME_DOMAIN> window(32, 0.5f);
	REQUIRE(EnergyGain(window) == Approx(0.25f));
}

template <class T, eSignalDomain Domain>
static bool IsSymmetric(const Signal<T, Domain>& window) {
	auto it = window.begin();
	auto rev = window.rbegin();
	for (; it != rev.base(); ++it, ++rev) {
		auto diff = std::abs(*it - *rev);
		if (diff > 0.001f) {
			return false;
		}
	}
	return true;
}

template <class T, eSignalDomain Domain>
static bool IsPeakCentered(const Signal<T, Domain>& window) {
	return std::abs(Max(Abs(window)) - std::abs(window[window.Size() / 2])) < 0.01f;
}

TEST_CASE("Hamming window", "[WindowFunctions]") {
	auto window = HammingWindow<float>(256);

	REQUIRE(window.Length() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(CoherentGain(window) == Approx(0.54f).margin(0.01f));
}


TEST_CASE("Hamming window complex", "[WindowFunctions]") {
	auto window = HammingWindow<std::complex<float>>(256);

	REQUIRE(window.Length() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(std::abs(CoherentGain(window)) == Approx(0.54f).margin(0.01f));
	REQUIRE(Sum(Abs(Imag(window))) == Approx(0.0f));
}


TEST_CASE("Flat top window", "[WindowFunctions]") {
	auto window = FlatTopWindow<float>(256);

	REQUIRE(window.Length() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(CoherentGain(window) == Approx(0.22f).margin(0.01f));
}


TEST_CASE("Flat top complex", "[WindowFunctions]") {
	auto window = FlatTopWindow<std::complex<float>>(256);

	REQUIRE(window.Length() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(std::abs(CoherentGain(window)) == Approx(0.22f).margin(0.01f));
	REQUIRE(Sum(Abs(Imag(window))) == Approx(0.0f));
}


TEST_CASE("Rectangular window", "[WindowFunctions]") {
	auto window = RectangularWindow<float>(256);

	REQUIRE(window.Length() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(CoherentGain(window) == Approx(1.0f).margin(0.01f));
}


TEST_CASE("Rectangular complex", "[WindowFunctions]") {
	auto window = RectangularWindow<std::complex<float>>(256);

	REQUIRE(window.Length() == 256);
	REQUIRE(IsPeakCentered(window));
	REQUIRE(IsSymmetric(window));
	REQUIRE(Max(Abs(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(std::abs(CoherentGain(window)) == Approx(1.0f).margin(0.01f));
	REQUIRE(Sum(Abs(Imag(window))) == Approx(0.0f));
}