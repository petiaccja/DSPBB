#include <catch2/catch.hpp>
#include <dspbb/Filtering/FilterParameters.hpp>
#include <dspbb/Generators/Spaces.hpp>
#include <random>

using namespace dspbb;



//------------------------------------------------------------------------------
// Helpers to define example responses
//------------------------------------------------------------------------------

float Transition(float x, float from, float to) {
	const auto x_ = (x - from) / (to - from);
	return 3 * x_ * x_ - 2 * x_ * x_ * x_;
}

float Passband() {
	return 1.0f;
}

float Stopband() {
	return 0.0f;
}

float Ripple(float x, float scale, float amplitude, float limit = std::numeric_limits<float>::infinity()) {
	if (x > limit) {
		return 0.0f;
	}
	const float p = 1.655f / scale;
	const float a = amplitude * 2.325f;
	const float pxs = p * x + 1;
	if (pxs >= 0.0f) {
		return a * std::sin(pxs) * std::pow(pxs, 3) / (std::pow(pxs, 4) + 3.0f);
	}
	return 0.0f;
}

struct Band {
	float lower = 0.0f;
	float upper = 1.0f;
	bool pass = true;
	float ripple = 0.0f;
};

float Response(float x, const std::vector<Band>& bands) {
	const auto withinIt = std::find_if(bands.begin(), bands.end(), [&x](const Band& b) {
		return b.lower <= x && x <= b.upper;
	});
	const auto betweenIt = std::adjacent_find(bands.begin(), bands.end(), [&x](const Band& b1, const Band& b2) {
		return b1.upper < x && x < b2.lower;
	});

	float y = 0.0f;
	if (withinIt != bands.end()) {
		y += float(withinIt->pass);
	}
	if (betweenIt != bands.end()) {
		auto nextIt = betweenIt;
		++nextIt;
		float from = betweenIt->upper;
		float to = nextIt->lower;
		if (betweenIt->pass) {
			std::swap(from, to);
		}
		y += Transition(x, from, to);
	}
	for (auto& b : bands) {
		const float xl = x - b.lower;
		const float xu = b.upper - x;
		const float w = b.upper - b.lower;
		const float sign = 2.0f * float(b.pass) - 1.0f;
		y += sign * Ripple(xl, 0.06f * w, b.ripple, w / 2.0f);
		y += sign * Ripple(xu, 0.06f * w, b.ripple, w / 2.0f);
	}
	return std::abs(y);
}


//------------------------------------------------------------------------------
// Define example responses
//------------------------------------------------------------------------------

constexpr float transitionLower = 0.35f;
constexpr float transitionUpper = 0.45f;
constexpr float ripplePass = 0.05f;
constexpr float rippleStop = 0.03f;

const std::vector lowpassFlat = {
	Band{ 0.0f, transitionLower, true },
	Band{ transitionUpper, 1.0f, false },
};

const std::vector lowpassRipple = {
	Band{ 0.0f, transitionLower, true, ripplePass },
	Band{ transitionUpper, 1.0f, false, rippleStop },
};

const std::vector highpassFlat = {
	Band{ 0.0f, transitionLower, false },
	Band{ transitionUpper, 1.0f, true },
};

const std::vector highpassRipple = {
	Band{ 0.0f, transitionLower, false, rippleStop },
	Band{ transitionUpper, 1.0f, true, ripplePass },
};

constexpr float transitionLower1 = 0.25f;
constexpr float transitionUpper1 = 0.35f;
constexpr float transitionLower2 = 0.55f;
constexpr float transitionUpper2 = 0.65f;
constexpr float ripplePass1 = 0.06f;
constexpr float ripplePass2 = 0.05f;
constexpr float rippleStop1 = 0.03f;
constexpr float rippleStop2 = 0.04f;

const std::vector bandpassFlat = {
	Band{ 0.0f, transitionLower1, false },
	Band{ transitionUpper1, transitionLower2, true },
	Band{ transitionUpper2, 1.0f, false },
};

const std::vector bandpassRipple = {
	Band{ 0.0f, transitionLower1, false, rippleStop1 },
	Band{ transitionUpper1, transitionLower2, true, ripplePass1 },
	Band{ transitionUpper2, 1.0f, false, rippleStop2 },
};

const std::vector bandstopFlat = {
	Band{ 0.0f, transitionLower1, true },
	Band{ transitionUpper1, transitionLower2, false },
	Band{ transitionUpper2, 1.0f, true },
};

const std::vector bandstopRipple = {
	Band{ 0.0f, transitionLower1, true, ripplePass1 },
	Band{ transitionUpper1, transitionLower2, false, rippleStop1 },
	Band{ transitionUpper2, 1.0f, true, ripplePass2 },
};


//------------------------------------------------------------------------------
// Helpers for tests
//------------------------------------------------------------------------------

template <class Func>
Spectrum<float> MockSpectrum(size_t size, Func func) {
	auto spectrum = LinSpace<float, FREQUENCY_DOMAIN>(0.0f, 1.0f, size, true);
	std::for_each(spectrum.begin(), spectrum.end(), [&func](auto& v) { v = func(v); });
	return spectrum;
}

Spectrum<float> MockSpectrum(size_t size, const std::vector<Band>& bands) {
	return MockSpectrum(size, [&](float x) { return Response(x, bands); });
}

//------------------------------------------------------------------------------
// Verify classification
//------------------------------------------------------------------------------

TEST_CASE("Classify flat low-pass", "[FilterParameters]") {
	const auto response = MockSpectrum(1000, lowpassFlat);
	REQUIRE_NOTHROW(ParametrizeLowpassFilter(response));
	REQUIRE_THROWS(ParametrizeHighpassFilter(response));
	REQUIRE_THROWS(ParametrizeBandpassFilter(response));
	REQUIRE_THROWS(ParametrizeBandstopFilter(response));
}

TEST_CASE("Classify ripple low-pass", "[FilterParameters]") {
	const auto response = MockSpectrum(1000, lowpassRipple);
	REQUIRE_NOTHROW(ParametrizeLowpassFilter(response));
	REQUIRE_THROWS(ParametrizeHighpassFilter(response));
	REQUIRE_THROWS(ParametrizeBandpassFilter(response));
	REQUIRE_THROWS(ParametrizeBandstopFilter(response));
}


TEST_CASE("Classify flat high-pass", "[FilterParameters]") {
	const auto response = MockSpectrum(1000, highpassFlat);
	REQUIRE_THROWS(ParametrizeLowpassFilter(response));
	REQUIRE_NOTHROW(ParametrizeHighpassFilter(response));
	REQUIRE_THROWS(ParametrizeBandpassFilter(response));
	REQUIRE_THROWS(ParametrizeBandstopFilter(response));
}

TEST_CASE("Classify ripple high-pass", "[FilterParameters]") {
	const auto response = MockSpectrum(1000, highpassRipple);
	REQUIRE_THROWS(ParametrizeLowpassFilter(response));
	REQUIRE_NOTHROW(ParametrizeHighpassFilter(response));
	REQUIRE_THROWS(ParametrizeBandpassFilter(response));
	REQUIRE_THROWS(ParametrizeBandstopFilter(response));
}


TEST_CASE("Classify flat band-pass", "[FilterParameters]") {
	const auto response = MockSpectrum(1000, bandpassFlat);
	REQUIRE_THROWS(ParametrizeLowpassFilter(response));
	REQUIRE_THROWS(ParametrizeHighpassFilter(response));
	REQUIRE_NOTHROW(ParametrizeBandpassFilter(response));
	REQUIRE_THROWS(ParametrizeBandstopFilter(response));
}

TEST_CASE("Classify ripple band-pass", "[FilterParameters]") {
	const auto response = MockSpectrum(1000, bandpassRipple);
	REQUIRE_THROWS(ParametrizeLowpassFilter(response));
	REQUIRE_THROWS(ParametrizeHighpassFilter(response));
	REQUIRE_NOTHROW(ParametrizeBandpassFilter(response));
	REQUIRE_THROWS(ParametrizeBandstopFilter(response));
}


TEST_CASE("Classify flat band-stop", "[FilterParameters]") {
	const auto response = MockSpectrum(1000, bandstopFlat);
	REQUIRE_THROWS(ParametrizeLowpassFilter(response));
	REQUIRE_THROWS(ParametrizeHighpassFilter(response));
	REQUIRE_THROWS(ParametrizeBandpassFilter(response));
	REQUIRE_NOTHROW(ParametrizeBandstopFilter(response));
}

TEST_CASE("Classify ripple band-stop", "[FilterParameters]") {
	const auto response = MockSpectrum(1000, bandstopRipple);
	REQUIRE_THROWS(ParametrizeLowpassFilter(response));
	REQUIRE_THROWS(ParametrizeHighpassFilter(response));
	REQUIRE_THROWS(ParametrizeBandpassFilter(response));
	REQUIRE_NOTHROW(ParametrizeBandstopFilter(response));
}


//------------------------------------------------------------------------------
// Verify parametrizations
//------------------------------------------------------------------------------

TEST_CASE("Parametrize flat low-pass", "[FilterParameters]") {
	const auto response = MockSpectrum(1000, lowpassFlat);
	const auto params = ParametrizeLowpassFilter(response);
	REQUIRE(params.passbandEdge == Approx(transitionLower).margin(0.005f));
	REQUIRE(params.stopbandEdge == Approx(transitionUpper).margin(0.005f));
	REQUIRE(params.passbandRipple == Approx(0.0f).margin(1e-5f));
	REQUIRE(params.stopbandAtten == Approx(0.0f).margin(1e-5f));
}

TEST_CASE("Parametrize ripple low-pass", "[FilterParameters]") {
	const auto response = MockSpectrum(1000, lowpassRipple);
	const auto params = ParametrizeLowpassFilter(response);
	REQUIRE(params.passbandEdge == Approx(transitionLower).margin(0.005f));
	REQUIRE(params.stopbandEdge == Approx(transitionUpper).margin(0.005f));
	REQUIRE(params.passbandRipple == Approx(ripplePass).margin(1e-4f));
	REQUIRE(params.stopbandAtten == Approx(rippleStop).margin(1e-4f));
}


TEST_CASE("Parametrize flat high-pass", "[FilterParameters]") {
	const auto response = MockSpectrum(1000, highpassFlat);
	const auto params = ParametrizeHighpassFilter(response);
	REQUIRE(params.stopbandEdge == Approx(transitionLower).margin(0.005f));
	REQUIRE(params.passbandEdge == Approx(transitionUpper).margin(0.005f));
	REQUIRE(params.stopbandAtten == Approx(0.0f).margin(1e-5f));
	REQUIRE(params.passbandRipple == Approx(0.0f).margin(1e-5f));
}

TEST_CASE("Parametrize ripple high-pass", "[FilterParameters]") {
	const auto response = MockSpectrum(1000, highpassRipple);
	const auto params = ParametrizeHighpassFilter(response);
	REQUIRE(params.stopbandEdge == Approx(transitionLower).margin(0.005f));
	REQUIRE(params.passbandEdge == Approx(transitionUpper).margin(0.005f));
	REQUIRE(params.stopbandAtten == Approx(rippleStop).margin(1e-4f));
	REQUIRE(params.passbandRipple == Approx(ripplePass).margin(1e-4f));
}


TEST_CASE("Parametrize flat band-pass", "[FilterParameters]") {
	const auto response = MockSpectrum(1000, bandpassFlat);
	const auto params = ParametrizeBandpassFilter(response);
	REQUIRE(params.lowerStopbandEdge == Approx(transitionLower1).margin(0.005f));
	REQUIRE(params.passbandLowerEdge == Approx(transitionUpper1).margin(0.005f));
	REQUIRE(params.passbandUpperEdge == Approx(transitionLower2).margin(0.005f));
	REQUIRE(params.upperStopbandEdge == Approx(transitionUpper2).margin(0.005f));
	REQUIRE(params.lowerStopbandAtten == Approx(0.0f).margin(1e-5f));
	REQUIRE(params.passbandRipple == Approx(0.0f).margin(1e-5f));
	REQUIRE(params.upperStopbandAtten == Approx(0.0f).margin(1e-5f));
}

TEST_CASE("Parametrize ripple band-pass", "[FilterParameters]") {
	const auto response = MockSpectrum(1000, bandpassRipple);
	const auto params = ParametrizeBandpassFilter(response);
	REQUIRE(params.lowerStopbandEdge == Approx(transitionLower1).margin(0.005f));
	REQUIRE(params.passbandLowerEdge == Approx(transitionUpper1).margin(0.005f));
	REQUIRE(params.passbandUpperEdge == Approx(transitionLower2).margin(0.005f));
	REQUIRE(params.upperStopbandEdge == Approx(transitionUpper2).margin(0.005f));
	REQUIRE(params.lowerStopbandAtten == Approx(rippleStop1).margin(3e-4f));
	REQUIRE(params.passbandRipple == Approx(ripplePass1).margin(3e-4f));
	REQUIRE(params.upperStopbandAtten == Approx(rippleStop2).margin(3e-4f));
}


TEST_CASE("Parametrize flat band-stop", "[FilterParameters]") {
	const auto response = MockSpectrum(1000, bandstopFlat);
	const auto params = ParametrizeBandstopFilter(response);
	REQUIRE(params.lowerPassbandEdge == Approx(transitionLower1).margin(0.005f));
	REQUIRE(params.stopbandLowerEdge == Approx(transitionUpper1).margin(0.005f));
	REQUIRE(params.stopbandUpperEdge == Approx(transitionLower2).margin(0.005f));
	REQUIRE(params.upperPassbandEdge == Approx(transitionUpper2).margin(0.005f));
	REQUIRE(params.lowerPassbandRipple == Approx(0.0f).margin(1e-5f));
	REQUIRE(params.stopbandAtten == Approx(0.0f).margin(1e-5f));
	REQUIRE(params.upperPassbandRipple == Approx(0.0f).margin(1e-5f));
}

TEST_CASE("Parametrize ripple band-stop", "[FilterParameters]") {
	const auto response = MockSpectrum(1000, bandstopRipple);
	const auto params = ParametrizeBandstopFilter(response);
	REQUIRE(params.lowerPassbandEdge == Approx(transitionLower1).margin(0.005f));
	REQUIRE(params.stopbandLowerEdge == Approx(transitionUpper1).margin(0.005f));
	REQUIRE(params.stopbandUpperEdge == Approx(transitionLower2).margin(0.005f));
	REQUIRE(params.upperPassbandEdge == Approx(transitionUpper2).margin(0.005f));
	REQUIRE(params.lowerPassbandRipple == Approx(ripplePass1).margin(3e-4f));
	REQUIRE(params.stopbandAtten == Approx(rippleStop1).margin(3e-4f));
	REQUIRE(params.upperPassbandRipple == Approx(ripplePass2).margin(3e-4f));
}
