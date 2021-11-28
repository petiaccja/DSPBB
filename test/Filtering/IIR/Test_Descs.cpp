#include <dspbb/Filtering/IIR/Descs.hpp>

#include <catch2/catch.hpp>

using namespace dspbb;

//------------------------------------------------------------------------------
// Butterworth
//------------------------------------------------------------------------------

constexpr float butterCutoff = 0.3f;
constexpr float butterBandLow = 0.4f;
constexpr float butterBandHigh = 0.6f;

TEST_CASE("Low pass butterworth", "[IIR Descs]") {
	const auto desc = Lowpass(BUTTERWORTH).Cutoff(butterCutoff);
	REQUIRE(desc.cutoff == butterCutoff);
}

TEST_CASE("High pass butterworth", "[IIR Descs]") {
	const auto desc = Highpass(BUTTERWORTH).Cutoff(butterCutoff);
	REQUIRE(desc.cutoff == butterCutoff);
}

TEST_CASE("Band pass butterworth", "[IIR Descs]") {
	const auto desc = Bandpass(BUTTERWORTH).Band(butterBandLow, butterBandHigh);
	REQUIRE(desc.lower == butterBandLow);
	REQUIRE(desc.upper == butterBandHigh);
}

TEST_CASE("Band stop butterworth", "[IIR Descs]") {
	const auto desc = Bandstop(BUTTERWORTH).Band(butterBandLow, butterBandHigh);
	REQUIRE(desc.lower == butterBandLow);
	REQUIRE(desc.upper == butterBandHigh);
}


//------------------------------------------------------------------------------
// Chebyshev 1
//------------------------------------------------------------------------------

constexpr float cheby1Cutoff = 0.31f;
constexpr float cheby1Ripple = 0.05f;
constexpr float cheby1BandLow = 0.41f;
constexpr float cheby1BandHigh = 0.61f;

TEST_CASE("Low pass chebyshev 1", "[IIR Descs]") {
	const auto desc = Lowpass(CHEBYSHEV1).Cutoff(cheby1Cutoff).PassbandRipple(cheby1Ripple);
	REQUIRE(desc.cutoff == cheby1Cutoff);
	REQUIRE(desc.passbandRipple == cheby1Ripple);
}

TEST_CASE("High pass chebyshev 1", "[IIR Descs]") {
	const auto desc = Highpass(CHEBYSHEV1).Cutoff(cheby1Cutoff).PassbandRipple(cheby1Ripple);
	REQUIRE(desc.cutoff == cheby1Cutoff);
	REQUIRE(desc.passbandRipple == cheby1Ripple);
}

TEST_CASE("Band pass chebyshev 1", "[IIR Descs]") {
	const auto desc = Bandpass(CHEBYSHEV1).Band(cheby1BandLow, cheby1BandHigh).PassbandRipple(cheby1Ripple);
	REQUIRE(desc.lower == cheby1BandLow);
	REQUIRE(desc.upper == cheby1BandHigh);
	REQUIRE(desc.passbandRipple == cheby1Ripple);
}

TEST_CASE("Band stop chebyshev 1", "[IIR Descs]") {
	const auto desc = Bandstop(CHEBYSHEV1).Band(cheby1BandLow, cheby1BandHigh).PassbandRipple(cheby1Ripple);
	REQUIRE(desc.lower == cheby1BandLow);
	REQUIRE(desc.upper == cheby1BandHigh);
	REQUIRE(desc.passbandRipple == cheby1Ripple);
}


//------------------------------------------------------------------------------
// Chebyshev 2
//------------------------------------------------------------------------------

constexpr float cheby2Cutoff = 0.32f;
constexpr float cheby2Ripple = 0.04f;
constexpr float cheby2BandLow = 0.42f;
constexpr float cheby2BandHigh = 0.62f;

TEST_CASE("Low pass chebyshev 2", "[IIR Descs]") {
	const auto desc = Lowpass(CHEBYSHEV2).Cutoff(cheby2Cutoff).StopbandRipple(cheby2Ripple);
	REQUIRE(desc.cutoff == cheby2Cutoff);
	REQUIRE(desc.stopbandRipple == cheby2Ripple);
}

TEST_CASE("High pass chebyshev 2", "[IIR Descs]") {
	const auto desc = Highpass(CHEBYSHEV2).Cutoff(cheby2Cutoff).StopbandRipple(cheby2Ripple);
	REQUIRE(desc.cutoff == cheby2Cutoff);
	REQUIRE(desc.stopbandRipple == cheby2Ripple);
}

TEST_CASE("Band pass chebyshev 2", "[IIR Descs]") {
	const auto desc = Bandpass(CHEBYSHEV2).Band(cheby2BandLow, cheby2BandHigh).StopbandRipple(cheby2Ripple);
	REQUIRE(desc.lower == cheby2BandLow);
	REQUIRE(desc.upper == cheby2BandHigh);
	REQUIRE(desc.stopbandRipple == cheby2Ripple);
}

TEST_CASE("Band stop chebyshev 2", "[IIR Descs]") {
	const auto desc = Bandstop(CHEBYSHEV2).Band(cheby2BandLow, cheby2BandHigh).StopbandRipple(cheby2Ripple);
	REQUIRE(desc.lower == cheby2BandLow);
	REQUIRE(desc.upper == cheby2BandHigh);
	REQUIRE(desc.stopbandRipple == cheby2Ripple);
}


//------------------------------------------------------------------------------
// Elliptic
//------------------------------------------------------------------------------

constexpr float ellipticCutoff = 0.33f;
constexpr float ellipticPassRipple = 0.03f;
constexpr float ellipticStopRipple = 0.07f;
constexpr float ellipticBandLow = 0.43f;
constexpr float ellipticBandHigh = 0.63f;

TEST_CASE("Low pass elliptic", "[IIR Descs]") {
	const auto desc = Lowpass(ELLIPTIC).Cutoff(ellipticCutoff).PassbandRipple(ellipticPassRipple).StopbandRipple(ellipticStopRipple);
	REQUIRE(desc.cutoff == ellipticCutoff);
	REQUIRE(desc.passbandRipple == ellipticPassRipple);
	REQUIRE(desc.stopbandRipple == ellipticStopRipple);
}

TEST_CASE("High pass elliptic", "[IIR Descs]") {
	const auto desc = Highpass(ELLIPTIC).Cutoff(ellipticCutoff).PassbandRipple(ellipticPassRipple).StopbandRipple(ellipticStopRipple);
	REQUIRE(desc.cutoff == ellipticCutoff);
	REQUIRE(desc.passbandRipple == ellipticPassRipple);
	REQUIRE(desc.stopbandRipple == ellipticStopRipple);
}

TEST_CASE("Band pass elliptic", "[IIR Descs]") {
	const auto desc = Bandpass(ELLIPTIC).Band(ellipticBandLow, ellipticBandHigh).PassbandRipple(ellipticPassRipple).StopbandRipple(ellipticStopRipple);
	REQUIRE(desc.lower == ellipticBandLow);
	REQUIRE(desc.upper == ellipticBandHigh);
	REQUIRE(desc.passbandRipple == ellipticPassRipple);
	REQUIRE(desc.stopbandRipple == ellipticStopRipple);
}

TEST_CASE("Band stop elliptic", "[IIR Descs]") {
	const auto desc = Bandstop(ELLIPTIC).Band(ellipticBandLow, ellipticBandHigh).PassbandRipple(ellipticPassRipple).StopbandRipple(ellipticStopRipple);
	REQUIRE(desc.lower == ellipticBandLow);
	REQUIRE(desc.upper == ellipticBandHigh);
	REQUIRE(desc.passbandRipple == ellipticPassRipple);
	REQUIRE(desc.stopbandRipple == ellipticStopRipple);
}
