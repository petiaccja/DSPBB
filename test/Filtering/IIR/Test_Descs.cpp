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
	REQUIRE(desc.low == butterBandLow);
	REQUIRE(desc.high == butterBandHigh);
}

TEST_CASE("Band stop butterworth", "[IIR Descs]") {
	const auto desc = Bandstop(BUTTERWORTH).Band(butterBandLow, butterBandHigh);
	REQUIRE(desc.low == butterBandLow);
	REQUIRE(desc.high == butterBandHigh);
}