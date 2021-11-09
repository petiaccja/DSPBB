#include <catch2/catch.hpp>
#include <dspbb/Filtering/FilterParameters.hpp>
#include <dspbb/Filtering/IIR.hpp>
#include <dspbb/Math/Statistics.hpp>

using namespace dspbb;

//------------------------------------------------------------------------------
// Butterworth method
//------------------------------------------------------------------------------

constexpr float butterCutoff = 0.75f;
constexpr float butterLower = 0.35f;
constexpr float butterUpper = 0.65f;
constexpr float butterRippleTolerance = 1e-4f;

TEST_CASE("Butterworth lowpass", "[IIR]") {
	constexpr int order = 7;
	const auto filter = TransferFunction(IirFilter<float>(order, Lowpass(BUTTERWORTH).Cutoff(butterCutoff)));
	const auto [amplitude, phase] = FrequencyResponse(filter);
	const auto params = ParametrizeLowpassFilter(amplitude);
	REQUIRE(params.passbandEdge < butterCutoff);
	REQUIRE(params.stopbandEdge > butterCutoff);
	REQUIRE(params.stopbandAtten < butterRippleTolerance);
	REQUIRE(params.passbandRipple < butterRippleTolerance);
}

TEST_CASE("Butterworth highpass", "[IIR]") {
	constexpr int order = 7;
	const auto filter = TransferFunction(IirFilter<float>(order, Highpass(BUTTERWORTH).Cutoff(0.75f)));
	const auto [amplitude, phase] = FrequencyResponse(filter);
	const auto params = ParametrizeHighpassFilter(amplitude);
	REQUIRE(params.passbandEdge > butterCutoff);
	REQUIRE(params.stopbandEdge < butterCutoff);
	REQUIRE(params.stopbandAtten < butterRippleTolerance);
	REQUIRE(params.passbandRipple < butterRippleTolerance);
}

TEST_CASE("Butterworth bandpass", "[IIR]") {
	constexpr int order = 8;
	const auto filter = TransferFunction(IirFilter<float>(order, Bandpass(BUTTERWORTH).Band(butterLower, butterUpper)));
	const auto [amplitude, phase] = FrequencyResponse(filter);
	const auto params = ParametrizeBandpassFilter(amplitude);
	REQUIRE(params.lowerStopbandEdge < butterLower);
	REQUIRE(params.passbandLowerEdge > butterLower);
	REQUIRE(params.passbandUpperEdge < butterUpper);
	REQUIRE(params.upperStopbandEdge > butterUpper);
	REQUIRE(params.lowerStopbandAtten < butterRippleTolerance);
	REQUIRE(params.passbandRipple < butterRippleTolerance);
	REQUIRE(params.upperStopbandAtten < butterRippleTolerance);
}

TEST_CASE("Butterworth bandstop", "[IIR]") {
	constexpr int order = 8;
	const auto filter = TransferFunction(IirFilter<float>(order, Bandstop(BUTTERWORTH).Band(butterLower, butterUpper)));
	const auto [amplitude, phase] = FrequencyResponse(filter);
	const auto params = ParametrizeBandstopFilter(amplitude);
	REQUIRE(params.lowerPassbandEdge < butterLower);
	REQUIRE(params.stopbandLowerEdge > butterLower);
	REQUIRE(params.stopbandUpperEdge < butterUpper);
	REQUIRE(params.upperPassbandEdge > butterUpper);
	REQUIRE(params.lowerPassbandRipple < butterRippleTolerance);
	REQUIRE(params.stopbandAtten < butterRippleTolerance);
	REQUIRE(params.upperPassbandRipple < butterRippleTolerance);
}

TEST_CASE("Butterworth bandpass odd order", "[IIR]") {
	constexpr int order = 7;
	REQUIRE_THROWS(IirFilter<float>(order, Bandpass(BUTTERWORTH).Band(butterLower, butterUpper)));
}

TEST_CASE("Butterworth bandstop odd order", "[IIR]") {
	constexpr int order = 7;
	REQUIRE_THROWS(IirFilter<float>(order, Bandstop(BUTTERWORTH).Band(butterLower, butterUpper)));
}