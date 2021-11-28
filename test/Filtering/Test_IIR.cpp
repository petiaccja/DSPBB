#include "../TestUtils.hpp"

#include <dspbb/Filtering/FilterParameters.hpp>
#include <dspbb/Filtering/IIR.hpp>
#include <dspbb/Math/Statistics.hpp>

#include <catch2/catch.hpp>

using namespace dspbb;


//------------------------------------------------------------------------------
// Filter application helpers
//------------------------------------------------------------------------------

TEST_CASE("Filter direct form I", "[IIR]") {
	constexpr int order = 7;
	const auto filter = TransferFunction(IirFilter<float>(order, Lowpass(BUTTERWORTH).Cutoff(0.3f)));
	DirectFormI<float> state{ order };
	const Signal<float, TIME_DOMAIN> signal(64, 1.0f);
	const auto filtered = Filter(signal, filter, state);
	REQUIRE(filtered.Size() == signal.Size());
}

TEST_CASE("Filter direct form II", "[IIR]") {
	constexpr int order = 7;
	const auto filter = TransferFunction(IirFilter<float>(order, Lowpass(BUTTERWORTH).Cutoff(0.3f)));
	DirectFormII<float> state{ order };
	const Signal<float, TIME_DOMAIN> signal(64, 1.0f);
	const auto filtered = Filter(signal, filter, state);
	REQUIRE(filtered.Size() == signal.Size());
}

TEST_CASE("Filter cascaded form", "[IIR]") {
	constexpr int order = 7;
	const auto filter = CascadedBiquad(IirFilter<float>(order, Lowpass(BUTTERWORTH).Cutoff(0.3f)));
	CascadedForm<float> state{ order };
	const Signal<float, TIME_DOMAIN> signal(64, 1.0f);
	const auto filtered = Filter(signal, filter, state);
	REQUIRE(filtered.Size() == signal.Size());
}

TEST_CASE("Filter overrun", "[IIR]") {
	constexpr int order = 7;
	const auto filter = TransferFunction(IirFilter<float>(order, Lowpass(BUTTERWORTH).Cutoff(0.3f)));
	DirectFormI<float> state{ order };
	const Signal<float, TIME_DOMAIN> signal(64, 1.0f);
	Signal<float, TIME_DOMAIN> out(1000, 1.0f);
	REQUIRE_THROWS(Filter(out, signal, filter, state));
}

TEST_CASE("Filter continuity", "[IIR]") {
	constexpr int order = 7;
	const auto filter = TransferFunction(IirFilter<float>(order, Lowpass(BUTTERWORTH).Cutoff(0.3f)));
	DirectFormI<float> state{ order };

	constexpr int length = 80;
	static_assert(length % 2 == 0);

	const auto signal = RandomSignal<double, TIME_DOMAIN>(length);
	const auto expected = Filter(signal, filter, state);
	state.Reset();

	TimeSignal<double> result(length);
	Filter(AsView(result).SubSignal(0, length / 2), AsView(signal).SubSignal(0, length / 2), filter, state);
	Filter(AsView(result).SubSignal(length / 2, length / 2), AsView(signal).SubSignal(length / 2, length / 2), filter, state);

	REQUIRE(Max(Abs(result - expected)) < 1e-4f);
}


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


//------------------------------------------------------------------------------
// Chebyshev 1 method
//------------------------------------------------------------------------------

constexpr float cheby1Cutoff = 0.75f;
constexpr float cheby1Lower = 0.35f;
constexpr float cheby1Upper = 0.65f;
constexpr float cheby1Ripple = 0.05f;
constexpr float cheby1RippleTolerance = 5e-4f;

TEST_CASE("Chebyshev 1 lowpass", "[IIR]") {
	constexpr int order = 7;
	const auto filter = TransferFunction(IirFilter<float>(order, Lowpass(CHEBYSHEV1).Cutoff(cheby1Cutoff).PassbandRipple(cheby1Ripple)));
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = ParametrizeLowpassFilter(amplitude);
	REQUIRE(params.passbandEdge < cheby1Cutoff);
	REQUIRE(params.stopbandEdge > cheby1Cutoff);
	REQUIRE(params.stopbandAtten < cheby1RippleTolerance);
	REQUIRE(params.passbandRipple == Approx(cheby1Ripple).margin(0.005));
}

TEST_CASE("Chebyshev 1 highpass", "[IIR]") {
	constexpr int order = 7;
	const auto filter = TransferFunction(IirFilter<float>(order, Highpass(CHEBYSHEV1).Cutoff(0.75f).PassbandRipple(cheby1Ripple)));
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = ParametrizeHighpassFilter(amplitude);
	REQUIRE(params.passbandEdge > cheby1Cutoff);
	REQUIRE(params.stopbandEdge < cheby1Cutoff);
	REQUIRE(params.stopbandAtten < cheby1RippleTolerance);
	REQUIRE(params.passbandRipple == Approx(cheby1Ripple).margin(0.005));
}

TEST_CASE("Chebyshev 1 bandpass", "[IIR]") {
	constexpr int order = 8;
	const auto filter = TransferFunction(IirFilter<float>(order, Bandpass(CHEBYSHEV1).Band(cheby1Lower, cheby1Upper).PassbandRipple(cheby1Ripple)));
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = ParametrizeBandpassFilter(amplitude);
	REQUIRE(params.lowerStopbandEdge < cheby1Lower);
	REQUIRE(params.passbandLowerEdge > cheby1Lower);
	REQUIRE(params.passbandUpperEdge < cheby1Upper);
	REQUIRE(params.upperStopbandEdge > cheby1Upper);
	REQUIRE(params.lowerStopbandAtten < cheby1RippleTolerance);
	REQUIRE(params.passbandRipple == Approx(cheby1Ripple).margin(0.005));
	REQUIRE(params.upperStopbandAtten < cheby1RippleTolerance);
}

TEST_CASE("Chebyshev 1 bandstop", "[IIR]") {
	constexpr int order = 8;
	const auto filter = TransferFunction(IirFilter<float>(order, Bandstop(CHEBYSHEV1).Band(cheby1Lower, cheby1Upper).PassbandRipple(cheby1Ripple)));
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = ParametrizeBandstopFilter(amplitude);
	REQUIRE(params.lowerPassbandEdge < cheby1Lower);
	REQUIRE(params.stopbandLowerEdge > cheby1Lower);
	REQUIRE(params.stopbandUpperEdge < cheby1Upper);
	REQUIRE(params.upperPassbandEdge > cheby1Upper);
	REQUIRE(params.lowerPassbandRipple == Approx(cheby1Ripple).margin(0.005));
	REQUIRE(params.stopbandAtten < cheby1RippleTolerance);
	REQUIRE(params.upperPassbandRipple == Approx(cheby1Ripple).margin(0.005));
}

TEST_CASE("Chebyshev 1 bandpass odd order", "[IIR]") {
	constexpr int order = 7;
	REQUIRE_THROWS(IirFilter<float>(order, Bandpass(CHEBYSHEV1).Band(cheby1Lower, cheby1Upper)));
}

TEST_CASE("Chebyshev 1 bandstop odd order", "[IIR]") {
	constexpr int order = 7;
	REQUIRE_THROWS(IirFilter<float>(order, Bandstop(CHEBYSHEV1).Band(cheby1Lower, cheby1Upper)));
}


//------------------------------------------------------------------------------
// Chebyshev 2 method
//------------------------------------------------------------------------------

constexpr float cheby2Cutoff = 0.75f;
constexpr float cheby2Lower = 0.35f;
constexpr float cheby2Upper = 0.65f;
constexpr float cheby2Ripple = 0.05f;
// TODO: make this tighter. Response seems good, maybe it's the parametrization.
constexpr float cheby2RippleTolerance = 2e-3f;

TEST_CASE("Chebyshev 2 lowpass", "[IIR]") {
	constexpr int order = 7;
	const auto filter = TransferFunction(IirFilter<float>(order, Lowpass(CHEBYSHEV2).Cutoff(cheby2Cutoff).StopbandRipple(cheby2Ripple)));
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = ParametrizeLowpassFilter(amplitude);
	REQUIRE(params.passbandEdge < cheby2Cutoff);
	REQUIRE(params.stopbandEdge > cheby2Cutoff);
	REQUIRE(params.stopbandAtten == Approx(cheby2Ripple).margin(0.005));
	REQUIRE(params.passbandRipple < cheby2RippleTolerance);
}

TEST_CASE("Chebyshev 2 highpass", "[IIR]") {
	constexpr int order = 7;
	const auto filter = TransferFunction(IirFilter<float>(order, Highpass(CHEBYSHEV2).Cutoff(0.75f).StopbandRipple(cheby2Ripple)));
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = ParametrizeHighpassFilter(amplitude);
	REQUIRE(params.passbandEdge > cheby2Cutoff);
	REQUIRE(params.stopbandEdge < cheby2Cutoff);
	REQUIRE(params.stopbandAtten == Approx(cheby2Ripple).margin(0.005));
	REQUIRE(params.passbandRipple < cheby2RippleTolerance);
}

TEST_CASE("Chebyshev 2 bandpass", "[IIR]") {
	constexpr int order = 8;
	const auto filter = TransferFunction(IirFilter<float>(order, Bandpass(CHEBYSHEV2).Band(cheby2Lower, cheby2Upper).StopbandRipple(cheby2Ripple)));
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = ParametrizeBandpassFilter(amplitude);
	REQUIRE(params.lowerStopbandEdge < cheby2Lower);
	REQUIRE(params.passbandLowerEdge > cheby2Lower);
	REQUIRE(params.passbandUpperEdge < cheby2Upper);
	REQUIRE(params.upperStopbandEdge > cheby2Upper);
	REQUIRE(params.lowerStopbandAtten == Approx(cheby2Ripple).margin(0.005));
	REQUIRE(params.passbandRipple < cheby2RippleTolerance);
	REQUIRE(params.upperStopbandAtten == Approx(cheby2Ripple).margin(0.005));
}

TEST_CASE("Chebyshev 2 bandstop", "[IIR]") {
	constexpr int order = 8;
	const auto filter = TransferFunction(IirFilter<float>(order, Bandstop(CHEBYSHEV2).Band(cheby2Lower, cheby2Upper).StopbandRipple(cheby2Ripple)));
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = ParametrizeBandstopFilter(amplitude);
	REQUIRE(params.lowerPassbandEdge < cheby2Lower);
	REQUIRE(params.stopbandLowerEdge > cheby2Lower);
	REQUIRE(params.stopbandUpperEdge < cheby2Upper);
	REQUIRE(params.upperPassbandEdge > cheby2Upper);
	REQUIRE(params.lowerPassbandRipple < cheby2RippleTolerance);
	REQUIRE(params.stopbandAtten == Approx(cheby2Ripple).margin(0.005));
	REQUIRE(params.upperPassbandRipple < cheby2RippleTolerance);
}

TEST_CASE("Chebyshev 2 bandpass odd order", "[IIR]") {
	constexpr int order = 7;
	REQUIRE_THROWS(IirFilter<float>(order, Bandpass(CHEBYSHEV2).Band(cheby2Lower, cheby2Upper)));
}

TEST_CASE("Chebyshev 2 bandstop odd order", "[IIR]") {
	constexpr int order = 7;
	REQUIRE_THROWS(IirFilter<float>(order, Bandstop(CHEBYSHEV2).Band(cheby2Lower, cheby2Upper)));
}


//------------------------------------------------------------------------------
// Elliptic method
//------------------------------------------------------------------------------

constexpr float ellipticCutoff = 0.75f;
constexpr float ellipticLower = 0.35f;
constexpr float ellipticUpper = 0.65f;
constexpr float ellipticPassRipple = 0.05f;
constexpr float ellipticStopRipple = 0.05f;

TEST_CASE("Elliptic lowpass", "[IIR]") {
	constexpr int order = 5;
	const auto filter = TransferFunction(IirFilter<float>(order, Lowpass(ELLIPTIC).Cutoff(ellipticCutoff).PassbandRipple(ellipticPassRipple).StopbandRipple(ellipticStopRipple)));
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = ParametrizeLowpassFilter(amplitude);
	REQUIRE(params.passbandEdge < ellipticCutoff);
	REQUIRE(params.stopbandEdge > ellipticCutoff);
	REQUIRE(params.stopbandAtten == Approx(ellipticStopRipple).margin(0.005));
	REQUIRE(params.passbandRipple == Approx(ellipticPassRipple).margin(0.005));
}

TEST_CASE("Elliptic highpass", "[IIR]") {
	constexpr int order = 5;
	const auto filter = TransferFunction(IirFilter<float>(order, Highpass(ELLIPTIC).Cutoff(0.75f).PassbandRipple(ellipticPassRipple).StopbandRipple(ellipticStopRipple)));
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = ParametrizeHighpassFilter(amplitude);
	REQUIRE(params.passbandEdge > ellipticCutoff);
	REQUIRE(params.stopbandEdge < ellipticCutoff);
	REQUIRE(params.stopbandAtten == Approx(ellipticStopRipple).margin(0.005));
	REQUIRE(params.passbandRipple == Approx(ellipticPassRipple).margin(0.005));
}

TEST_CASE("Elliptic bandpass", "[IIR]") {
	constexpr int order = 6;
	const auto filter = TransferFunction(IirFilter<float>(order, Bandpass(ELLIPTIC).Band(ellipticLower, ellipticUpper).PassbandRipple(ellipticPassRipple).StopbandRipple(ellipticStopRipple)));
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = ParametrizeBandpassFilter(amplitude);
	REQUIRE(params.lowerStopbandEdge < ellipticLower);
	REQUIRE(params.passbandLowerEdge > ellipticLower);
	REQUIRE(params.passbandUpperEdge < ellipticUpper);
	REQUIRE(params.upperStopbandEdge > ellipticUpper);
	REQUIRE(params.lowerStopbandAtten == Approx(ellipticStopRipple).margin(0.005));
	REQUIRE(params.passbandRipple == Approx(ellipticPassRipple).margin(0.005));
	REQUIRE(params.upperStopbandAtten == Approx(ellipticStopRipple).margin(0.005));
}

TEST_CASE("Elliptic bandstop", "[IIR]") {
	constexpr int order = 6;
	const auto filter = TransferFunction(IirFilter<float>(order, Bandstop(ELLIPTIC).Band(ellipticLower, ellipticUpper).PassbandRipple(ellipticPassRipple).StopbandRipple(ellipticStopRipple)));
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = ParametrizeBandstopFilter(amplitude);
	REQUIRE(params.lowerPassbandEdge < ellipticLower);
	REQUIRE(params.stopbandLowerEdge > ellipticLower);
	REQUIRE(params.stopbandUpperEdge < ellipticUpper);
	REQUIRE(params.upperPassbandEdge > ellipticUpper);
	REQUIRE(params.lowerPassbandRipple == Approx(ellipticPassRipple).margin(0.005));
	REQUIRE(params.stopbandAtten == Approx(ellipticStopRipple).margin(0.005));
	REQUIRE(params.upperPassbandRipple == Approx(ellipticPassRipple).margin(0.005));
}

TEST_CASE("Elliptic bandpass odd order", "[IIR]") {
	constexpr int order = 7;
	REQUIRE_THROWS(IirFilter<float>(order, Bandpass(ELLIPTIC).Band(ellipticLower, ellipticUpper)));
}

TEST_CASE("Elliptic bandstop odd order", "[IIR]") {
	constexpr int order = 7;
	REQUIRE_THROWS(IirFilter<float>(order, Bandstop(ELLIPTIC).Band(ellipticLower, ellipticUpper)));
}