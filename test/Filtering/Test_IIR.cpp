#include "../TestUtils.hpp"

#include <dspbb/Filtering/IIR.hpp>
#include <dspbb/Filtering/MeasureFilter.hpp>
#include <dspbb/Math/Statistics.hpp>

#include <catch2/catch.hpp>

using namespace dspbb;


//------------------------------------------------------------------------------
// Filter application helpers
//------------------------------------------------------------------------------

TEST_CASE("Filter direct form I", "[IIR]") {
	constexpr int order = 7;
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Lowpass.Butterworth.Cutoff(0.3f)));
	DirectFormI<float> state{ order };
	const BasicSignal<float, TIME_DOMAIN> signal(64, 1.0f);
	const auto filtered = Filter(signal, filter, state);
	REQUIRE(filtered.size() == signal.size());
}

TEST_CASE("Filter direct form II", "[IIR]") {
	constexpr int order = 7;
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Lowpass.Butterworth.Cutoff(0.3f)));
	DirectFormII<float> state{ order };
	const BasicSignal<float, TIME_DOMAIN> signal(64, 1.0f);
	const auto filtered = Filter(signal, filter, state);
	REQUIRE(filtered.size() == signal.size());
}

TEST_CASE("Filter cascaded form", "[IIR]") {
	constexpr int order = 7;
	const auto filter = CascadedBiquad(DesignFilter<float>(order, Iir.Lowpass.Butterworth.Cutoff(0.3f)));
	CascadedForm<float> state{ order };
	const BasicSignal<float, TIME_DOMAIN> signal(64, 1.0f);
	const auto filtered = Filter(signal, filter, state);
	REQUIRE(filtered.size() == signal.size());
}

TEST_CASE("Filter continuity", "[IIR]") {
	constexpr int order = 7;
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Lowpass.Butterworth.Cutoff(0.3f)));
	DirectFormI<float> state{ order };

	constexpr int length = 80;
	static_assert(length % 2 == 0);

	const auto signal = RandomSignal<double, TIME_DOMAIN>(length);
	const auto expected = Filter(signal, filter, state);
	state.Reset();

	Signal<double> result(length);
	Filter(AsView(result).subsignal(0, length / 2), AsView(signal).subsignal(0, length / 2), filter, state);
	Filter(AsView(result).subsignal(length / 2, length / 2), AsView(signal).subsignal(length / 2, length / 2), filter, state);

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
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Lowpass.Butterworth.Cutoff(butterCutoff)));
	REQUIRE(filter.Order() == 7);
	const auto [amplitude, phase] = FrequencyResponse(filter);
	const auto params = MeasureLowpassFilter(amplitude);
	REQUIRE(params.passbandEdge < butterCutoff);
	REQUIRE(params.stopbandEdge > butterCutoff);
	REQUIRE(params.stopbandAtten < butterRippleTolerance);
	REQUIRE(params.passbandRipple < butterRippleTolerance);
}

TEST_CASE("Butterworth highpass", "[IIR]") {
	constexpr int order = 7;
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Highpass.Butterworth.Cutoff(0.75f)));
	REQUIRE(filter.Order() == 7);
	const auto [amplitude, phase] = FrequencyResponse(filter);
	const auto params = MeasureHighpassFilter(amplitude);
	REQUIRE(params.passbandEdge > butterCutoff);
	REQUIRE(params.stopbandEdge < butterCutoff);
	REQUIRE(params.stopbandAtten < butterRippleTolerance);
	REQUIRE(params.passbandRipple < butterRippleTolerance);
}

TEST_CASE("Butterworth bandpass", "[IIR]") {
	constexpr int order = 8;
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Bandpass.Butterworth.Band(butterLower, butterUpper)));
	REQUIRE(filter.Order() == 8);
	const auto [amplitude, phase] = FrequencyResponse(filter);
	const auto params = MeasureBandpassFilter(amplitude);
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
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Bandstop.Butterworth.Band(butterLower, butterUpper)));
	REQUIRE(filter.Order() == 8);
	const auto [amplitude, phase] = FrequencyResponse(filter);
	const auto params = MeasureBandstopFilter(amplitude);
	REQUIRE(params.lowerPassbandEdge < butterLower);
	REQUIRE(params.stopbandLowerEdge > butterLower);
	REQUIRE(params.stopbandUpperEdge < butterUpper);
	REQUIRE(params.upperPassbandEdge > butterUpper);
	REQUIRE(params.lowerPassbandRipple < butterRippleTolerance);
	REQUIRE(params.stopbandAtten < butterRippleTolerance);
	REQUIRE(params.upperPassbandRipple < butterRippleTolerance);
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
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Lowpass.Chebyshev1.Cutoff(cheby1Cutoff).PassbandRipple(cheby1Ripple)));
	REQUIRE(filter.Order() == 7);
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = MeasureLowpassFilter(amplitude);
	REQUIRE(params.passbandEdge < cheby1Cutoff);
	REQUIRE(params.stopbandEdge > cheby1Cutoff);
	REQUIRE(params.stopbandAtten < cheby1RippleTolerance);
	REQUIRE(params.passbandRipple == Approx(cheby1Ripple).margin(0.005));
}

TEST_CASE("Chebyshev 1 highpass", "[IIR]") {
	constexpr int order = 7;
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Highpass.Chebyshev1.Cutoff(0.75f).PassbandRipple(cheby1Ripple)));
	REQUIRE(filter.Order() == 7);
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = MeasureHighpassFilter(amplitude);
	REQUIRE(params.passbandEdge > cheby1Cutoff);
	REQUIRE(params.stopbandEdge < cheby1Cutoff);
	REQUIRE(params.stopbandAtten < cheby1RippleTolerance);
	REQUIRE(params.passbandRipple == Approx(cheby1Ripple).margin(0.005));
}

TEST_CASE("Chebyshev 1 bandpass", "[IIR]") {
	constexpr int order = 8;
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Bandpass.Chebyshev1.Band(cheby1Lower, cheby1Upper).PassbandRipple(cheby1Ripple)));
	REQUIRE(filter.Order() == 8);
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = MeasureBandpassFilter(amplitude);
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
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Bandstop.Chebyshev1.Band(cheby1Lower, cheby1Upper).PassbandRipple(cheby1Ripple)));
	REQUIRE(filter.Order() == 8);
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = MeasureBandstopFilter(amplitude);
	REQUIRE(params.lowerPassbandEdge < cheby1Lower);
	REQUIRE(params.stopbandLowerEdge > cheby1Lower);
	REQUIRE(params.stopbandUpperEdge < cheby1Upper);
	REQUIRE(params.upperPassbandEdge > cheby1Upper);
	REQUIRE(params.lowerPassbandRipple == Approx(cheby1Ripple).margin(0.005));
	REQUIRE(params.stopbandAtten < cheby1RippleTolerance);
	REQUIRE(params.upperPassbandRipple == Approx(cheby1Ripple).margin(0.005));
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
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Lowpass.Chebyshev2.Cutoff(cheby2Cutoff).StopbandRipple(cheby2Ripple)));
	REQUIRE(filter.Order() == 7);
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = MeasureLowpassFilter(amplitude);
	REQUIRE(params.passbandEdge < cheby2Cutoff);
	REQUIRE(params.stopbandEdge > cheby2Cutoff);
	REQUIRE(params.stopbandAtten == Approx(cheby2Ripple).margin(0.005));
	REQUIRE(params.passbandRipple < cheby2RippleTolerance);
}

TEST_CASE("Chebyshev 2 highpass", "[IIR]") {
	constexpr int order = 7;
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Highpass.Chebyshev2.Cutoff(0.75f).StopbandRipple(cheby2Ripple)));
	REQUIRE(filter.Order() == 7);
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = MeasureHighpassFilter(amplitude);
	REQUIRE(params.passbandEdge > cheby2Cutoff);
	REQUIRE(params.stopbandEdge < cheby2Cutoff);
	REQUIRE(params.stopbandAtten == Approx(cheby2Ripple).margin(0.005));
	REQUIRE(params.passbandRipple < cheby2RippleTolerance);
}

TEST_CASE("Chebyshev 2 bandpass", "[IIR]") {
	constexpr int order = 8;
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Bandpass.Chebyshev2.Band(cheby2Lower, cheby2Upper).StopbandRipple(cheby2Ripple)));
	REQUIRE(filter.Order() == 8);
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = MeasureBandpassFilter(amplitude);
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
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Bandstop.Chebyshev2.Band(cheby2Lower, cheby2Upper).StopbandRipple(cheby2Ripple)));
	REQUIRE(filter.Order() == 8);
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = MeasureBandstopFilter(amplitude);
	REQUIRE(params.lowerPassbandEdge < cheby2Lower);
	REQUIRE(params.stopbandLowerEdge > cheby2Lower);
	REQUIRE(params.stopbandUpperEdge < cheby2Upper);
	REQUIRE(params.upperPassbandEdge > cheby2Upper);
	REQUIRE(params.lowerPassbandRipple < cheby2RippleTolerance);
	REQUIRE(params.stopbandAtten == Approx(cheby2Ripple).margin(0.005));
	REQUIRE(params.upperPassbandRipple < cheby2RippleTolerance);
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
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Lowpass.Elliptic.Cutoff(ellipticCutoff).PassbandRipple(ellipticPassRipple).StopbandRipple(ellipticStopRipple)));
	REQUIRE(filter.Order() == 5);
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = MeasureLowpassFilter(amplitude);
	REQUIRE(params.passbandEdge < ellipticCutoff);
	REQUIRE(params.stopbandEdge > ellipticCutoff);
	REQUIRE(params.stopbandAtten == Approx(ellipticStopRipple).margin(0.005));
	REQUIRE(params.passbandRipple == Approx(ellipticPassRipple).margin(0.005));
}

TEST_CASE("Elliptic highpass", "[IIR]") {
	constexpr int order = 5;
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Highpass.Elliptic.Cutoff(0.75f).PassbandRipple(ellipticPassRipple).StopbandRipple(ellipticStopRipple)));
	REQUIRE(filter.Order() == 5);
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = MeasureHighpassFilter(amplitude);
	REQUIRE(params.passbandEdge > ellipticCutoff);
	REQUIRE(params.stopbandEdge < ellipticCutoff);
	REQUIRE(params.stopbandAtten == Approx(ellipticStopRipple).margin(0.005));
	REQUIRE(params.passbandRipple == Approx(ellipticPassRipple).margin(0.005));
}

TEST_CASE("Elliptic bandpass", "[IIR]") {
	constexpr int order = 6;
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Bandpass.Elliptic.Band(ellipticLower, ellipticUpper).PassbandRipple(ellipticPassRipple).StopbandRipple(ellipticStopRipple)));
	REQUIRE(filter.Order() == 6);
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = MeasureBandpassFilter(amplitude);
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
	const auto filter = TransferFunction(DesignFilter<float>(order, Iir.Bandstop.Elliptic.Band(ellipticLower, ellipticUpper).PassbandRipple(ellipticPassRipple).StopbandRipple(ellipticStopRipple)));
	REQUIRE(filter.Order() == 6);
	const auto [amplitude, phase] = FrequencyResponse(filter, 8192);
	const auto params = MeasureBandstopFilter(amplitude);
	REQUIRE(params.lowerPassbandEdge < ellipticLower);
	REQUIRE(params.stopbandLowerEdge > ellipticLower);
	REQUIRE(params.stopbandUpperEdge < ellipticUpper);
	REQUIRE(params.upperPassbandEdge > ellipticUpper);
	REQUIRE(params.lowerPassbandRipple == Approx(ellipticPassRipple).margin(0.005));
	REQUIRE(params.stopbandAtten == Approx(ellipticStopRipple).margin(0.005));
	REQUIRE(params.upperPassbandRipple == Approx(ellipticPassRipple).margin(0.005));
}