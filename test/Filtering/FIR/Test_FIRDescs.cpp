#include <dspbb/Filtering/FIR.hpp>

#include <catch2/catch.hpp>


using namespace dspbb;

//------------------------------------------------------------------------------
// Windowed
//------------------------------------------------------------------------------

const Signal<float, TIME_DOMAIN> winWindow = { 1, 2, 3 };
constexpr float winCutoff = 0.3f;
constexpr float winBandLow = 0.4f;
constexpr float winBandHigh = 0.6f;

TEST_CASE("Low pass windowed view", "[FIR Descs]") {
	const auto desc = Lowpass(WINDOWED).Cutoff(winCutoff).Window(winWindow);
	REQUIRE(desc.cutoff == winCutoff);
	REQUIRE(desc.window.Size() == 3);
	REQUIRE(desc.window[0] == 1);
}

TEST_CASE("Low pass windowed function", "[FIR Descs]") {
	const auto desc = Lowpass(WINDOWED).Cutoff(winCutoff).Window(windows::hamming);
	REQUIRE(desc.cutoff == winCutoff);
	Signal<float, TIME_DOMAIN> window(3);
	desc.window(window);
	REQUIRE(window.Size() == 3);
	REQUIRE(Max(window) == Approx(1));
}

TEST_CASE("High pass windowed view", "[FIR Descs]") {
	const auto desc = Highpass(WINDOWED).Cutoff(winCutoff).Window(winWindow);
	REQUIRE(desc.cutoff == winCutoff);
	REQUIRE(desc.window.Size() == 3);
	REQUIRE(desc.window[0] == 1);
}

TEST_CASE("High pass windowed function", "[FIR Descs]") {
	const auto desc = Highpass(WINDOWED).Cutoff(winCutoff).Window(windows::hamming);
	REQUIRE(desc.cutoff == winCutoff);
	Signal<float, TIME_DOMAIN> window(3);
	desc.window(window);
	REQUIRE(window.Size() == 3);
	REQUIRE(Max(window) == Approx(1));
}

TEST_CASE("Band pass windowed view", "[FIR Descs]") {
	const auto desc = Bandpass(WINDOWED).Band(winBandLow, winBandHigh).Window(winWindow);
	REQUIRE(desc.lower == winBandLow);
	REQUIRE(desc.upper == winBandHigh);
	REQUIRE(desc.window.Size() == 3);
	REQUIRE(desc.window[0] == 1);
}

TEST_CASE("Band pass windowed function", "[FIR Descs]") {
	const auto desc = Bandpass(WINDOWED).Band(winBandLow, winBandHigh).Window(windows::hamming);
	REQUIRE(desc.lower == winBandLow);
	REQUIRE(desc.upper == winBandHigh);
	Signal<float, TIME_DOMAIN> window(3);
	desc.window(window);
	REQUIRE(window.Size() == 3);
	REQUIRE(Max(window) == Approx(1));
}

TEST_CASE("Band stop windowed view", "[FIR Descs]") {
	const auto desc = Bandstop(WINDOWED).Band(winBandLow, winBandHigh).Window(winWindow);
	REQUIRE(desc.lower == winBandLow);
	REQUIRE(desc.upper == winBandHigh);
	REQUIRE(desc.window.Size() == 3);
	REQUIRE(desc.window[0] == 1);
}

TEST_CASE("Band stop windowed function", "[FIR Descs]") {
	const auto desc = Bandstop(WINDOWED).Band(winBandLow, winBandHigh).Window(windows::hamming);
	REQUIRE(desc.lower == winBandLow);
	REQUIRE(desc.upper == winBandHigh);
	Signal<float, TIME_DOMAIN> window(3);
	desc.window(window);
	REQUIRE(window.Size() == 3);
	REQUIRE(Max(window) == Approx(1));
}

TEST_CASE("Arbitrary windowed view", "[FIR Descs]") {
	const auto desc = Arbitrary(WINDOWED).Response([](float) { return 1.f; }).Window(winWindow);
	REQUIRE(desc.responseFunc(0.3f) == Approx(1.0f));
	REQUIRE(desc.window.Size() == 3);
	REQUIRE(desc.window[0] == 1);
}

TEST_CASE("Arbitrary windowed function", "[FIR Descs]") {
	const auto desc = Arbitrary(WINDOWED).Response([](float) { return 1.f; }).Window(windows::blackman);
	REQUIRE(desc.responseFunc(0.3f) == Approx(1.0f));
	Signal<float, TIME_DOMAIN> window(3);
	desc.window(window);
	REQUIRE(window.Size() == 3);
	REQUIRE(Max(window) == Approx(1));
}

TEST_CASE("Hilbert windowed view", "[FIR Descs]") {
	const auto desc = Hilbert(WINDOWED).Window(winWindow);
	REQUIRE(desc.window.Size() == 3);
	REQUIRE(desc.window[0] == 1);
}

TEST_CASE("Hilbert windowed function", "[FIR Descs]") {
	const auto desc = Hilbert(WINDOWED).Window(windows::blackman);
	Signal<float, TIME_DOMAIN> window(3);
	desc.window(window);
	REQUIRE(window.Size() == 3);
	REQUIRE(Max(window) == Approx(1));
}

//------------------------------------------------------------------------------
// Least Squares
//------------------------------------------------------------------------------

constexpr float lsBegin1 = 0.28f;
constexpr float lsEnd1 = 0.32f;
constexpr float lsBegin2 = 0.68f;
constexpr float lsEnd2 = 0.72f;
constexpr float lsWeightLow = 2.0f;
constexpr float lsWeightTr1 = 0.1f;
constexpr float lsWeightMid = 0.1f;
constexpr float lsWeightTr2 = 0.1f;
constexpr float lsWeightHigh = 0.1f;

TEST_CASE("Low pass least squares", "[FIR Descs]") {
	const auto desc = Lowpass(LEAST_SQUARES).Cutoff(lsBegin1, lsEnd1).Weight(lsWeightLow, lsWeightTr1, lsWeightHigh);
	REQUIRE(desc.cutoffBegin == lsBegin1);
	REQUIRE(desc.cutoffEnd == lsEnd1);
	REQUIRE(desc.weightLow == lsWeightLow);
	REQUIRE(desc.weightTransition == lsWeightTr1);
	REQUIRE(desc.weightHigh == lsWeightHigh);
}

TEST_CASE("High pass least squares", "[FIR Descs]") {
	const auto desc = Highpass(LEAST_SQUARES).Cutoff(lsBegin1, lsEnd1).Weight(lsWeightLow, lsWeightTr1, lsWeightHigh);
	REQUIRE(desc.cutoffBegin == lsBegin1);
	REQUIRE(desc.cutoffEnd == lsEnd1);
	REQUIRE(desc.weightLow == lsWeightLow);
	REQUIRE(desc.weightTransition == lsWeightTr1);
	REQUIRE(desc.weightHigh == lsWeightHigh);
}

TEST_CASE("Band pass least squares", "[FIR Descs]") {
	const auto desc = Bandpass(LEAST_SQUARES).Band(lsBegin1, lsEnd1, lsBegin2, lsEnd2).Weight(lsWeightLow, lsWeightTr1, lsWeightMid, lsWeightTr2, lsWeightHigh);
	REQUIRE(desc.lowerBegin == lsBegin1);
	REQUIRE(desc.lowerEnd == lsEnd1);
	REQUIRE(desc.upperBegin == lsBegin2);
	REQUIRE(desc.upperEnd == lsEnd2);
	REQUIRE(desc.weightLow == lsWeightLow);
	REQUIRE(desc.weightTransition1 == lsWeightTr1);
	REQUIRE(desc.weightMid == lsWeightMid);
	REQUIRE(desc.weightTransition2 == lsWeightTr2);
	REQUIRE(desc.weightHigh == lsWeightHigh);
}

TEST_CASE("Band stop least squares", "[FIR Descs]") {
	const auto desc = Bandstop(LEAST_SQUARES).Band(lsBegin1, lsEnd1, lsBegin2, lsEnd2).Weight(lsWeightLow, lsWeightTr1, lsWeightMid, lsWeightTr2, lsWeightHigh);
	REQUIRE(desc.lowerBegin == lsBegin1);
	REQUIRE(desc.lowerEnd == lsEnd1);
	REQUIRE(desc.upperBegin == lsBegin2);
	REQUIRE(desc.upperEnd == lsEnd2);
	REQUIRE(desc.weightLow == lsWeightLow);
	REQUIRE(desc.weightTransition1 == lsWeightTr1);
	REQUIRE(desc.weightMid == lsWeightMid);
	REQUIRE(desc.weightTransition2 == lsWeightTr2);
	REQUIRE(desc.weightHigh == lsWeightHigh);
}

TEST_CASE("Arbitrary least squares", "[FIR Descs]") {
	const auto desc = Arbitrary(LEAST_SQUARES).Response([](float) { return 1.f; }).Weight([](float f) { return f < 0.5f ? 1.0f : 0.5f; });
	REQUIRE(desc.responseFunc(0.3f) == Approx(1.0f));
	REQUIRE(desc.weightFunc(0.4f) == Approx(1.0f));
	REQUIRE(desc.weightFunc(0.6f) == Approx(0.5f));
}

TEST_CASE("Hilbert least squares", "[FIR Descs]") {
	const auto desc = Hilbert(LEAST_SQUARES).TransitionWidth(0.95f);
	REQUIRE(desc.transitionWidth == Approx(0.95f));
}