#include <catch2/catch.hpp>
#include <dspbb/Filtering/FIRDescs.hpp>


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
	REQUIRE(desc.low == winBandLow);
	REQUIRE(desc.high == winBandHigh);
	REQUIRE(desc.window.Size() == 3);
	REQUIRE(desc.window[0] == 1);
}

TEST_CASE("Band pass windowed function", "[FIR Descs]") {
	const auto desc = Bandpass(WINDOWED).Band(winBandLow, winBandHigh).Window(windows::hamming);
	REQUIRE(desc.low == winBandLow);
	REQUIRE(desc.high == winBandHigh);
	Signal<float, TIME_DOMAIN> window(3);
	desc.window(window);
	REQUIRE(window.Size() == 3);
	REQUIRE(Max(window) == Approx(1));
}

TEST_CASE("Band stop windowed view", "[FIR Descs]") {
	const auto desc = Bandstop(WINDOWED).Band(winBandLow, winBandHigh).Window(winWindow);
	REQUIRE(desc.low == winBandLow);
	REQUIRE(desc.high == winBandHigh);
	REQUIRE(desc.window.Size() == 3);
	REQUIRE(desc.window[0] == 1);
}

TEST_CASE("Band stop windowed function", "[FIR Descs]") {
	const auto desc = Bandstop(WINDOWED).Band(winBandLow, winBandHigh).Window(windows::hamming);
	REQUIRE(desc.low == winBandLow);
	REQUIRE(desc.high == winBandHigh);
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
constexpr bool lsSmooth1 = true;
constexpr bool lsSmooth2 = false;

TEST_CASE("Low pass least squares", "[FIR Descs]") {
	const auto desc = Lowpass(LEAST_SQUARES).Cutoff(lsBegin1, lsEnd1).Weight(lsWeightLow, lsWeightTr1, lsWeightHigh).Smooth(lsSmooth1);
	REQUIRE(desc.cutoffBegin == lsBegin1);
	REQUIRE(desc.cutoffEnd == lsEnd1);
	REQUIRE(desc.weightLow == lsWeightLow);
	REQUIRE(desc.weightTransition == lsWeightTr1);
	REQUIRE(desc.weightHigh == lsWeightHigh);
	REQUIRE(desc.smooth == lsSmooth1);
}

TEST_CASE("High pass least squares", "[FIR Descs]") {
	const auto desc = Highpass(LEAST_SQUARES).Cutoff(lsBegin1, lsEnd1).Weight(lsWeightLow, lsWeightTr1, lsWeightHigh).Smooth(lsSmooth1);
	REQUIRE(desc.cutoffBegin == lsBegin1);
	REQUIRE(desc.cutoffEnd == lsEnd1);
	REQUIRE(desc.weightLow == lsWeightLow);
	REQUIRE(desc.weightTransition == lsWeightTr1);
	REQUIRE(desc.weightHigh == lsWeightHigh);
	REQUIRE(desc.smooth == lsSmooth1);
}

TEST_CASE("Band pass least squares", "[FIR Descs]") {
	const auto desc = Bandpass(LEAST_SQUARES).Cutoff(lsBegin1, lsEnd1, lsBegin2, lsEnd2).Weight(lsWeightLow, lsWeightTr1, lsWeightMid, lsWeightTr2, lsWeightHigh).Smooth(lsSmooth1, lsSmooth2);
	REQUIRE(desc.cutoffBegin1 == lsBegin1);
	REQUIRE(desc.cutoffEnd1 == lsEnd1);
	REQUIRE(desc.cutoffBegin2 == lsBegin2);
	REQUIRE(desc.cutoffEnd2 == lsEnd2);
	REQUIRE(desc.weightLow == lsWeightLow);
	REQUIRE(desc.weightTransition1 == lsWeightTr1);
	REQUIRE(desc.weightMid == lsWeightMid);
	REQUIRE(desc.weightTransition2 == lsWeightTr2);
	REQUIRE(desc.weightHigh == lsWeightHigh);
	REQUIRE(desc.smooth1 == lsSmooth1);
	REQUIRE(desc.smooth2 == lsSmooth2);
}

TEST_CASE("Band stop least squares", "[FIR Descs]") {
	const auto desc = Bandstop(LEAST_SQUARES).Cutoff(lsBegin1, lsEnd1, lsBegin2, lsEnd2).Weight(lsWeightLow, lsWeightTr1, lsWeightMid, lsWeightTr2, lsWeightHigh).Smooth(lsSmooth1, lsSmooth2);
	REQUIRE(desc.cutoffBegin1 == lsBegin1);
	REQUIRE(desc.cutoffEnd1 == lsEnd1);
	REQUIRE(desc.cutoffBegin2 == lsBegin2);
	REQUIRE(desc.cutoffEnd2 == lsEnd2);
	REQUIRE(desc.weightLow == lsWeightLow);
	REQUIRE(desc.weightTransition1 == lsWeightTr1);
	REQUIRE(desc.weightMid == lsWeightMid);
	REQUIRE(desc.weightTransition2 == lsWeightTr2);
	REQUIRE(desc.weightHigh == lsWeightHigh);
	REQUIRE(desc.smooth1 == lsSmooth1);
	REQUIRE(desc.smooth2 == lsSmooth2);
}

TEST_CASE("Arbitrary least squares", "[FIR Descs]") {
	const auto desc = Arbitrary(LEAST_SQUARES).Response([](float) { return 1.f; }).Weight([](float f) { return f < 0.5f ? 1.0f : 0.5f; });
	REQUIRE(desc.responseFunc(0.3f) == Approx(1.0f));
	REQUIRE(desc.weightFunc(0.4f) == Approx(1.0f));
	REQUIRE(desc.weightFunc(0.6f) == Approx(0.5f));
}