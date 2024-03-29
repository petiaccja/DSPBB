#include <dspbb/Filtering/FIR.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>


using namespace dspbb;
using Catch::Approx;


//------------------------------------------------------------------------------
// Windowed
//------------------------------------------------------------------------------

const BasicSignal<float, TIME_DOMAIN> winWindow = { 1, 2, 3 };
constexpr float winCutoff = 0.3f;
constexpr float winBandLow = 0.4f;
constexpr float winBandHigh = 0.6f;

TEST_CASE("Low pass windowed view", "[FIR Descs]") {
	const auto desc = Fir.Lowpass.Windowed.Cutoff(winCutoff).Window(winWindow);
	REQUIRE(desc.cutoff == winCutoff);
	REQUIRE(desc.window.size() == 3);
	REQUIRE(desc.window[0] == 1);
}

TEST_CASE("Low pass windowed function", "[FIR Descs]") {
	const auto desc = Fir.Lowpass.Windowed.Cutoff(winCutoff).Window(windows::hamming);
	REQUIRE(desc.cutoff == winCutoff);
	BasicSignal<float, TIME_DOMAIN> window(3);
	desc.window(window);
	REQUIRE(window.size() == 3);
	REQUIRE(Max(window) == Approx(1));
}

TEST_CASE("High pass windowed view", "[FIR Descs]") {
	const auto desc = Fir.Highpass.Windowed.Cutoff(winCutoff).Window(winWindow);
	REQUIRE(desc.cutoff == winCutoff);
	REQUIRE(desc.window.size() == 3);
	REQUIRE(desc.window[0] == 1);
}

TEST_CASE("High pass windowed function", "[FIR Descs]") {
	const auto desc = Fir.Highpass.Windowed.Cutoff(winCutoff).Window(windows::hamming);
	REQUIRE(desc.cutoff == winCutoff);
	BasicSignal<float, TIME_DOMAIN> window(3);
	desc.window(window);
	REQUIRE(window.size() == 3);
	REQUIRE(Max(window) == Approx(1));
}

TEST_CASE("Band pass windowed view", "[FIR Descs]") {
	const auto desc = Fir.Bandpass.Windowed.Band(winBandLow, winBandHigh).Window(winWindow);
	REQUIRE(desc.lower == winBandLow);
	REQUIRE(desc.upper == winBandHigh);
	REQUIRE(desc.window.size() == 3);
	REQUIRE(desc.window[0] == 1);
}

TEST_CASE("Band pass windowed function", "[FIR Descs]") {
	const auto desc = Fir.Bandpass.Windowed.Band(winBandLow, winBandHigh).Window(windows::hamming);
	REQUIRE(desc.lower == winBandLow);
	REQUIRE(desc.upper == winBandHigh);
	BasicSignal<float, TIME_DOMAIN> window(3);
	desc.window(window);
	REQUIRE(window.size() == 3);
	REQUIRE(Max(window) == Approx(1));
}

TEST_CASE("Band stop windowed view", "[FIR Descs]") {
	const auto desc = Fir.Bandstop.Windowed.Band(winBandLow, winBandHigh).Window(winWindow);
	REQUIRE(desc.lower == winBandLow);
	REQUIRE(desc.upper == winBandHigh);
	REQUIRE(desc.window.size() == 3);
	REQUIRE(desc.window[0] == 1);
}

TEST_CASE("Band stop windowed function", "[FIR Descs]") {
	const auto desc = Fir.Bandstop.Windowed.Band(winBandLow, winBandHigh).Window(windows::hamming);
	REQUIRE(desc.lower == winBandLow);
	REQUIRE(desc.upper == winBandHigh);
	BasicSignal<float, TIME_DOMAIN> window(3);
	desc.window(window);
	REQUIRE(window.size() == 3);
	REQUIRE(Max(window) == Approx(1));
}

TEST_CASE("Arbitrary windowed view", "[FIR Descs]") {
	const auto desc = Fir.Arbitrary.Windowed.Response([](float) { return 1.f; }).Window(winWindow);
	REQUIRE(desc.responseFunc(0.3f) == Approx(1.0f));
	REQUIRE(desc.window.size() == 3);
	REQUIRE(desc.window[0] == 1);
}

TEST_CASE("Arbitrary windowed function", "[FIR Descs]") {
	const auto desc = Fir.Arbitrary.Windowed.Response([](float) { return 1.f; }).Window(windows::blackman);
	REQUIRE(desc.responseFunc(0.3f) == Approx(1.0f));
	BasicSignal<float, TIME_DOMAIN> window(3);
	desc.window(window);
	REQUIRE(window.size() == 3);
	REQUIRE(Max(window) == Approx(1));
}

TEST_CASE("Hilbert windowed view", "[FIR Descs]") {
	const auto desc = Fir.Hilbert.Windowed.Window(winWindow);
	REQUIRE(desc.window.size() == 3);
	REQUIRE(desc.window[0] == 1);
}

TEST_CASE("Hilbert windowed function", "[FIR Descs]") {
	const auto desc = Fir.Hilbert.Windowed.Window(windows::blackman);
	BasicSignal<float, TIME_DOMAIN> window(3);
	desc.window(window);
	REQUIRE(window.size() == 3);
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
constexpr size_t lsGrid = 234;

TEST_CASE("Fresh least squares", "[FIR Descs]") {
	const auto lpc = Fir.Lowpass.LeastSquares.Cutoff(lsBegin1, lsEnd1);
	REQUIRE(lpc.cutoffBegin == lsBegin1);
	REQUIRE(lpc.cutoffEnd == lsEnd1);

	const auto lpw = Fir.Lowpass.LeastSquares.Weight(lsWeightLow, lsWeightTr1, lsWeightHigh);
	REQUIRE(lpw.weightLow == lsWeightLow);
	REQUIRE(lpw.weightTransition == lsWeightTr1);
	REQUIRE(lpw.weightHigh == lsWeightHigh);

	const auto bpc = Fir.Bandpass.LeastSquares.Band(lsBegin1, lsEnd1, lsBegin2, lsEnd2);
	REQUIRE(bpc.lowerBegin == lsBegin1);
	REQUIRE(bpc.lowerEnd == lsEnd1);
	REQUIRE(bpc.upperBegin == lsBegin2);
	REQUIRE(bpc.upperEnd == lsEnd2);

	const auto bpw = Fir.Bandpass.LeastSquares.Weight(lsWeightLow, lsWeightTr1, lsWeightMid, lsWeightTr2, lsWeightHigh);
	REQUIRE(bpw.weightLow == lsWeightLow);
	REQUIRE(bpw.weightTransition1 == lsWeightTr1);
	REQUIRE(bpw.weightMid == lsWeightMid);
	REQUIRE(bpw.weightTransition2 == lsWeightTr2);
	REQUIRE(bpw.weightHigh == lsWeightHigh);
}

TEST_CASE("Low pass least squares", "[FIR Descs]") {
	const auto desc = Fir.Lowpass.LeastSquares.Cutoff(lsBegin1, lsEnd1).Weight(lsWeightLow, lsWeightTr1, lsWeightHigh).Grid(lsGrid);
	REQUIRE(desc.cutoffBegin == lsBegin1);
	REQUIRE(desc.cutoffEnd == lsEnd1);
	REQUIRE(desc.weightLow == lsWeightLow);
	REQUIRE(desc.weightTransition == lsWeightTr1);
	REQUIRE(desc.weightHigh == lsWeightHigh);
	REQUIRE(desc.grid == lsGrid);
}

TEST_CASE("High pass least squares", "[FIR Descs]") {
	const auto desc = Fir.Highpass.LeastSquares.Cutoff(lsBegin1, lsEnd1).Weight(lsWeightLow, lsWeightTr1, lsWeightHigh).Grid(lsGrid);
	REQUIRE(desc.cutoffBegin == lsBegin1);
	REQUIRE(desc.cutoffEnd == lsEnd1);
	REQUIRE(desc.weightLow == lsWeightLow);
	REQUIRE(desc.weightTransition == lsWeightTr1);
	REQUIRE(desc.weightHigh == lsWeightHigh);
	REQUIRE(desc.grid == lsGrid);
}

TEST_CASE("Band pass least squares", "[FIR Descs]") {
	const auto desc = Fir.Bandpass.LeastSquares.Band(lsBegin1, lsEnd1, lsBegin2, lsEnd2).Weight(lsWeightLow, lsWeightTr1, lsWeightMid, lsWeightTr2, lsWeightHigh).Grid(lsGrid);
	REQUIRE(desc.lowerBegin == lsBegin1);
	REQUIRE(desc.lowerEnd == lsEnd1);
	REQUIRE(desc.upperBegin == lsBegin2);
	REQUIRE(desc.upperEnd == lsEnd2);
	REQUIRE(desc.weightLow == lsWeightLow);
	REQUIRE(desc.weightTransition1 == lsWeightTr1);
	REQUIRE(desc.weightMid == lsWeightMid);
	REQUIRE(desc.weightTransition2 == lsWeightTr2);
	REQUIRE(desc.weightHigh == lsWeightHigh);
	REQUIRE(desc.grid == lsGrid);
}

TEST_CASE("Band stop least squares", "[FIR Descs]") {
	const auto desc = Fir.Bandstop.LeastSquares.Band(lsBegin1, lsEnd1, lsBegin2, lsEnd2).Weight(lsWeightLow, lsWeightTr1, lsWeightMid, lsWeightTr2, lsWeightHigh).Grid(lsGrid);
	REQUIRE(desc.lowerBegin == lsBegin1);
	REQUIRE(desc.lowerEnd == lsEnd1);
	REQUIRE(desc.upperBegin == lsBegin2);
	REQUIRE(desc.upperEnd == lsEnd2);
	REQUIRE(desc.weightLow == lsWeightLow);
	REQUIRE(desc.weightTransition1 == lsWeightTr1);
	REQUIRE(desc.weightMid == lsWeightMid);
	REQUIRE(desc.weightTransition2 == lsWeightTr2);
	REQUIRE(desc.weightHigh == lsWeightHigh);
	REQUIRE(desc.grid == lsGrid);
}

TEST_CASE("Arbitrary least squares", "[FIR Descs]") {
	const auto desc = Fir.Arbitrary.LeastSquares.Response([](float) { return 1.f; }).Weight([](float f) { return f < 0.5f ? 1.0f : 0.5f; }).Grid(lsGrid);
	REQUIRE(desc.responseFunc(0.3f) == Approx(1.0f));
	REQUIRE(desc.weightFunc(0.4f) == Approx(1.0f));
	REQUIRE(desc.weightFunc(0.6f) == Approx(0.5f));
	REQUIRE(desc.grid == lsGrid);
}

TEST_CASE("Hilbert least squares", "[FIR Descs]") {
	const auto desc = Fir.Hilbert.LeastSquares.TransitionWidth(0.95f).TransitionWeight(0.3f).Grid(lsGrid);
	REQUIRE(desc.transitionWidth == Approx(0.95f));
	REQUIRE(desc.transitionWeight == Approx(0.3f));
	REQUIRE(desc.grid == lsGrid);
}