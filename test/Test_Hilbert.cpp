#include <catch2/catch.hpp>
#include <dspbb/Filtering/Hilbert.hpp>
#include <dspbb/Filtering/Interpolation.hpp>

using namespace dspbb;


TEST_CASE("Form Type III", "[Hilbert]") {
	const auto filter = HilbertFirWinIII<float, TIME_DOMAIN>(247, dspbb::windows::hamming);
	REQUIRE(filter.Size() == 247);
	const auto nonZeroSamples = Decimate(filter, 2);
	const auto zeroSamples = Decimate(AsView(filter).SubSignal(1), 2);
	REQUIRE(Max(zeroSamples) == 0.0f);
	REQUIRE(Min(Abs(nonZeroSamples)) > 0.0f);
	const auto firstHalf = AsView(nonZeroSamples).SubSignal(0, nonZeroSamples.Size() / 2);
	const auto secondHalf = AsView(nonZeroSamples).SubSignal(nonZeroSamples.Size() / 2);
	REQUIRE(Max(firstHalf) < 0.0f);
	REQUIRE(Min(secondHalf) > 0.0f);
}