#include "dspbb/Math/Statistics.hpp"

#include <catch2/catch.hpp>
#include <dspbb/Filtering/IIR.hpp>
#include <dspbb/Math/FFT.hpp>
#include <iostream>

using namespace dspbb;

TEST_CASE("IIR test", "[IIR]") {
	constexpr int order = 7;
	const auto butter = IirFilter<float>(order, Lowpass(BUTTERWORTH).Cutoff(0.32f));
	const auto tf = TransferFunction(butter);
	
	Signal<float, TIME_DOMAIN> dirac(250, 0.0f);
	dirac[0] = 1.0f;

	DirectFormI<float> state{order};

	const auto out = Filter(dirac, tf, state);
	const auto normalization = Sum(out);
	auto padded = out;
	padded.Resize(2048, 0.0f);
	const auto spectrum = Abs(FourierTransform(padded, false));

	REQUIRE(tf.denominator.Size() == order + 1);
	REQUIRE(normalization == Approx(1.0f));
	REQUIRE(Max(spectrum) == Approx(1.0f));
}