#include "dspbb/Math/Statistics.hpp"

#include <catch2/catch.hpp>
#include <dspbb/Filtering/IIR.hpp>
#include <dspbb/Math/DotProduct.hpp>
#include <dspbb/Math/FFT.hpp>
#include <iostream>

using namespace dspbb;

TEST_CASE("IIR test", "[IIR]") {
	constexpr int order = 12;
	const auto system = Butterworth<float>(order);
	const auto discrete = BilinearTransform(system, 1.0f);
	const auto tf = TransferFunction(discrete);

	for (auto& pole : system.Poles()) {
		std::cout << std::arg(pole) / pi_v<float> << std::endl;
	}

	Signal<float, TIME_DOMAIN> dirac(250, 0.0f);
	dirac[0] = 1.0f;

	DirectForm1<float> state;
	state.poleState.Resize(tf.Denominator().size() - 1, 0.0f);
	state.zeroState.Resize(tf.Numerator().size(), 0.0f);

	const auto out = Filter(dirac, tf, state);
	const auto normalization = Sum(out);
	auto padded = out;
	padded.Resize(2048, 0.0f);
	const auto spectrum = Abs(FourierTransform(padded, false));

	REQUIRE(tf.Denominator().size() == order + 1);
	REQUIRE(normalization == Approx(1.0f));
	REQUIRE(Max(spectrum) == Approx(1.0f));
}