#include <catch2/catch.hpp>
#include <dspbb/Filtering/FilterParameters.hpp>
#include <dspbb/Filtering/IIR.hpp>
#include <dspbb/Math/FFT.hpp>
#include <dspbb/Math/Statistics.hpp>
#include <iostream>

using namespace dspbb;

template <class T>
bool IsStable(const DiscretePoleZeroSystem<T>& system) {
	std::vector<T> lengths;
	std::transform(system.poles.RealRoots().begin(), system.poles.RealRoots().end(), std::back_inserter(lengths), [](const auto& arg) { return std::abs(arg); });
	std::transform(system.poles.ComplexRoots().begin(), system.poles.ComplexRoots().end(), std::back_inserter(lengths), [](const auto& arg) { return std::abs(arg); });
	return std::all_of(lengths.begin(), lengths.end(), [](auto len) { return len < 1.0f; });
}

TEST_CASE("IIR test", "[IIR]") {
	constexpr int order = 3;
	const auto butter = IirFilter<float>(order, Lowpass(BUTTERWORTH).Cutoff(0.5f));
	const auto butter2 = Halfband2Lowpass(butter, 0.5f);
	const auto tf = TransferFunction(butter2);

	Signal<float, TIME_DOMAIN> dirac(250, 0.0f);
	dirac[0] = 1.0f;

	DirectFormI<float> state{ butter2.poles.NumRoots() };

	const auto out = Filter(dirac, tf, state);
	const auto normalization = Sum(out);
	auto padded = out;
	padded.Resize(2048, 0.0f);
	const auto spectrum = Abs(FourierTransform(padded, false));

	const auto [amplitude, phase] = FrequencyResponse(tf, 1024);

	REQUIRE(IsStable(butter2));
	REQUIRE(tf.denominator.Size() == order + 1);
	REQUIRE(normalization == Approx(1.0f));
	REQUIRE(Max(spectrum) == Approx(1.0f));
}