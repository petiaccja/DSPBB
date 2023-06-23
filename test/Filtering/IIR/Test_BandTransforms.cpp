#include <dspbb/Filtering/IIR/BandTransforms.hpp>
#include <dspbb/Filtering/MeasureFilter.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace dspbb;
using namespace std::complex_literals;
using Catch::Approx;


const DiscreteZeroPoleGain<float> prototype = {
	0.166666672f,
	{ -1.f, -1.f, -1.f },
	{ 0.f, 0.f + 0.577350259if, 0.f - 0.577350259if }
};
constexpr float minus3dB = 0.70794578438414f;


template <class Iter, class T>
std::optional<T> CrossoverPoint(Iter first, Iter last, T crossoverLevel) {
	const auto crossoverIt = std::adjacent_find(first, last, [crossoverLevel](const auto& v1, const auto& v2) {
		return (v1 - crossoverLevel) * (v2 - crossoverLevel) <= T(0); // Opposite signs or zero
	});
	if (crossoverIt != last) {
		const size_t crossoverIndex = crossoverIt - first;
		auto nextIt = crossoverIt;
		++nextIt;
		const T d1 = std::abs(*crossoverIt - crossoverLevel);
		const T d2 = std::abs(*nextIt - crossoverLevel);
		return { T(crossoverIndex) + d1 / (d1 + d2) };
	}
	return {};
}


TEST_CASE("Verify prototype filter", "[IIR band transforms]") {
	const auto [amplitude, phase] = FrequencyResponse(prototype, 1024);
	const float crossover = CrossoverPoint(amplitude.begin(), amplitude.end(), minus3dB).value() / float(amplitude.size() - 1);
	REQUIRE(crossover == Approx(0.5).epsilon(5e-3f));
}

TEST_CASE("Lowpass to lowpass", "[IIR band transforms]") {
	const auto lp = Halfband2Lowpass(prototype, 0.3f);
	const auto [amplitude, phase] = FrequencyResponse(lp, 1024);
	const float crossover = CrossoverPoint(amplitude.begin(), amplitude.end(), minus3dB).value() / float(amplitude.size() - 1);
	REQUIRE(crossover == Approx(0.3).epsilon(5e-3f));
	REQUIRE_NOTHROW(MeasureLowpassFilter(amplitude));
}

TEST_CASE("Lowpass to highpass", "[IIR band transforms]") {
	const auto hp = Halfband2Highpass(prototype, 0.4f);
	const auto [amplitude, phase] = FrequencyResponse(hp, 1024);
	const float crossover = CrossoverPoint(amplitude.begin(), amplitude.end(), minus3dB).value() / float(amplitude.size() - 1);
	REQUIRE(crossover == Approx(0.4).epsilon(5e-3f));
	REQUIRE_NOTHROW(MeasureHighpassFilter(amplitude));
}

TEST_CASE("Lowpass to bandpass", "[IIR band transforms]") {
	const auto bp = Halfband2Bandpass(prototype, 0.35f, 0.6f);
	const auto [amplitude, phase] = FrequencyResponse(bp, 1024);
	const float crossover1 = CrossoverPoint(amplitude.begin(), amplitude.end(), minus3dB).value() / float(amplitude.size() - 1);
	const float crossover2 = (500 + CrossoverPoint(amplitude.begin() + 500, amplitude.end(), minus3dB).value()) / float(amplitude.size() - 1);
	REQUIRE(crossover1 == Approx(0.35).epsilon(5e-3f));
	REQUIRE(crossover2 == Approx(0.6).epsilon(5e-3f));
	REQUIRE_NOTHROW(MeasureBandpassFilter(amplitude));
}

TEST_CASE("Lowpass to bandstop", "[IIR band transforms]") {
	const auto hp = Halfband2Bandstop(prototype, 0.45f, 0.65f);
	const auto [amplitude, phase] = FrequencyResponse(hp, 1024);
	const float crossover1 = CrossoverPoint(amplitude.begin(), amplitude.end(), minus3dB).value() / float(amplitude.size() - 1);
	const float crossover2 = (500 + CrossoverPoint(amplitude.begin() + 500, amplitude.end(), minus3dB).value()) / float(amplitude.size() - 1);
	REQUIRE(crossover1 == Approx(0.45).epsilon(5e-3f));
	REQUIRE(crossover2 == Approx(0.65).epsilon(5e-3f));
	REQUIRE_NOTHROW(MeasureBandstopFilter(amplitude));
}