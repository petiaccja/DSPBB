#include <catch2/catch.hpp>
#include <dspbb/LTISystems/DiscretizationTransforms.hpp>
#include <dspbb/Utility/Numbers.hpp>

using namespace dspbb;
using namespace std::complex_literals;


TEST_CASE("Bilinear C->D", "[DiscreteizationTransforms]") {
	constexpr float sampleRate = 6.0f;
	const ContinuousPoleZeroSystem<float> c{
		1.5f,
		{ -10000.f, -2.3f + 0.3if, -2.3f - 0.3if },
		{ -0.0f + 0.7if, -0.0f - 0.7if },
	};
	const DiscretePoleZeroSystem<float> d = BilinearTransform(c, sampleRate);
	// Number of poles and zeros.
	REQUIRE(d.poles.RealRoots().Size() == 1);
	REQUIRE(d.poles.ComplexRoots().Size() == 1);
	REQUIRE(d.zeros.RealRoots().Size() == 1);
	REQUIRE(d.zeros.ComplexRoots().Size() == 1);
	// -INF s maps to -1 z
	REQUIRE(std::real(d.zeros.RealRoots()[0]) == Approx(-1).margin(0.01));
	// Points on the left plane map to points inside the unit circle.
	REQUIRE(std::abs(d.zeros.ComplexRoots()[0]) < 0.9f);
	// jw axis maps to unit circle.
	REQUIRE(std::abs(d.poles.ComplexRoots()[0]) == Approx(1));
	// Frequency warping is correct.
	REQUIRE(std::arg(d.poles.ComplexRoots()[0]) * sampleRate == Approx(2 * sampleRate * std::atan(c.poles.ComplexRoots()[0].imag() / sampleRate / 2)));
}


TEST_CASE("Bilinear C->D Prewarp", "[DiscreteizationTransforms]") {
	constexpr float sampleRate = 6.0f;
	constexpr float nyquistLimit = sampleRate / 2.0f;
	constexpr float angularLimit = 2 * pi_v<float> * nyquistLimit;
	constexpr float cutoff = 0.65f;
	const ContinuousPoleZeroSystem<float> c{
		1.5f,
		{ -2.3f + 0.3if, -2.3f - 0.3if },
		{ -0.0f + 1if * angularLimit * cutoff, -0.0f - 1if * angularLimit * cutoff },
	};
	const DiscretePoleZeroSystem<float> d = BilinearTransform(c, sampleRate, { cutoff * angularLimit });

	REQUIRE(arg(d.poles.ComplexRoots()[0]) == Approx(cutoff * pi_v<float>));
}