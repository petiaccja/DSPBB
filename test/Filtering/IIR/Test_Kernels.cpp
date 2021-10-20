#include <catch2/catch.hpp>
#include <dspbb/Filtering/IIR/Butterworth.hpp>

using namespace dspbb;

TEST_CASE("Butterworth even", "[IIR Kernels]") {
	auto tf = Butterworth<float>(6);
	REQUIRE(tf.gain == 1.0f);
	REQUIRE(tf.zeros.NumRoots() == 0);
	REQUIRE(tf.poles.NumRoots() == 6);
	REQUIRE(std::arg(tf.poles.ComplexRoots()[0]) == Approx(pi_v<float> * 0.58333333f));
	REQUIRE(std::arg(tf.poles.ComplexRoots()[1]) == Approx(pi_v<float> * 0.75f));
	REQUIRE(std::arg(tf.poles.ComplexRoots()[2]) == Approx(pi_v<float> * 0.91666666f));
}

TEST_CASE("Butterworth odd", "[IIR Kernels]") {
	auto tf = Butterworth<float>(7);
	REQUIRE(tf.gain == 1.0f);
	REQUIRE(tf.zeros.NumRoots() == 0);
	REQUIRE(tf.poles.NumRoots() == 7);
	REQUIRE(std::arg(tf.poles.RealRoots()[0]) == Approx(pi_v<float> * 1.0f));
	REQUIRE(std::arg(tf.poles.ComplexRoots()[0]) == Approx(pi_v<float> * 0.57142857142857f));
	REQUIRE(std::arg(tf.poles.ComplexRoots()[1]) == Approx(pi_v<float> * 0.71428571428572f));
	REQUIRE(std::arg(tf.poles.ComplexRoots()[2]) == Approx(pi_v<float> * 0.85714285714286f));
}