#include <dspbb/Filtering/IIR/Butterworth.hpp>
#include <dspbb/Filtering/IIR/Chebyshev.hpp>

#include <catch2/catch.hpp>

using namespace dspbb;
using namespace std::complex_literals;


TEST_CASE("Butterworth even", "[IIR Kernels]") {
	auto tf = Butterworth<double>(6);
	REQUIRE(tf.gain == 1.0);
	REQUIRE(tf.zeros.NumRoots() == 0);
	REQUIRE(tf.poles.NumRoots() == 6);
	REQUIRE(std::arg(tf.poles.ComplexPairs()[0]) == Approx(pi_v<double> * 0.58333333));
	REQUIRE(std::arg(tf.poles.ComplexPairs()[1]) == Approx(pi_v<double> * 0.75));
	REQUIRE(std::arg(tf.poles.ComplexPairs()[2]) == Approx(pi_v<double> * 0.91666666));
}

TEST_CASE("Butterworth odd", "[IIR Kernels]") {
	auto tf = Butterworth<double>(7);
	REQUIRE(tf.gain == 1.0);
	REQUIRE(tf.zeros.NumRoots() == 0);
	REQUIRE(tf.poles.NumRoots() == 7);
	REQUIRE(std::arg(tf.poles.RealRoots()[0]) == Approx(pi_v<double> * 1.0));
	REQUIRE(std::arg(tf.poles.ComplexPairs()[0]) == Approx(pi_v<double> * 0.57142857142857));
	REQUIRE(std::arg(tf.poles.ComplexPairs()[1]) == Approx(pi_v<double> * 0.71428571428572));
	REQUIRE(std::arg(tf.poles.ComplexPairs()[2]) == Approx(pi_v<double> * 0.85714285714286));
}

TEST_CASE("Chebyshev type I even", "[IIR Kernels]") {
	constexpr auto ripple = 0.005;
	auto tf = Chebyshev1(8, ripple);
	REQUIRE(tf.zeros.NumRoots() == 0);
	REQUIRE(tf.poles.NumRoots() == 8);
	REQUIRE(std::abs(tf(0.0)) == Approx(1.0 - ripple).margin(1e-9));
	REQUIRE(std::abs(tf(1.0i)) == Approx(1.0 - ripple));
}

TEST_CASE("Chebyshev type I odd", "[IIR Kernels]") {
	constexpr auto ripple = 0.005;
	auto tf = Chebyshev1(9, ripple);
	REQUIRE(tf.zeros.NumRoots() == 0);
	REQUIRE(tf.poles.NumRoots() == 9);
	REQUIRE(std::abs(tf(0.0)) == Approx(1.0).margin(1e-9));
	REQUIRE(std::abs(tf(1.0i)) == Approx(1.0 - ripple));
}

TEST_CASE("Chebyshev type II even", "[IIR Kernels]") {
	constexpr auto ripple = 0.005;
	auto tf = Chebyshev2(8, ripple);
	REQUIRE(tf.zeros.NumRoots() == 8);
	REQUIRE(tf.poles.NumRoots() == 8);
	REQUIRE(std::abs(tf(0.0)) == Approx(1.0).margin(1e-9));
	REQUIRE(std::abs(tf(1.0i)) == Approx(ripple).margin(1e-9));
}

TEST_CASE("Chebyshev type II odd", "[IIR Kernels]") {
	constexpr auto ripple = 0.05;
	auto tf = Chebyshev2(9, ripple);
	REQUIRE(tf.zeros.NumRoots() == 8);
	REQUIRE(tf.poles.NumRoots() == 9);
	REQUIRE(std::abs(tf(0.0)) == Approx(1.0).margin(1e-9));
	REQUIRE(std::abs(tf(1.0i)) == Approx(ripple).margin(1e-9));
}