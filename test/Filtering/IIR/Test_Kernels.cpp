#include <dspbb/Filtering/IIR/Butterworth.hpp>
#include <dspbb/Filtering/IIR/Chebyshev.hpp>
#include <dspbb/Filtering/IIR/Elliptic.hpp>

#include <catch2/catch.hpp>

using namespace dspbb;
using namespace std::complex_literals;


TEST_CASE("Butterworth even", "[IIR Kernels]") {
	auto tf = Butterworth<double>(6);
	REQUIRE(tf.gain == 1.0);
	REQUIRE(tf.zeros.num_roots() == 0);
	REQUIRE(tf.poles.num_roots() == 6);
	REQUIRE(std::arg(tf.poles.complex_pairs()[0]) == Approx(pi_v<double> * 0.58333333));
	REQUIRE(std::arg(tf.poles.complex_pairs()[1]) == Approx(pi_v<double> * 0.75));
	REQUIRE(std::arg(tf.poles.complex_pairs()[2]) == Approx(pi_v<double> * 0.91666666));
}

TEST_CASE("Butterworth odd", "[IIR Kernels]") {
	auto tf = Butterworth<double>(7);
	REQUIRE(tf.gain == 1.0);
	REQUIRE(tf.zeros.num_roots() == 0);
	REQUIRE(tf.poles.num_roots() == 7);
	REQUIRE(std::arg(tf.poles.real_roots()[0]) == Approx(pi_v<double> * 1.0));
	REQUIRE(std::arg(tf.poles.complex_pairs()[0]) == Approx(pi_v<double> * 0.57142857142857));
	REQUIRE(std::arg(tf.poles.complex_pairs()[1]) == Approx(pi_v<double> * 0.71428571428572));
	REQUIRE(std::arg(tf.poles.complex_pairs()[2]) == Approx(pi_v<double> * 0.85714285714286));
}

TEST_CASE("Chebyshev type I even", "[IIR Kernels]") {
	constexpr auto ripple = 0.005;
	for (auto order : { 2, 4, 6, 8 }) {
		auto tf = Chebyshev1(order, ripple);
		REQUIRE(tf.zeros.num_roots() == 0);
		REQUIRE(tf.poles.num_roots() == order);
		REQUIRE(std::abs(tf(0.0)) == Approx(1.0 - ripple).margin(1e-9));
		REQUIRE(std::abs(tf(1.0i)) == Approx(1.0 - ripple));
	}
}

TEST_CASE("Chebyshev type I odd", "[IIR Kernels]") {
	constexpr auto ripple = 0.005;
	for (auto order : { 1, 3, 5, 7 }) {
		auto tf = Chebyshev1(order, ripple);
		REQUIRE(tf.zeros.num_roots() == 0);
		REQUIRE(tf.poles.num_roots() == order);
		REQUIRE(std::abs(tf(0.0)) == Approx(1.0).margin(1e-9));
		REQUIRE(std::abs(tf(1.0i)) == Approx(1.0 - ripple));
	}
}

TEST_CASE("Chebyshev type II even", "[IIR Kernels]") {
	constexpr auto ripple = 0.005;
	for (auto order : { 2, 4, 6, 8 }) {
		auto tf = Chebyshev2(order, ripple);
		REQUIRE(tf.zeros.num_roots() == order);
		REQUIRE(tf.poles.num_roots() == order);
		REQUIRE(std::abs(tf(0.0)) == Approx(1.0).margin(1e-9));
		REQUIRE(std::abs(tf(1.0i)) == Approx(ripple).margin(1e-9));
	}
}

TEST_CASE("Chebyshev type II odd", "[IIR Kernels]") {
	constexpr auto ripple = 0.05;
	for (auto order : { 1, 3, 5, 7 }) {
		auto tf = Chebyshev2(order, ripple);
		REQUIRE(tf.zeros.num_roots() == order - 1);
		REQUIRE(tf.poles.num_roots() == order);
		REQUIRE(std::abs(tf(0.0)) == Approx(1.0).margin(1e-9));
		REQUIRE(std::abs(tf(1.0i)) == Approx(ripple).margin(1e-9));
	}
}

TEST_CASE("Elliptic even", "[IIR Kernels]") {
	constexpr auto passbandRipple = 0.05;
	constexpr auto stopbandRipple = 0.1;

	for (auto order : { 2, 4, 6, 8 }) {
		const auto [k, kp, K, Kp, k1, k1p, K1, K1p, epsilon] = impl::EllipticOrderRipples(order, passbandRipple, stopbandRipple);
		const auto ws = 1.0 / k;

		auto tf = Elliptic(order, passbandRipple, stopbandRipple);

		REQUIRE(tf.zeros.num_roots() == order);
		REQUIRE(tf.poles.num_roots() == order);
		REQUIRE(tf(0.0) == Approx(1.0 - passbandRipple).margin(1e-9));
		REQUIRE(abs(tf(1.0i)) == Approx(1.0 - passbandRipple).margin(1e-9));
		REQUIRE(abs(tf(1.0i * ws)) == Approx(stopbandRipple).margin(1e-9));
	}
}

TEST_CASE("Elliptic odd", "[IIR Kernels]") {
	constexpr auto passbandRipple = 0.05;
	constexpr auto stopbandRipple = 0.1;

	for (auto order : { 1, 3, 5, 7 }) {
		const auto [k, kp, K, Kp, k1, k1p, K1, K1p, epsilon] = impl::EllipticOrderRipples(order, passbandRipple, stopbandRipple);
		const auto ws = 1.0 / k;

		auto tf = Elliptic(order, passbandRipple, stopbandRipple);

		REQUIRE(tf.zeros.num_roots() == order - 1);
		REQUIRE(tf.poles.num_roots() == order);
		REQUIRE(tf(0.0) == Approx(1.0).margin(1e-9));
		REQUIRE(abs(tf(1.0i)) == Approx(1.0 - passbandRipple).margin(1e-9));
		REQUIRE(abs(tf(1.0i * ws)) == Approx(stopbandRipple).margin(1e-9));
	}
}