#include "../TestUtils.hpp"

#include <dspbb/Generators/Spaces.hpp>
#include <dspbb/Math/JacobiFunctions.hpp>

#include <catch2/catch.hpp>

using namespace dspbb;
using namespace std::complex_literals;


TEST_CASE("Shift scalar", "[Jacobi functions]") {
	std::array<float, 5> inputs = { -1234567891011.0f, -23.15f, 0.86f, 3.14f, 987654321.0f };
	std::array<float, 5> remainders = { 0.0f, -0.15f, -0.14f, 0.14f, 0.0f };
	std::array<int, 5> counts = { 0, 23, -1, -3, 0 };
	for (int i = 0; i < 5; ++i) {
		const auto [remainder, count] = ShiftScalar(inputs[i]);
		REQUIRE(remainder == Approx(remainders[i]).epsilon(1e-5f));
		REQUIRE(count % 8 == counts[i] % 8);
	}
}

TEST_CASE("Shift variant", "[Jacobi functions]") {
	REQUIRE(ShiftVariant(1, 5) == 1);
	REQUIRE(ShiftVariant(1, 6) == 1);
	REQUIRE(ShiftVariant(2, 5) == 2);
	REQUIRE(ShiftVariant(2, 6) == 2);
	REQUIRE(ShiftVariant(3, 5) == 4);
	REQUIRE(ShiftVariant(3, 6) == 3);
	REQUIRE(ShiftVariant(4, 5) == 3);
	REQUIRE(ShiftVariant(4, 6) == 4);
}

TEST_CASE("Shift multiplier", "[Jacobi functions]") {
	REQUIRE(ShiftMultiplier<float>(1, 5) == ApproxComplex(std::polar(1.0f, 5.0f / 4.0f * pi_v<float>)));
	REQUIRE(ShiftMultiplier<float>(1, 14) == ApproxComplex(std::polar(1.0f, 6.0f / 4.0f * pi_v<float>)));
	REQUIRE(ShiftMultiplier<float>(2, 5) == ApproxComplex(std::polar(1.0f, 5.0f / 4.0f * pi_v<float>)));
	REQUIRE(ShiftMultiplier<float>(2, 14) == ApproxComplex(std::polar(1.0f, 6.0f / 4.0f * pi_v<float>)));
	REQUIRE(ShiftMultiplier<float>(3, 5) == ApproxComplex(std::polar(1.0f, 0.0f)));
	REQUIRE(ShiftMultiplier<float>(3, 6) == ApproxComplex(std::polar(1.0f, 0.0f)));
	REQUIRE(ShiftMultiplier<float>(4, 5) == ApproxComplex(std::polar(1.0f, 0.0f)));
	REQUIRE(ShiftMultiplier<float>(4, 6) == ApproxComplex(std::polar(1.0f, 0.0f)));
}

TEST_CASE("Shift tau negative direction var 1", "[Jacobi functions]") {
	const std::complex<float> tau = { 12.99345f, 0.1f };
	const int variant = 1;
	const auto [_, newTau, multiplier, __, newVariant] = ShiftTau(tau, variant);
	REQUIRE(real(newTau) == Approx(12.99345f - 13.0f).epsilon(1e-6f));
	REQUIRE(imag(newTau) == Approx(0.1f).epsilon(1e-6f));
	REQUIRE(multiplier == ApproxComplex(std::exp(1if * pi_v<float> * -13.0f / 4.0f)).epsilon(1e-6f));
	REQUIRE(newVariant == variant);
}

TEST_CASE("Shift tau positive direction var 2", "[Jacobi functions]") {
	const std::complex<float> tau = { -12.99345f, 0.1f };
	const int variant = 2;
	const auto [_, newTau, multiplier, __, newVariant] = ShiftTau(tau, variant);
	REQUIRE(real(newTau) == Approx(-12.99345f + 13.0f).epsilon(1e-6f));
	REQUIRE(imag(newTau) == Approx(0.1f).epsilon(1e-6f));
	REQUIRE(multiplier == ApproxComplex(std::exp(1if * pi_v<float> * +13.0f / 4.0f)).epsilon(1e-6f));
	REQUIRE(newVariant == variant);
}

TEST_CASE("Shift tau negative direction var 31", "[Jacobi functions]") {
	const std::complex<float> tau = { 12.99345f, 0.1f };
	const int variant = 3;
	const auto [_, newTau, multiplier, __, newVariant] = ShiftTau(tau, variant);
	REQUIRE(real(newTau) == Approx(12.99345f - 13.0f).epsilon(1e-6f));
	REQUIRE(imag(newTau) == Approx(0.1f).epsilon(1e-6f));
	REQUIRE(multiplier == ApproxComplex(1.0f).epsilon(1e-6f));
	REQUIRE(newVariant == 4);
}

TEST_CASE("Shift tau positive direction var 4", "[Jacobi functions]") {
	const std::complex<float> tau = { -12.99345f, 0.1f };
	const int variant = 4;
	const auto [_, newTau, multiplier, __, newVariant] = ShiftTau(tau, variant);
	REQUIRE(real(newTau) == Approx(-12.99345f + 13.0f).epsilon(1e-6f));
	REQUIRE(imag(newTau) == Approx(0.1f).epsilon(1e-6f));
	REQUIRE(multiplier == ApproxComplex(1.0f).epsilon(1e-6f));
	REQUIRE(newVariant == 3);
}

TEST_CASE("Invert variant") {
	REQUIRE(InvertVariant(1) == 1);
	REQUIRE(InvertVariant(2) == 4);
	REQUIRE(InvertVariant(3) == 3);
	REQUIRE(InvertVariant(4) == 2);
}

TEST_CASE("Evaluate theta_1", "[Jacobi functions]") {
	const std::complex<float> q = -0.1f + 0.3if;
	const auto tau = -i_v<float> * std::log(q) / pi_v<float>;
	const auto q2 = std::exp(i_v<float> * pi_v<float> * tau);
	REQUIRE(q == ApproxComplex(q2));

	const auto z = LinSpace<float, DOMAINLESS>(-5.0f, 5.0f, 400);
	std::vector<float> r;
	std::vector<float> i;
	for (auto& v : z) {
		const auto c = theta_1(std::complex<float>(v, 0.7f), tau);
		r.push_back(real(c));
		i.push_back(imag(c));
	}
	REQUIRE(false);
}