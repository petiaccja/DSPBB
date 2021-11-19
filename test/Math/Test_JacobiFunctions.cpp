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
	REQUIRE(ShiftMultiplier<float>(1, 5) == ApproxComplex(std::polar(1.0f, -5.0f / 4.0f * pi_v<float>)));
	REQUIRE(ShiftMultiplier<float>(1, 14) == ApproxComplex(std::polar(1.0f, -6.0f / 4.0f * pi_v<float>)));
	REQUIRE(ShiftMultiplier<float>(2, 5) == ApproxComplex(std::polar(1.0f, -5.0f / 4.0f * pi_v<float>)));
	REQUIRE(ShiftMultiplier<float>(2, 14) == ApproxComplex(std::polar(1.0f, -6.0f / 4.0f * pi_v<float>)));
	REQUIRE(ShiftMultiplier<float>(3, 5) == ApproxComplex(std::polar(1.0f, 0.0f)));
	REQUIRE(ShiftMultiplier<float>(3, 6) == ApproxComplex(std::polar(1.0f, 0.0f)));
	REQUIRE(ShiftMultiplier<float>(4, 5) == ApproxComplex(std::polar(1.0f, 0.0f)));
	REQUIRE(ShiftMultiplier<float>(4, 6) == ApproxComplex(std::polar(1.0f, 0.0f)));
}

TEST_CASE("Shift tau negative direction var 1", "[Jacobi functions]") {
	const std::complex<float> tau = { 12.99345f, 0.1f };
	const int variant = 1;
	const auto [_, newTau, multiplier, __, newVariant] = ShiftTau({}, tau, variant);
	REQUIRE(real(newTau) == Approx(12.99345f - 13.0f).epsilon(1e-6f));
	REQUIRE(imag(newTau) == Approx(0.1f).epsilon(1e-6f));
	REQUIRE(multiplier == ApproxComplex(std::exp(1if * pi_v<float> * 13.0f / 4.0f)).epsilon(1e-6f));
	REQUIRE(newVariant == variant);
}

TEST_CASE("Shift tau positive direction var 2", "[Jacobi functions]") {
	const std::complex<float> tau = { -12.99345f, 0.1f };
	const int variant = 2;
	const auto [_, newTau, multiplier, __, newVariant] = ShiftTau({}, tau, variant);
	REQUIRE(real(newTau) == Approx(-12.99345f + 13.0f).epsilon(1e-6f));
	REQUIRE(imag(newTau) == Approx(0.1f).epsilon(1e-6f));
	REQUIRE(multiplier == ApproxComplex(std::exp(1if * pi_v<float> * -13.0f / 4.0f)).epsilon(1e-6f));
	REQUIRE(newVariant == variant);
}

TEST_CASE("Shift tau negative direction var 3", "[Jacobi functions]") {
	const std::complex<float> tau = { 12.99345f, 0.1f };
	const int variant = 3;
	const auto [_, newTau, multiplier, __, newVariant] = ShiftTau({} , tau, variant);
	REQUIRE(real(newTau) == Approx(12.99345f - 13.0f).epsilon(1e-6f));
	REQUIRE(imag(newTau) == Approx(0.1f).epsilon(1e-6f));
	REQUIRE(multiplier == ApproxComplex(1.0f).epsilon(1e-6f));
	REQUIRE(newVariant == 4);
}

TEST_CASE("Shift tau positive direction var 4", "[Jacobi functions]") {
	const std::complex<float> tau = { -12.99345f, 0.1f };
	const int variant = 4;
	const auto [_, newTau, multiplier, __, newVariant] = ShiftTau({}, tau, variant);
	REQUIRE(real(newTau) == Approx(-12.99345f + 13.0f).epsilon(1e-6f));
	REQUIRE(imag(newTau) == Approx(0.1f).epsilon(1e-6f));
	REQUIRE(multiplier == ApproxComplex(1.0f).epsilon(1e-6f));
	REQUIRE(newVariant == 3);
}

TEST_CASE("Invert variant", "[Jacobi functions]") {
	REQUIRE(InvertVariant(1) == 1);
	REQUIRE(InvertVariant(2) == 4);
	REQUIRE(InvertVariant(3) == 3);
	REQUIRE(InvertVariant(4) == 2);
}

TEST_CASE("Invert multiplier #1", "[Jacobi functions]") {
	const auto z = 2.0f + 0.0if;
	const auto tau = 0.5if;
	const auto [factor, exponent] = InvertMultiplier(z, tau, 1);
	REQUIRE(factor == ApproxComplex(-i_v<float> * std::sqrt(2.0f)));
	REQUIRE(exponent == ApproxComplex(-8.0f / pi_v<float>));
}

TEST_CASE("Invert multiplier #2,3,4", "[Jacobi functions]") {
	const auto z = 2.0f + 0.0if;
	const auto tau = 0.5if;
	for (int variant = 2; variant <= 4; ++variant) {
		const auto [factor, exponent] = InvertMultiplier(z, tau, variant);
		REQUIRE(factor == ApproxComplex(std::sqrt(2.0f)));
		REQUIRE(exponent == ApproxComplex(-8.0f / pi_v<float>));
	}
}

TEST_CASE("Invert tau #1", "[Jacobi functions]") {
	const auto z = 2.0f + 0.0if;
	const auto tau = 0.5if;
	const auto [newZ, newTau, factor, exponent, newVariant] = InvertTau(z, tau, 1);
	REQUIRE(newZ == ApproxComplex(4if));
	REQUIRE(newTau == ApproxComplex(2if));
	REQUIRE(factor == ApproxComplex(-i_v<float> * std::sqrt(2.0f)));
	REQUIRE(exponent == ApproxComplex(-8.0f / pi_v<float>));
	REQUIRE(newVariant == 1);
}

TEST_CASE("Invert tau #2,3,4", "[Jacobi functions]") {
	const auto z = 2.0f + 0.0if;
	const auto tau = 0.5if;
	for (int variant = 2; variant <= 4; ++variant) {
		const auto [newZ, newTau, factor, exponent, newVariant] = InvertTau(z, tau, variant);
		REQUIRE(newZ == ApproxComplex(4if));
		REQUIRE(newTau == ApproxComplex(2if));
		REQUIRE(factor == ApproxComplex(std::sqrt(2.0f)));
		REQUIRE(exponent == ApproxComplex(-8.0f / pi_v<float>));
		REQUIRE(newVariant == 2 + 4 - variant);
	}
}

TEST_CASE("Lattice shift identity", "[Jacobi functions]") {
	const auto z = 0.7f + 0.3if;
	const auto tau = -1.55f + 0.4if;
	const auto shiftedTau = tau + 1.0f;
	REQUIRE(theta(1, z, shiftedTau, 0) == ApproxComplex(std::exp(pi_v<float> * i_v<float> / 4.0f) * theta(1, z, tau, 0)));
	REQUIRE(theta(2, z, shiftedTau, 0) == ApproxComplex(std::exp(pi_v<float> * i_v<float> / 4.0f) * theta(2, z, tau, 0)));
	REQUIRE(theta(3, z, shiftedTau, 0) == ApproxComplex(theta(4, z, tau, 0)));
	REQUIRE(theta(4, z, shiftedTau, 0) == ApproxComplex(theta(3, z, tau, 0)));
}

TEST_CASE("Lattice shift identity double shift", "[Jacobi functions]") {
	const auto z = 0.7f + 0.3if;
	const auto tau = -1.55f + 0.4if;
	const auto shiftedTau = tau + 2.0f;

	REQUIRE(theta(1, z, shiftedTau, 0) == ApproxComplex(std::exp(2.0f * pi_v<float> * i_v<float> / 4.0f) * theta(1, z, tau, 0)));
	REQUIRE(theta(2, z, shiftedTau, 0) == ApproxComplex(std::exp(2.0f * pi_v<float> * i_v<float> / 4.0f) * theta(2, z, tau, 0)));
	REQUIRE(theta(3, z, shiftedTau, 0) == ApproxComplex(theta(3, z, tau, 0)));
	REQUIRE(theta(4, z, shiftedTau, 0) == ApproxComplex(theta(4, z, tau, 0)));
}

TEST_CASE("Lattice inversion identity", "[Jacobi functions]") {
	const auto z = 0.7f + 0.3if;
	const auto tau = -1.55f + 0.4if;
	const auto invTau = -1.0f / tau;
	[[maybe_unused]] const auto invQ = std::exp(i_v<float> * invTau * pi_v<float>);

	const auto factor1 = std::sqrt(tau / i_v<float>);
	const auto factor2 = std::exp(i_v<float> * tau * z * z / pi_v<float>);

	REQUIRE(theta(1, z, invTau, 0) == ApproxComplex(-i_v<float> * factor1 * factor2 * theta(1, tau * z, tau, 0)));
	REQUIRE(theta(2, z, invTau, 0) == ApproxComplex(factor1 * factor2 * theta(4, tau * z, tau, 0)));
	REQUIRE(theta(3, z, invTau, 0) == ApproxComplex(factor1 * factor2 * theta(3, tau * z, tau, 0)));
	REQUIRE(theta(4, z, invTau, 0) == ApproxComplex(factor1 * factor2 * theta(2, tau * z, tau, 0)));
}



TEST_CASE("Shift tau identity", "[Jacobi functions]") {
	const auto z = 0.7f + 0.3if;

	for (int shift : { 0, 1, 3, 4, -7, -8, 11 }) {
		const auto tau = 0.45f + float(shift) + 0.4if;
		for (int variant = 1; variant <= 4; ++variant) {
			const auto [newZ, newTau, factor, exponent, newVariant] = ShiftTau(z, tau, variant);
			REQUIRE(theta(variant, z, tau, 0) == ApproxComplex(factor * theta(newVariant, newZ, newTau, 0)));
		}
	}
}


TEST_CASE("Invert tau identity", "[Jacobi functions]") {
	const auto z = 0.7 + 0.3i;
	const std::array taus = {
		0.73 + 1.49i,
		-3.17 + 0.49i,
		0.77i,
		-0.11 + 1.03i,
	};

	for (const auto& tau : taus) {
		for (int variant = 1; variant <= 4; ++variant) {
			UNSCOPED_INFO("variant=" << variant << "  tau=" << tau);
			const auto [newZ, newTau, factor, exponent, newVariant] = InvertTau(z, tau, variant);
			const auto trial = factor * std::exp(exponent) * theta(newVariant, newZ, newTau, 0);
			const auto control = theta(variant, z, tau, 0);
			REQUIRE(control == ApproxComplex(trial).epsilon(1e-3f));
		}
	}
}


TEST_CASE("Rotate tau identity", "[Jacobi functions]") {
	const auto z = 0.7 + 0.3i;
	const std::array taus = {
		0.73 + 1.49i,
		-3.17 + 0.49i,
		0.77i,
		-0.11 + 1.03i,
	};

	for (const auto& tau : taus) {
		for (int variant = 1; variant <= 4; ++variant) {
			UNSCOPED_INFO("variant=" << variant << "  tau=" << tau);
			const auto [newZ, newTau, factor, exponent, newVariant] = RotateTau(z, tau, variant);
			const auto trial = factor * std::exp(exponent) * theta(newVariant, newZ, newTau, 0);
			const auto control = theta(variant, z, tau, 0);
			REQUIRE(control == ApproxComplex(trial).epsilon(1e-3f));
		}
	}
}


TEST_CASE("Rotate tau series identity", "[Jacobi functions]") {
	const auto z = 0.7 + 0.3i;
	const std::array taus = {
		0.73 + 1.49i,
		-3.17 + 0.49i,
		0.77i,
		-0.11 + 1.03i,
	};

	for (int variant = 1; variant <= 4; ++variant) {
		for (const auto& tau : taus) {
			const auto [newZ, newTau, factor, exponent, newVariant] = RotateTau(z, tau, variant);
			const auto trial = factor * ThetaSeries(newVariant, newZ, newTau, exponent);
			const auto control = ThetaSeries(variant, z, tau, 0.0i);
			UNSCOPED_INFO("variant=" << variant << " -> " << newVariant
									 << "  tau=" << tau << " -> " << newTau);
			REQUIRE(control == ApproxComplex(trial).epsilon(1e-3f));
		}
	}
}

TEST_CASE("Evaluate theta_1", "[Jacobi functions]") {
	using real_t = double;

	const std::complex<real_t> q = -0.1f + 0.3if;
	const auto tau = -i_v<real_t> * std::log(q) / pi_v<real_t>;
	const auto q2 = std::exp(i_v<real_t> * pi_v<real_t> * tau);
	REQUIRE(q == ApproxComplex(q2));

	const auto z = LinSpace<real_t, DOMAINLESS>(-2*pi_v<float>, 2*pi_v<float>, 400);
	std::vector<real_t> r;
	std::vector<real_t> i;
	for (auto& v : z) {
		const auto c = theta(2, std::complex<real_t>(v, 0.7f), tau, 0);
		r.push_back(real(c));
		i.push_back(imag(c));
	}
	REQUIRE(false);
}