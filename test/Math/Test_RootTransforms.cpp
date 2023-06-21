#include "../TestUtils.hpp"

#include <dspbb/Math/RootTransforms.hpp>

#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>


using namespace dspbb;
using namespace std::complex_literals;
using Catch::Approx;


const FactoredPolynomial<float> rootsDouble = {
	1.0f, 2.0f, 3.0f + 2.0if, 3.0f - 2.0if
};

const auto double1st = [](const auto& root) {
	using T = remove_complex_t<std::decay_t<decltype(root)>>;
	const std::complex<T> z = 2.0f * root;
	return std::array{ z };
};

const FactoredPolynomial<float> rootsSpin = {
	0.0f, 2.0f, 3.0f + 2.0if, 3.0f - 2.0if
};

const auto spin2nd = [](const auto& root) {
	using T = remove_complex_t<std::decay_t<decltype(root)>>;
	const std::complex<T> z1 = root * std::complex<T>(T(0), T(1));
	const std::complex<T> z2 = root * std::complex<T>(T(0), -T(1));
	return std::array{ z1, z2 };
};

const auto faulty2nd = [](const auto&) {
	return std::array{ -1.0f + 1.0if, -1.0f + 1.0if }; // Not a conjugate pair.
};

TEST_CASE("First order default root count", "[Root transforms]") {
	const auto transformedRoots = TransformRoots<float, 1>(rootsDouble, double1st);
	REQUIRE(transformedRoots.num_real_roots() == rootsDouble.num_real_roots());
	REQUIRE(transformedRoots.num_complex_pairs() == rootsDouble.num_complex_pairs());
	REQUIRE(transformedRoots.real_roots()[0] == Approx(2.0f));
	REQUIRE(transformedRoots.real_roots()[1] == Approx(4.0f));
	REQUIRE(transformedRoots.complex_pairs()[0] == ApproxComplex(6.0f + 4.0if));
}

TEST_CASE("First order padding", "[Root transforms]") {
	const auto transformedRoots = TransformRoots(rootsDouble, double1st, 8, std::array<std::complex<float>, 1>{ -1.0f });
	REQUIRE(transformedRoots.num_roots() == 8);
	REQUIRE(transformedRoots.num_complex_pairs() == rootsDouble.num_complex_pairs());
	REQUIRE(transformedRoots.num_real_roots() == 6);
	REQUIRE(transformedRoots.real_roots()[0] == Approx(2.0f));
	REQUIRE(transformedRoots.real_roots()[1] == Approx(4.0f));
	REQUIRE(transformedRoots.real_roots()[2] == Approx(-1.0f));
	REQUIRE(transformedRoots.real_roots()[3] == Approx(-1.0f));
	REQUIRE(transformedRoots.real_roots()[4] == Approx(-1.0f));
	REQUIRE(transformedRoots.real_roots()[5] == Approx(-1.0f));
	REQUIRE(transformedRoots.complex_pairs()[0] == ApproxComplex(6.0f + 4.0if));
}

TEST_CASE("second order default root count", "[Root transforms]") {
	const auto transformedRoots = TransformRoots<float, 2>(rootsSpin, spin2nd);
	REQUIRE(transformedRoots.num_real_roots() == 2);
	REQUIRE(transformedRoots.num_complex_roots() == 6);
	REQUIRE(transformedRoots.real_roots()[0] == Approx(0.0f));
	REQUIRE(transformedRoots.real_roots()[1] == Approx(0.0f));
	REQUIRE(transformedRoots.complex_pairs()[0] == ApproxComplex(2.0if));
	REQUIRE(transformedRoots.complex_pairs()[1] == ApproxComplex(3.0if - 2.0f));
	REQUIRE(transformedRoots.complex_pairs()[2] == ApproxComplex(-3.0if + 2.0f));
}

TEST_CASE("Second order padding with real", "[Root transforms]") {
	const auto transformedRoots = TransformRoots(rootsSpin, spin2nd, 8, std::array<std::complex<float>, 2>{ -1.0f, -2.0f });
	REQUIRE(transformedRoots.num_roots() == 16);
	REQUIRE(transformedRoots.num_real_roots() == 10);
	REQUIRE(transformedRoots.num_complex_roots() == 6);
	REQUIRE(transformedRoots.real_roots()[0] == Approx(0.0f));
	REQUIRE(transformedRoots.real_roots()[0] == Approx(0.0f));
	REQUIRE(transformedRoots.real_roots()[2] == Approx(-1.0f));
	REQUIRE(transformedRoots.real_roots()[3] == Approx(-2.0f));
	REQUIRE(transformedRoots.real_roots()[4] == Approx(-1.0f));
	REQUIRE(transformedRoots.real_roots()[5] == Approx(-2.0f));
	REQUIRE(transformedRoots.real_roots()[6] == Approx(-1.0f));
	REQUIRE(transformedRoots.real_roots()[7] == Approx(-2.0f));
	REQUIRE(transformedRoots.real_roots()[8] == Approx(-1.0f));
	REQUIRE(transformedRoots.real_roots()[9] == Approx(-2.0f));
	REQUIRE(transformedRoots.complex_pairs()[0] == ApproxComplex(2.0if));
	REQUIRE(transformedRoots.complex_pairs()[1] == ApproxComplex(3.0if - 2.0f));
	REQUIRE(transformedRoots.complex_pairs()[2] == ApproxComplex(-3.0if + 2.0f));
}

TEST_CASE("Second order padding with complex pair", "[Root transforms]") {
	const auto transformedRoots = TransformRoots(rootsSpin, spin2nd, 8, std::array<std::complex<float>, 2>{ -1.0f + 1.0if, -1.0f - 1.0if });
	REQUIRE(transformedRoots.num_roots() == 16);
	REQUIRE(transformedRoots.num_real_roots() == 2);
	REQUIRE(transformedRoots.num_complex_roots() == 14);
	REQUIRE(transformedRoots.real_roots()[0] == Approx(0.0f));
	REQUIRE(transformedRoots.real_roots()[0] == Approx(0.0f));
	REQUIRE(transformedRoots.complex_pairs()[0] == ApproxComplex(2.0if));
	REQUIRE(transformedRoots.complex_pairs()[1] == ApproxComplex(3.0if - 2.0f));
	REQUIRE(transformedRoots.complex_pairs()[2] == ApproxComplex(-3.0if + 2.0f));
	REQUIRE(transformedRoots.complex_pairs()[3] == ApproxComplex(-1.0f + 1.0if));
	REQUIRE(transformedRoots.complex_pairs()[4] == ApproxComplex(-1.0f + 1.0if));
	REQUIRE(transformedRoots.complex_pairs()[5] == ApproxComplex(-1.0f + 1.0if));
	REQUIRE(transformedRoots.complex_pairs()[6] == ApproxComplex(-1.0f + 1.0if));
}

TEST_CASE("Second order padding not complex pair", "[Root transforms]") {
	REQUIRE_THROWS(TransformRoots(rootsSpin, spin2nd, 8, std::array<std::complex<float>, 2>{ -1.0f + 1.0if, -1.0f + 1.0if }));
}

TEST_CASE("Second order transform not complex pair", "[Root transforms]") {
	REQUIRE_THROWS(TransformRoots<float, 2>(rootsSpin, faulty2nd));
}

TEST_CASE("Requesting too few roots", "[Root transforms]") {
	REQUIRE_THROWS(TransformRoots<float, 2>(rootsSpin, faulty2nd, 1));
}