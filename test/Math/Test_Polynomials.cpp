#include "../TestUtils.hpp"

#include <dspbb/Math/Polynomials.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>


using namespace dspbb;
using namespace std::complex_literals;
using Catch::Approx;


TEST_CASE("empty polynomial", "[Polynomials]") {
	Polynomial<float> poly;
	REQUIRE(poly.size() == 0);
}

TEST_CASE("Non-empty polynomial", "[Polynomials]") {
	Polynomial<float> poly = { 1, 2, 3 };
	REQUIRE(poly.size() == 3);
	REQUIRE(poly.coefficients()[0] == Approx(1));
	REQUIRE(poly.coefficients()[1] == Approx(2));
	REQUIRE(poly.coefficients()[2] == Approx(3));
}

TEST_CASE("Polynomial resize", "[Polynomials]") {
	Polynomial<float> poly;
	poly.resize(3, 1.0f);
	REQUIRE(poly.coefficients().size() == 3);
	REQUIRE(std::all_of(poly.coefficients().begin(), poly.coefficients().end(), [](auto v) { return v == Approx(1); }));
}

TEST_CASE("Polynomial real evaluate", "[Polynomials]") {
	Polynomial<float> poly{ 3, 4, -2 }; // 3 + 4x -2x^2
	const float x = 1.5f;
	REQUIRE(poly(x) == Approx(3.f + 4.f * x - 2.f * x * x));
}

TEST_CASE("Polynomial complex evaluate", "[Polynomials]") {
	Polynomial<float> poly{ 3, 4, -2 }; // 3 + 4x -2x^2
	const std::complex<float> x = 1.5f + 2.3if;
	REQUIRE(poly(x) == ApproxComplex(3.f + 4.f * x - 2.f * x * x));
}



TEST_CASE("empty factored polynomial", "[Polynomials]") {
	FactoredPolynomial<float> poly;
	REQUIRE(poly.num_roots() == 0);
}

TEST_CASE("Factored polynomial missing pairs", "[Polynomials]") {
	REQUIRE_THROWS(FactoredPolynomial<float>{ 1.0f, 2.0f + 1if });
}

TEST_CASE("Factored polynomial no pairs", "[Polynomials]") {
	REQUIRE_THROWS(FactoredPolynomial<float>{ 1.0f, 2.0f + 1if, 2.0f - 0.9if });
}

TEST_CASE("Non-empty factored polynomial", "[Polynomials]") {
	FactoredPolynomial<float> poly = { 1.0f, 3.0f, 2.0f + 1if, 2.0f - 1if };
	REQUIRE(poly.num_roots() == 4);
	REQUIRE(poly.num_real_roots() == 2);
	REQUIRE(poly.num_complex_pairs() == 1);
	REQUIRE(poly.num_complex_roots() == 2);

	REQUIRE(poly.real_roots().size() == 2);
	REQUIRE(poly.complex_pairs().size() == 1);

	REQUIRE(poly.real_roots()[0] == Approx(1));
	REQUIRE(poly.real_roots()[1] == Approx(3));
	REQUIRE(poly.complex_pairs()[0] == ApproxComplex(2.0f + 1if));
}

TEST_CASE("Factored polynomial resize initial", "[Polynomials]") {
	FactoredPolynomial<float> poly;
	poly.resize(2, 1, 1.0f, 2.f + 1.0if);
	REQUIRE(poly.real_roots().size() == 2);
	REQUIRE(poly.complex_pairs().size() == 1);
	REQUIRE(std::all_of(poly.real_roots().begin(), poly.real_roots().end(), [](auto v) { return v == 1.0f; }));
	REQUIRE(std::all_of(poly.complex_pairs().begin(), poly.complex_pairs().end(), [](auto v) { return v == 2.f + 1.0if; }));
}

TEST_CASE("Factored polynomial resize shrink/grow", "[Polynomials]") {
	FactoredPolynomial<float> poly;
	const float r1 = 1.0f;
	const float r2 = 2.0f;
	const std::complex<float> c1 = 10.f + 10.if;
	const std::complex<float> c2 = 20.f + 20.if;
	poly.resize(6, 4, r1, c1);
	poly.resize(4, 6, r2, c2);
	REQUIRE(poly.real_roots().size() == 4);
	REQUIRE(poly.complex_pairs().size() == 6);
	REQUIRE(std::all_of(poly.real_roots().begin(), poly.real_roots().end(), [&](auto v) { return v == r1; }));
	REQUIRE(std::all_of(poly.complex_pairs().begin(), poly.complex_pairs().begin() + 4, [&](auto v) { return v == c1; }));
	REQUIRE(std::all_of(poly.complex_pairs().begin() + 4, poly.complex_pairs().begin() + 6, [&](auto v) { return v == c2; }));
}

TEST_CASE("Factored polynomial resize grow/shrink", "[Polynomials]") {
	FactoredPolynomial<float> poly;
	const float r1 = 1.0f;
	const float r2 = 2.0f;
	const std::complex<float> c1 = 10.f + 11.if;
	const std::complex<float> c2 = 20.f + 22.if;
	poly.resize(4, 6, r1, c1);
	poly.resize(6, 4, r2, c2);
	REQUIRE(poly.real_roots().size() == 6);
	REQUIRE(poly.complex_pairs().size() == 4);
	REQUIRE(std::all_of(poly.real_roots().begin(), poly.real_roots().begin() + 4, [&](auto v) { return v == r1; }));
	REQUIRE(std::all_of(poly.real_roots().begin() + 4, poly.real_roots().begin() + 6, [&](auto v) { return v == r2; }));
	REQUIRE(std::all_of(poly.complex_pairs().begin(), poly.complex_pairs().end(), [&](auto v) { return v == c1; }));
}


TEST_CASE("Factored polynomial regroup shrink/grow", "[Polynomials]") {
	FactoredPolynomial<float> poly;
	const float r1 = 1.0f;
	const float r2 = 2.0f;
	const std::complex<float> c1 = 10.f + 10.if;
	const std::complex<float> c2 = 20.f + 20.if;
	poly.resize(6, 4, r1, c1);
	poly.regroup(4, r2, c2);
	REQUIRE(poly.real_roots().size() == 4);
	REQUIRE(poly.complex_pairs().size() == 5);
	REQUIRE(std::all_of(poly.real_roots().begin(), poly.real_roots().end(), [&](auto v) { return v == r1; }));
	REQUIRE(std::all_of(poly.complex_pairs().begin(), poly.complex_pairs().begin() + 4, [&](auto v) { return v == c1; }));
	REQUIRE(std::all_of(poly.complex_pairs().begin() + 4, poly.complex_pairs().begin() + 5, [&](auto v) { return v == c2; }));
}

TEST_CASE("Factored polynomial regroup grow/shrink", "[Polynomials]") {
	FactoredPolynomial<float> poly;
	const float r1 = 1.0f;
	const float r2 = 2.0f;
	const std::complex<float> c1 = 10.f + 11.if;
	const std::complex<float> c2 = 20.f + 22.if;
	poly.resize(4, 6, r1, c1);
	poly.regroup(6, r2, c2);
	REQUIRE(poly.real_roots().size() == 6);
	REQUIRE(poly.complex_pairs().size() == 5);
	REQUIRE(std::all_of(poly.real_roots().begin(), poly.real_roots().begin() + 4, [&](auto v) { return v == r1; }));
	REQUIRE(std::all_of(poly.real_roots().begin() + 4, poly.real_roots().begin() + 6, [&](auto v) { return v == r2; }));
	REQUIRE(std::all_of(poly.complex_pairs().begin(), poly.complex_pairs().end(), [&](auto v) { return v == c1; }));
}

TEST_CASE("Factored polynomial regroup oversize", "[Polynomials]") {
	FactoredPolynomial<float> poly;
	poly.resize(1, 3);
	REQUIRE_NOTHROW(poly.regroup(7));
	REQUIRE_THROWS(poly.regroup(9));
}

TEST_CASE("Factored polynomial regroup no pair", "[Polynomials]") {
	FactoredPolynomial<float> poly;
	poly.resize(1, 3);
	REQUIRE_NOTHROW(poly.regroup(3));
	REQUIRE_THROWS(poly.regroup(0));
}

TEST_CASE("Factored polynomial real evaluate even", "[Polynomials]") {
	FactoredPolynomial<float> poly{ 1.0f, 3.0f, 2.0f + 1if, 2.0f - 1if }; // 15 - 32 x + 24 x^2 - 8 x^3 + x^4
	const float x = 1.5f;
	REQUIRE(poly(x) == Approx(15.f - 32.f * x + 24.f * x * x - 8.f * x * x * x + x * x * x * x));
}

TEST_CASE("Factored polynomial complex evaluate even", "[Polynomials]") {
	FactoredPolynomial<float> poly{ 1.0f, 3.0f, 2.0f + 1if, 2.0f - 1if }; // 15 - 32 x + 24 x^2 - 8 x^3 + x^4
	const std::complex<float> x = 1.5f + 2.3if;
	REQUIRE(poly(x) == ApproxComplex(15.f - 32.f * x + 24.f * x * x - 8.f * x * x * x + x * x * x * x));
}

TEST_CASE("Factored polynomial real evaluate odd", "[Polynomials]") {
	FactoredPolynomial<float> poly{ 3.0f, -4.0f + 3if, -4.0f - 3if }; // x^3 + 5 x^2 + x - 75
	const float x = 1.5f;
	REQUIRE(poly(x) == Approx(x * x * x + 5 * x * x + x - 75));
}

TEST_CASE("Factored polynomial complex evaluate odd", "[Polynomials]") {
	FactoredPolynomial<float> poly{ 3.0f, -4.0f + 3if, -4.0f - 3if }; // x^3 + 5 x^2 + x - 75
	const std::complex<float> x = 1.5f + 2.3if;
	REQUIRE(poly(x) == ApproxComplex(x * x * x + 5.0f * x * x + x - 75.0f));
}

TEST_CASE("Expand polynomials", "[Polynomials]") {
	const FactoredPolynomial<float> factored{ 1.0f, 3.0f, 2.0f + 1if, 2.0f - 1if };
	const Polynomial<float> expanded = ExpandPolynomial(factored);
	REQUIRE(expanded.coefficients()[0] == Approx(15.0f));
	REQUIRE(expanded.coefficients()[1] == Approx(-32.0f));
	REQUIRE(expanded.coefficients()[2] == Approx(24.0f));
	REQUIRE(expanded.coefficients()[3] == Approx(-8.0f));
	REQUIRE(expanded.coefficients()[4] == Approx(1.0f));
}