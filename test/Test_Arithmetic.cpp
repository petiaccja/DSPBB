#include <catch2/catch.hpp>
#include <dspbb/Math/Arithmetic.hpp>

using namespace dspbb;

//------------------------------------------------------------------------------
// Test for array behaviour and compilation problems.
//------------------------------------------------------------------------------

TEST_CASE("Multiply float", "[Arithmetic]") {
	const std::array<float, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const std::array<float, 9> b = { 9, 8, 7, 6, 5, 4, 3, 2, 1 };
	const std::array<float, 9> e = { 9, 16, 21, 24, 25, 24, 21, 16, 9 };
	std::array<float, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	Multiply(r.data(), a.data(), b.data(), r.size());
	REQUIRE(r == e);
}


TEST_CASE("Multiply int32", "[Arithmetic]") {
	const std::array<int32_t, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const std::array<int32_t, 9> b = { 9, 8, 7, 6, 5, 4, 3, 2, 1 };
	const std::array<int32_t, 9> e = { 9, 16, 21, 24, 25, 24, 21, 16, 9 };
	std::array<int32_t, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	Multiply(r.data(), a.data(), b.data(), r.size());
	REQUIRE(r == e);
}


TEST_CASE("Multiply float x scalar", "[Arithmetic]") {
	const std::array<float, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const float b = 2;
	const std::array<float, 9> e = { 2, 4, 6, 8, 10, 12, 14, 16, 18 };
	std::array<float, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	Multiply(r.data(), a.data(), b, r.size());
	REQUIRE(r == e);
	std::fill(r.begin(), r.end(), 0.f);
	Multiply(r.data(), b, a.data(), r.size());
	REQUIRE(r == e);
}


TEST_CASE("Multiply int32 x scalar", "[Arithmetic]") {
	const std::array<int32_t, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const int32_t b = 2;
	const std::array<int32_t, 9> e = { 2, 4, 6, 8, 10, 12, 14, 16, 18 };
	std::array<int32_t, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	Multiply(r.data(), a.data(), b, r.size());
	REQUIRE(r == e);
	std::fill(r.begin(), r.end(), 0);
	Multiply(r.data(), b, a.data(), r.size());
	REQUIRE(r == e);
}


TEST_CASE("Multiply complex<float>", "[Arithmetic]") {
	const std::array<std::complex<float>, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const std::array<std::complex<float>, 9> b = { 9, 8, 7, 6, 5, 4, 3, 2, 1 };
	const std::array<std::complex<float>, 9> e = { 9, 16, 21, 24, 25, 24, 21, 16, 9 };
	std::array<std::complex<float>, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	Multiply(r.data(), a.data(), b.data(), r.size());
	REQUIRE(r == e);
}


TEST_CASE("Multiply complex<float> x float", "[Arithmetic]") {
	const std::array<std::complex<float>, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const float b = 2.f;
	const std::array<std::complex<float>, 9> e = { 2, 4, 6, 8, 10, 12, 14, 16, 18 };
	std::array<std::complex<float>, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	Multiply(r.data(), a.data(), b, r.size());
	REQUIRE(r == e);
}


//------------------------------------------------------------------------------
// Test if the operators are doing the correct thing.
//------------------------------------------------------------------------------

TEST_CASE("Divide float", "[Arithmetic]") {
	const std::array<float, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const std::array<float, 9> b = { 0.5f, 0.25f, 0.5f, 0.25f, 0.5f, 0.25f, 0.5f, 0.25f, 1.f / 3.f };
	const std::array<float, 9> e = { 2, 8, 6, 16, 10, 24, 14, 32, 27 };
	std::array<float, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	Divide(r.data(), a.data(), b.data(), r.size());
	REQUIRE(r == e);
}

TEST_CASE("Divide float x scalar", "[Arithmetic]") {
	const std::array<float, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const float b = 0.5f;
	const std::array<float, 9> e = { 2, 4, 6, 8, 10, 12, 14, 16, 18 };
	std::array<float, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	Divide(r.data(), a.data(), b, r.size());
	REQUIRE(r == e);
	std::fill(r.begin(), r.end(), 0.f);
	Divide(r.data(), b, a.data(), r.size());
	Divide(r.data(), 1.0f, r.data(), r.size());
	for (size_t i = 0; i < 9; ++i)
		REQUIRE(r[i] == Approx(e[i]));
}

TEST_CASE("Add float", "[Arithmetic]") {
	const std::array<float, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const std::array<float, 9> b = { 2, 3, 4, 5, 6, 7, 8, 9, 10 };
	const std::array<float, 9> e = { 3, 5, 7, 9, 11, 13, 15, 17, 19 };
	std::array<float, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	Add(r.data(), a.data(), b.data(), r.size());
	REQUIRE(r == e);
}

TEST_CASE("Add float x scalar", "[Arithmetic]") {
	const std::array<float, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const float b = 2.0f;
	const std::array<float, 9> e = { 3, 4, 5, 6, 7, 8, 9, 10, 11 };
	std::array<float, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	Add(r.data(), a.data(), b, r.size());
	REQUIRE(r == e);
	std::fill(r.begin(), r.end(), 0.f);
	Add(r.data(), b, a.data(), r.size());
	REQUIRE(r == e);
}


TEST_CASE("Sub float", "[Arithmetic]") {
	const std::array<float, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const std::array<float, 9> b = { 9, 8, 7, 6, 5, 4, 3, 2, 1 };
	const std::array<float, 9> e = { -8, -6, -4, -2, 0, 2, 4, 6, 8 };
	std::array<float, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	Subtract(r.data(), a.data(), b.data(), r.size());
	REQUIRE(r == e);
}

TEST_CASE("Sub float x scalar", "[Arithmetic]") {
	const std::array<float, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const float b = -2.0f;
	const std::array<float, 9> e = { 3, 4, 5, 6, 7, 8, 9, 10, 11 };
	std::array<float, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	Subtract(r.data(), a.data(), b, r.size());
	REQUIRE(r == e);
	std::fill(r.begin(), r.end(), 0.f);
	Subtract(r.data(), b, a.data(), r.size());
	Multiply(r.data(), -1.f, r.data(), r.size());
	REQUIRE(r == e);
}
