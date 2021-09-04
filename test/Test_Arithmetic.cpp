#include "TestUtils.hpp"

#include <catch2/catch.hpp>
#include <dspbb/Math/Arithmetic.hpp>

using namespace dspbb;


TEMPLATE_PRODUCT_TEST_CASE("Multiply array", "[Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;
	using ResultType = BinaryProd<TestType0, TestType1>;

	const std::array<TestType0, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const std::array<TestType1, 9> b = { 9, 8, 7, 6, 5, 4, 3, 2, 1 };
	const std::array<ResultType, 9> e = { 9, 16, 21, 24, 25, 24, 21, 16, 9 };
	std::array<ResultType, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	Multiply(r.data(), a.data(), b.data(), r.size());
	REQUIRE(r == e);
}

TEMPLATE_PRODUCT_TEST_CASE("Multiply scalar", "[Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;
	using ResultType = BinaryProd<TestType0, TestType1>;
	
	const std::array<TestType0, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const TestType1 b = 2;
	const std::array<ResultType, 9> e = { 2, 4, 6, 8, 10, 12, 14, 16, 18 };
	std::array<ResultType, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	SECTION("array x scalar") {
		Multiply(r.data(), a.data(), b, r.size());
		REQUIRE(r == e);
	}	

	SECTION("scalar x array") {
		Multiply(r.data(), b, a.data(), r.size());
		REQUIRE(r == e);
	}
}

TEMPLATE_PRODUCT_TEST_CASE("Divide array", "[Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;
	using ResultType = BinaryQuot<TestType0, TestType1>;

	const std::array<TestType0, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const std::array<TestType1, 9> b = { 0.5f, 0.25f, 0.5f, 0.25f, 0.5f, 0.25f, 0.5f, 0.25f, 1.f / 3.f };
	const std::array<ResultType, 9> e = { 2, 8, 6, 16, 10, 24, 14, 32, 27 };
	std::array<ResultType, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	
	Divide(r.data(), a.data(), b.data(), r.size());
	for (size_t i = 0; i < 9; ++i) {
		REQUIRE(r[i] == ApproxComplex(e[i]));
	}
}

TEMPLATE_PRODUCT_TEST_CASE("Divide scalar", "[Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;
	using ResultType = BinaryQuot<TestType0, TestType1>;

	const std::array<TestType0, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const TestType1 b = 0.5f;
	const std::array<ResultType, 9> e = { 2, 4, 6, 8, 10, 12, 14, 16, 18 };
	std::array<ResultType, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	SECTION("array x scalar") {
		Divide(r.data(), a.data(), b, r.size());
		for (size_t i = 0; i < 9; ++i) {
			REQUIRE(r[i] == ApproxComplex(e[i]));
		}
	}

	SECTION("scalar x array") {
		Divide(r.data(), b, a.data(), r.size());
		Divide(r.data(), ResultType(1.0f), r.data(), r.size());
		for (size_t i = 0; i < 9; ++i) {
			REQUIRE(r[i] == ApproxComplex(e[i]));
		}
	}
}

TEMPLATE_PRODUCT_TEST_CASE("Add array", "[Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;
	using ResultType = BinarySum<TestType0, TestType1>;

	const std::array<TestType0, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const std::array<TestType1, 9> b = { 2, 3, 4, 5, 6, 7, 8, 9, 10 };
	const std::array<ResultType, 9> e = { 3, 5, 7, 9, 11, 13, 15, 17, 19 };
	std::array<ResultType, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	
	Add(r.data(), a.data(), b.data(), r.size());
	REQUIRE(r == e);
}

TEMPLATE_PRODUCT_TEST_CASE("Add scalar", "[Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;
	using ResultType = BinarySum<TestType0, TestType1>;

	const std::array<TestType0, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const TestType1 b = 2.0f;
	const std::array<ResultType, 9> e = { 3, 4, 5, 6, 7, 8, 9, 10, 11 };
	std::array<ResultType, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	SECTION("array x scalar") {
		Add(r.data(), a.data(), b, r.size());
		REQUIRE(r == e);
	}	

	SECTION("scalar x array") {
		Add(r.data(), b, a.data(), r.size());
		REQUIRE(r == e);
	}
}

TEMPLATE_PRODUCT_TEST_CASE("Sub array", "[Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;
	using ResultType = BinaryDiff<TestType0, TestType1>;

	const std::array<TestType0, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const std::array<TestType1, 9> b = { 9, 8, 7, 6, 5, 4, 3, 2, 1 };
	const std::array<ResultType, 9> e = { -8, -6, -4, -2, 0, 2, 4, 6, 8 };
	std::array<ResultType, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	
	Subtract(r.data(), a.data(), b.data(), r.size());
	REQUIRE(r == e);
}

TEMPLATE_PRODUCT_TEST_CASE("Sub scalar", "[Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;
	using ResultType = BinaryDiff<TestType0, TestType1>;

	const std::array<TestType0, 9> a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	const TestType1 b = -2.0f;
	const std::array<ResultType, 9> e = { 3, 4, 5, 6, 7, 8, 9, 10, 11 };
	std::array<ResultType, 9> r = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	SECTION("array x scalar") {
		Subtract(r.data(), a.data(), b, r.size());
		REQUIRE(r == e);
	}
	SECTION("scalar x array") {
		Subtract(r.data(), b, a.data(), r.size());
		Multiply(r.data(), ResultType(-1.f), r.data(), r.size());
		REQUIRE(r == e);
	}
}
