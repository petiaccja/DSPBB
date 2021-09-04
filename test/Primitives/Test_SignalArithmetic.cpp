#include "../TestUtils.hpp"

#include <catch2/catch.hpp>
#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/Primitives/SignalArithmetic.hpp>
#include <dspbb/Primitives/SignalView.hpp>

using namespace dspbb;

static_assert(is_mutable_signal_v<Signal<float, TIME_DOMAIN>>, "");
static_assert(is_mutable_signal_v<Signal<float, TIME_DOMAIN>&>, "");
static_assert(is_mutable_signal_v<Signal<float, TIME_DOMAIN>&&>, "");
static_assert(!is_mutable_signal_v<const Signal<float, TIME_DOMAIN>>, "");
static_assert(!is_mutable_signal_v<const Signal<float, TIME_DOMAIN>&>, "");
static_assert(!is_mutable_signal_v<const Signal<float, TIME_DOMAIN>&&>, "");

static_assert(is_mutable_signal_v<SignalView<float, TIME_DOMAIN>>, "");
static_assert(is_mutable_signal_v<SignalView<float, TIME_DOMAIN>&>, "");
static_assert(is_mutable_signal_v<SignalView<float, TIME_DOMAIN>&&>, "");
static_assert(is_mutable_signal_v<const SignalView<float, TIME_DOMAIN>>, "");
static_assert(is_mutable_signal_v<const SignalView<float, TIME_DOMAIN>&>, "");
static_assert(is_mutable_signal_v<const SignalView<float, TIME_DOMAIN>&&>, "");
static_assert(!is_mutable_signal_v<SignalView<const float, TIME_DOMAIN>>, "");
static_assert(!is_mutable_signal_v<SignalView<const float, TIME_DOMAIN>&>, "");
static_assert(!is_mutable_signal_v<SignalView<const float, TIME_DOMAIN>&&>, "");
static_assert(!is_mutable_signal_v<const SignalView<const float, TIME_DOMAIN>>, "");
static_assert(!is_mutable_signal_v<const SignalView<const float, TIME_DOMAIN>&>, "");
static_assert(!is_mutable_signal_v<const SignalView<const float, TIME_DOMAIN>&&>, "");



//------------------------------------------------------------------------------
// Operators
//------------------------------------------------------------------------------
/*
const TimeSignal<float> s1 = { 1, 2, 3 };
const TimeSignal<float> s2 = { 7, 4, 5 };
const float c1 = 5;

#define TEST_CASE_VIEW_OPERATOR(NAME, OPERATOR)                             \
	TEST_CASE("SignalView operator " NAME, "[SignalView]") { \
		const auto expected = s1 OPERATOR s2;                               \
		const auto v1 = s1 OPERATOR AsConstView(s2);                        \
		const auto v2 = AsConstView(s1) OPERATOR s2;                        \
		const auto v3 = AsConstView(s1) OPERATOR AsConstView(s2);           \
		for (size_t i = 0; i < s1.Size(); ++i) {                            \
			REQUIRE(expected[i] == Approx(v1[i]));                          \
			REQUIRE(expected[i] == Approx(v2[i]));                          \
			REQUIRE(expected[i] == Approx(v3[i]));                          \
		}                                                                   \
	}

#define TEST_CASE_VIEW_OPERATOR_SCALAR(NAME, OPERATOR)                                    \
	TEST_CASE("SignalView operator scalar " NAME, "[SignalView]") { \
		const auto expected1 = c1 OPERATOR s1;                                     \
		const auto v1 = c1 OPERATOR AsConstView(s1);                               \
		const auto expected2 = s1 OPERATOR c1;                                     \
		const auto v2 = AsConstView(s1) OPERATOR c1;                               \
		for (size_t i = 0; i < s1.Size(); ++i) {                                   \
			REQUIRE(expected1[i] == Approx(v1[i]));                                \
			REQUIRE(expected2[i] == Approx(v2[i]));                                \
		}                                                                          \
	}


TEST_CASE_VIEW_OPERATOR("multiply", *)
TEST_CASE_VIEW_OPERATOR("add", +)
TEST_CASE_VIEW_OPERATOR("subtract", -)
TEST_CASE_VIEW_OPERATOR("divide", /)

TEST_CASE_VIEW_OPERATOR_SCALAR("multiply", *)
TEST_CASE_VIEW_OPERATOR_SCALAR("add", +)
TEST_CASE_VIEW_OPERATOR_SCALAR("subtract", -)
TEST_CASE_VIEW_OPERATOR_SCALAR("divide", /)
*/



#define NAME "multiply"
#define OPERATOR *

#define TEST_SIGNAL_BINARY_OPERATOR(NAME, OPERATOR)                                                        \
	TEMPLATE_PRODUCT_TEST_CASE("Signal binary", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) { \
		constexpr size_t count = 1;                                                                        \
                                                                                                           \
		using TestType0 = std::tuple_element_t<0, TestType>;                                               \
		using TestType1 = std::tuple_element_t<1, TestType>;                                               \
                                                                                                           \
		const auto a = RandomPositiveSignal<TestType0>(count);                                             \
		const auto b = RandomPositiveSignal<TestType1>(count);                                             \
		const auto r = a OPERATOR b;                                                                       \
                                                                                                           \
		size_t numEqual = 0;                                                                               \
		for (size_t i = 0; i < count; ++i) {                                                               \
			numEqual += size_t(r[i] == ApproxComplex(a[i] OPERATOR b[i]));                                 \
		}                                                                                                  \
		REQUIRE(numEqual == count);                                                                        \
	}


template <class T>
constexpr T a0 = static_cast<T>(3.5f);
template <class T>
constexpr T a1 = static_cast<T>(2.9f);
template <class T>
constexpr T b0 = static_cast<T>(9.3f);
template <class T>
constexpr T b1 = static_cast<T>(2.5f);
template <class T>
constexpr T bs = static_cast<T>(2.63f);
template <class T>
constexpr T as = static_cast<T>(1.75f);

TEMPLATE_PRODUCT_TEST_CASE("Multiply", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	TimeSignal<TestType1> b = { b0<TestType1>, b1<TestType1> };
	auto r = a * b;
	REQUIRE(r[0] == ApproxComplex(a0<TestType0> * b0<TestType1>));
	REQUIRE(r[1] == ApproxComplex(a1<TestType0> * b1<TestType1>));
}

TEMPLATE_PRODUCT_TEST_CASE("Divide", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	TimeSignal<TestType1> b = { b0<TestType1>, b1<TestType1> };
	auto r = a / b;
	;
	REQUIRE(r[1] == ApproxComplex(a1<TestType0> / b1<TestType1>));
	REQUIRE(r[0] == ApproxComplex(a0<TestType0> / b0<TestType1>));
}

TEMPLATE_PRODUCT_TEST_CASE("Add", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	TimeSignal<TestType1> b = { b0<TestType1>, b1<TestType1> };
	auto r = a + b;
	auto rv = AsView(a) * b;
	REQUIRE(r[0] == ApproxComplex(a0<TestType0> + b0<TestType1>));
	REQUIRE(r[1] == ApproxComplex(a1<TestType0> + b1<TestType1>));
}

TEMPLATE_PRODUCT_TEST_CASE("Subtract", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	TimeSignal<TestType1> b = { b0<TestType1>, b1<TestType1> };
	auto r = a - b;
	REQUIRE(r[0] == ApproxComplex(a0<TestType0> - b0<TestType1>));
	REQUIRE(r[1] == ApproxComplex(a1<TestType0> - b1<TestType1>));
}

TEMPLATE_PRODUCT_TEST_CASE("Multiply mix", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	TimeSignal<TestType1> b = { b0<TestType1>, b1<TestType1> };
	auto r1 = AsView(a) * b;
	auto r2 = AsView(a) * AsView(b);
	auto r3 = a * AsView(b);
	auto r4 = AsView(a) * AsConstView(b);
	REQUIRE(r1[0] == r2[0]);
	REQUIRE(r1[0] == r3[0]);
	REQUIRE(r1[0] == r4[0]);
	REQUIRE(r1[1] == r2[1]);
	REQUIRE(r1[1] == r3[1]);
	REQUIRE(r1[1] == r4[1]);
}

TEMPLATE_PRODUCT_TEST_CASE("Compound multiply", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_REAL) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	TimeSignal<TestType1> b = { b0<TestType1>, b1<TestType1> };
	a *= b;
	REQUIRE(a[0] == ApproxComplex(a0<TestType0> * b0<TestType1>));
	REQUIRE(a[1] == ApproxComplex(a1<TestType0> * b1<TestType1>));
}

TEMPLATE_PRODUCT_TEST_CASE("Compound divide", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_REAL) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	TimeSignal<TestType1> b = { b0<TestType1>, b1<TestType1> };
	a /= b;
	REQUIRE(a[0] == ApproxComplex(a0<TestType0> / b0<TestType1>));
	REQUIRE(a[1] == ApproxComplex(a1<TestType0> / b1<TestType1>));
}

TEMPLATE_PRODUCT_TEST_CASE("Compound add", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_REAL) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	TimeSignal<TestType1> b = { b0<TestType1>, b1<TestType1> };
	a += b;
	REQUIRE(a[0] == ApproxComplex(a0<TestType0> + b0<TestType1>));
	REQUIRE(a[1] == ApproxComplex(a1<TestType0> + b1<TestType1>));
}

TEMPLATE_PRODUCT_TEST_CASE("Compound subtract", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_REAL) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	TimeSignal<TestType1> b = { b0<TestType1>, b1<TestType1> };
	a -= b;
	REQUIRE(a[0] == ApproxComplex(a0<TestType0> - b0<TestType1>));
	REQUIRE(a[1] == ApproxComplex(a1<TestType0> - b1<TestType1>));
}


TEMPLATE_PRODUCT_TEST_CASE("Scalar multiply", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	auto r = a * bs<TestType1>;
	REQUIRE(r[0] == ApproxComplex(a0<TestType0> * bs<TestType1>));
	REQUIRE(r[1] == ApproxComplex(a1<TestType0> * bs<TestType1>));
}

TEMPLATE_PRODUCT_TEST_CASE("Scalar divide", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	auto r = a / bs<TestType1>;
	REQUIRE(r[0] == ApproxComplex(a0<TestType0> / bs<TestType1>));
	REQUIRE(r[1] == ApproxComplex(a1<TestType0> / bs<TestType1>));
}

TEMPLATE_PRODUCT_TEST_CASE("Scalar add", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	auto r = a + TestType1(2);
	REQUIRE(r[0] == ApproxComplex(a0<TestType0> + TestType1(2)));
	REQUIRE(r[1] == ApproxComplex(a1<TestType0> + TestType1(2)));
}

TEMPLATE_PRODUCT_TEST_CASE("Scalar subs<TestType1>tract", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	auto r = a - bs<TestType1>;
	REQUIRE(r[0] == ApproxComplex(a0<TestType0> - bs<TestType1>));
	REQUIRE(r[1] == ApproxComplex(a1<TestType0> - bs<TestType1>));
}


TEMPLATE_PRODUCT_TEST_CASE("Scalar reverse multiply", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	auto r = bs<TestType1> * a;
	REQUIRE(r[0] == ApproxComplex(bs<TestType1> * a0<TestType0>));
	REQUIRE(r[1] == ApproxComplex(bs<TestType1> * a1<TestType0>));
}

TEMPLATE_PRODUCT_TEST_CASE("Scalar reverse divide", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	auto r = bs<TestType1> / a;
	REQUIRE(r[0] == ApproxComplex(bs<TestType1> / a0<TestType0>));
	REQUIRE(r[1] == ApproxComplex(bs<TestType1> / a1<TestType0>));
}

TEMPLATE_PRODUCT_TEST_CASE("Scalar reverse add", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	auto r = TestType1(2) + a;
	REQUIRE(r[0] == ApproxComplex(TestType1(2) + a0<TestType0>));
	REQUIRE(r[1] == ApproxComplex(TestType1(2) + a1<TestType0>));
}

TEMPLATE_PRODUCT_TEST_CASE("Scalar reverse subs<TestType1>tract", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	auto r = bs<TestType1> - a;
	REQUIRE(r[0] == ApproxComplex(bs<TestType1> - a0<TestType0>));
	REQUIRE(r[1] == ApproxComplex(bs<TestType1> - a1<TestType0>));
}


TEMPLATE_PRODUCT_TEST_CASE("Scalar compound multiply", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_REAL) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	a *= bs<TestType1>;
	REQUIRE(a[0] == ApproxComplex(a0<TestType0> * bs<TestType1>));
	REQUIRE(a[1] == ApproxComplex(a1<TestType0> * bs<TestType1>));
}

TEMPLATE_PRODUCT_TEST_CASE("Scalar compound divide", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_REAL) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	a /= bs<TestType1>;
	REQUIRE(a[0] == ApproxComplex(a0<TestType0> / bs<TestType1>));
	REQUIRE(a[1] == ApproxComplex(a1<TestType0> / bs<TestType1>));
}

TEMPLATE_PRODUCT_TEST_CASE("Scalar compound add", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_REAL) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	a += TestType1(2);
	REQUIRE(a[0] == ApproxComplex(a0<TestType0> + TestType1(2)));
	REQUIRE(a[1] == ApproxComplex(a1<TestType0> + TestType1(2)));
}

TEMPLATE_PRODUCT_TEST_CASE("Scalar compound subs<TestType1>tract", "[Signal Arithmetic]", std::tuple, TYPES_BINARY_REAL) {
	using TestType0 = std::tuple_element_t<0, TestType>;
	using TestType1 = std::tuple_element_t<1, TestType>;

	TimeSignal<TestType0> a = { a0<TestType0>, a1<TestType0> };
	a -= bs<TestType1>;
	REQUIRE(a[0] == ApproxComplex(a0<TestType0> - bs<TestType1>));
	REQUIRE(a[1] == ApproxComplex(a1<TestType0> - bs<TestType1>));
}