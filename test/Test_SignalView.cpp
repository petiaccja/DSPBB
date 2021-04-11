#include <catch2/catch.hpp>
#include <complex>
#include <dspbb/Primitives/SignalView.hpp>

using namespace dspbb;
using namespace std::complex_literals;

TEST_CASE("Default construct", "[SignalView]") {
	SignalView<float, TIME_DOMAIN> span;
	REQUIRE(span.Empty());
	REQUIRE(span.Size() == 0);

	SignalView<std::complex<float>, TIME_DOMAIN> cspan;
	REQUIRE(cspan.Empty());
	REQUIRE(cspan.Size() == 0);
}


TEST_CASE("Whole span", "[SignalView]") {
	TimeSignalF signal = { 1, 2, 3, 4, 5, 6 };

	SignalView<float, TIME_DOMAIN> span{ signal };
	REQUIRE(span.Size() == signal.Size());
	REQUIRE(span[0] == 1);
	REQUIRE(span[5] == 6);
}


TEST_CASE("Partial span size", "[SignalView]") {
	TimeSignalF signal = { 1, 2, 3, 4, 5, 6 };

	SignalView<float, TIME_DOMAIN> span{ signal.begin() + 3, 2 };
	REQUIRE(span.Size() == 2);
	REQUIRE(span[0] == 4);
	REQUIRE(span[1] == 5);
}


TEST_CASE("Partial span iterators", "[SignalView]") {
	TimeSignalF signal = { 1, 2, 3, 4, 5, 6 };

	SignalView<float, TIME_DOMAIN> span{ signal.begin() + 2, signal.begin() + 4 };
	REQUIRE(span.Size() == 2);
	REQUIRE(span[0] == 3);
	REQUIRE(span[1] == 4);
}


TEST_CASE("Data pointer", "[SignalView]") {
	TimeSignalF signal = { 1, 2, 3, 4, 5, 6 };

	SignalView<float, TIME_DOMAIN> span{ signal.begin() + 2, signal.begin() + 4 };
	REQUIRE(*span.Data() == 3);
}


TEST_CASE("Real/Imag pointer", "[SignalView]") {
	using namespace std::complex_literals;
	TimeSignalCF signal = { 1.f + 2.if,
							2.f + 3.if,
							3.f + 6.if,
							4.f + 7.if,
							5.f + 8.if,
							6.f + 9.if };

	SignalView<std::complex<float>, TIME_DOMAIN> span{ signal.begin() + 2, signal.begin() + 4 };
	REQUIRE(span.Data()->real() == 3);
	REQUIRE(span.Data()->imag() == 6);
}


TEST_CASE("Constant span", "[SignalView]") {
	TimeSignalF signal = { 1, 2, 3, 4, 5, 6 };

	SignalView<const float, TIME_DOMAIN> span{ signal.begin() + 2, signal.begin() + 4 };
	REQUIRE(*span.Data() == 3);

	TimeSignalCF csignal = { 1.f + 2.if,
							 2.f + 3.if,
							 3.f + 6.if,
							 4.f + 7.if,
							 5.f + 8.if,
							 6.f + 9.if };

	SignalView<const std::complex<float>, TIME_DOMAIN> cspan{ csignal.begin() + 2, csignal.begin() + 4 };
	REQUIRE(cspan.Data()->real() == 3);
}

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