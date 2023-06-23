#include <dspbb/Primitives/SignalView.hpp>

#include <catch2/catch_test_macros.hpp>
#include <complex>

using namespace dspbb;
using namespace std::complex_literals;


TEST_CASE("Conversion construct", "[SignalView]") {
	BasicSignal<float, TIME_DOMAIN> smut(5);
	const BasicSignal<float, TIME_DOMAIN> s(5);

	BasicSignalView<float, TIME_DOMAIN> v1{ smut };
	BasicSignalView<const float, TIME_DOMAIN> v2{ s };
	BasicSignalView<const float, TIME_DOMAIN> v3{ smut };
	BasicSignalView<const float, TIME_DOMAIN> v4{ v1 };
	BasicSignalView<float, TIME_DOMAIN> v5{ smut.begin(), smut.end() };
	BasicSignalView<const float, TIME_DOMAIN> v6{ smut.begin(), smut.end() };
	BasicSignalView<const float, TIME_DOMAIN> v7{ s.begin(), s.end() };

	REQUIRE(v1.size() == smut.size());
	REQUIRE(v2.size() == s.size());
	REQUIRE(v3.size() == smut.size());
	REQUIRE(v4.size() == v1.size());
	REQUIRE(v5.size() == smut.size());
	REQUIRE(v6.size() == smut.size());
	REQUIRE(v6.size() == s.size());
}

TEST_CASE("View of", "[SignalView]") {
	BasicSignal<float, TIME_DOMAIN> smut(5);
	const BasicSignal<float, TIME_DOMAIN> s(5);

	BasicSignalView<float, TIME_DOMAIN> v1 = AsView(smut);
	BasicSignalView<const float, TIME_DOMAIN> v2 = AsConstView(smut);
	BasicSignalView<float, TIME_DOMAIN> v3 = AsView<TIME_DOMAIN>(smut.begin(), smut.end());
	BasicSignalView<const float, TIME_DOMAIN> v4 = AsConstView<TIME_DOMAIN>(smut.begin(), smut.end());
	BasicSignalView<const float, TIME_DOMAIN> v5 = AsView(s);
	BasicSignalView<const float, TIME_DOMAIN> v6 = AsConstView(s);
	BasicSignalView<const float, TIME_DOMAIN> v7 = AsConstView(s);

	REQUIRE(v1.size() == smut.size());
	REQUIRE(v2.size() == smut.size());
	REQUIRE(v3.size() == smut.size());
	REQUIRE(v4.size() == smut.size());
	REQUIRE(v5.size() == s.size());
	REQUIRE(v6.size() == s.size());
	REQUIRE(v7.size() == s.size());
}

TEST_CASE("Default construct", "[SignalView]") {
	BasicSignalView<float, TIME_DOMAIN> span;
	REQUIRE(span.empty());
	REQUIRE(span.size() == 0);

	BasicSignalView<std::complex<float>, TIME_DOMAIN> cspan;
	REQUIRE(cspan.empty());
	REQUIRE(cspan.size() == 0);
}


TEST_CASE("Whole span", "[SignalView]") {
	SignalF signal = { 1, 2, 3, 4, 5, 6 };

	BasicSignalView<float, TIME_DOMAIN> span{ signal };
	REQUIRE(span.size() == signal.size());
	REQUIRE(span[0] == 1);
	REQUIRE(span[5] == 6);
}


TEST_CASE("Partial span size", "[SignalView]") {
	SignalF signal = { 1, 2, 3, 4, 5, 6 };

	BasicSignalView<float, TIME_DOMAIN> span{ signal.begin() + 3, 2 };
	REQUIRE(span.size() == 2);
	REQUIRE(span[0] == 4);
	REQUIRE(span[1] == 5);
}


TEST_CASE("Partial span iterators", "[SignalView]") {
	SignalF signal = { 1, 2, 3, 4, 5, 6 };

	BasicSignalView<float, TIME_DOMAIN> span{ signal.begin() + 2, signal.begin() + 4 };
	REQUIRE(span.size() == 2);
	REQUIRE(span[0] == 3);
	REQUIRE(span[1] == 4);
}


TEST_CASE("data pointer", "[SignalView]") {
	SignalF signal = { 1, 2, 3, 4, 5, 6 };

	BasicSignalView<float, TIME_DOMAIN> span{ signal.begin() + 2, signal.begin() + 4 };
	REQUIRE(*span.data() == 3);
}


TEST_CASE("Real/Imag pointer", "[SignalView]") {
	using namespace std::complex_literals;
	SignalCF signal = { 1.f + 2.if,
						2.f + 3.if,
						3.f + 6.if,
						4.f + 7.if,
						5.f + 8.if,
						6.f + 9.if };

	BasicSignalView<std::complex<float>, TIME_DOMAIN> span{ signal.begin() + 2, signal.begin() + 4 };
	REQUIRE(span.data()->real() == 3);
	REQUIRE(span.data()->imag() == 6);
}


TEST_CASE("Constant span", "[SignalView]") {
	SignalF signal = { 1, 2, 3, 4, 5, 6 };

	BasicSignalView<const float, TIME_DOMAIN> span{ signal.begin() + 2, signal.begin() + 4 };
	REQUIRE(*span.data() == 3);

	SignalCF csignal = { 1.f + 2.if,
						 2.f + 3.if,
						 3.f + 6.if,
						 4.f + 7.if,
						 5.f + 8.if,
						 6.f + 9.if };

	BasicSignalView<const std::complex<float>, TIME_DOMAIN> cspan{ csignal.begin() + 2, csignal.begin() + 4 };
	REQUIRE(cspan.data()->real() == 3);
}
