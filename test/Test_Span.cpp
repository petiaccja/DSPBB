#include <dspbb/Primitives/Span.hpp>

#include <Catch2/catch.hpp>
#include <complex>

using namespace enl;
using namespace std::complex_literals;

TEST_CASE("Default construct", "[AudioFramework:Span]") {
	Span<float, TIME_DOMAIN> span;
	REQUIRE(span.Empty());
	REQUIRE(span.Size() == 0);

	Span<std::complex<float>, TIME_DOMAIN> cspan;
	REQUIRE(cspan.Empty());
	REQUIRE(cspan.Size() == 0);
}


TEST_CASE("Whole span", "[AudioFramework:Span]") {
	TimeSignalF signal = { 1, 2, 3, 4, 5, 6 };

	Span<float, TIME_DOMAIN> span{ signal };
	REQUIRE(span.Size() == signal.Size());
	REQUIRE(span[0] == 1);
	REQUIRE(span[5] == 6);
}


TEST_CASE("Partial span size", "[AudioFramework:Span]") {
	TimeSignalF signal = { 1, 2, 3, 4, 5, 6 };

	Span<float, TIME_DOMAIN> span{ signal.begin() + 3, 2 };
	REQUIRE(span.Size() == 2);
	REQUIRE(span[0] == 4);
	REQUIRE(span[1] == 5);
}


TEST_CASE("Partial span iterators", "[AudioFramework:Span]") {
	TimeSignalF signal = { 1, 2, 3, 4, 5, 6 };

	Span<float, TIME_DOMAIN> span{ signal.begin() + 2, signal.begin() + 4 };
	REQUIRE(span.Size() == 2);
	REQUIRE(span[0] == 3);
	REQUIRE(span[1] == 4);
}


TEST_CASE("Data pointer", "[AudioFramework:Span]") {
	TimeSignalF signal = { 1, 2, 3, 4, 5, 6 };

	Span<float, TIME_DOMAIN> span{ signal.begin() + 2, signal.begin() + 4 };
	auto* v = span.Data();
	REQUIRE(*span.Data() == 3);
}


TEST_CASE("Real/Imag pointer", "[AudioFramework:Span]") {
	using namespace std::complex_literals;
	TimeSignalCF signal = { 1.f + 2.if,
							2.f + 3.if,
							3.f + 6.if,
							4.f + 7.if,
							5.f + 8.if,
							6.f + 9.if };

	Span<std::complex<float>, TIME_DOMAIN> span{ signal.begin() + 2, signal.begin() + 4 };
	REQUIRE(span.Data()->real() == 3);
	REQUIRE(span.Data()->imag() == 6);
}


TEST_CASE("Constant span", "[AudioFramework:Span]") {
	TimeSignalF signal = { 1, 2, 3, 4, 5, 6 };

	Span<const float, TIME_DOMAIN> span{ signal.begin() + 2, signal.begin() + 4 };
	REQUIRE(*span.Data() == 3);

	TimeSignalCF csignal = { 1.f + 2.if,
							2.f + 3.if,
							3.f + 6.if,
							4.f + 7.if,
							5.f + 8.if,
							6.f + 9.if };

	Span<const std::complex<float>, TIME_DOMAIN> cspan{ csignal.begin() + 2, csignal.begin() + 4 };
	REQUIRE(cspan.Data()->real() == 3);
}

