#include <Catch2/catch.hpp>
#include <complex>
#include <dspbb/Math/Statistics.hpp>

using namespace dspbb;
using namespace std::complex_literals;




TEST_CASE("Sum", "[Statistics]") {
	TimeSignal<float> s = { 1, 3, 2 };
	REQUIRE(Approx(6) == Sum(s));
}


TEST_CASE("Mean", "[Statistics]") {
	TimeSignal<float> s = { 1, 3, 2 };
	REQUIRE(Approx(2) == Mean(s));
}


TEST_CASE("SumSquare", "[Statistics]") {
	TimeSignal<float> s = { std::sqrt(2.0f), 3, 4 };
	REQUIRE(Approx(27) == SumSquare(s));
}


TEST_CASE("MeanSquare", "[Statistics]") {
	TimeSignal<float> s = { std::sqrt(2.0f), 3, 4 };
	REQUIRE(Approx(9) == MeanSquare(s));
}


TEST_CASE("RootMeanSquare", "[Statistics]") {
	TimeSignal<float> s = { std::sqrt(2.0f), 3, 4 };
	REQUIRE(Approx(3) == RootMeanSquare(s));
}


TEST_CASE("Standard deviation", "[Statistics]") {
	TimeSignal<float> s = { 2, 4, 4, 4, 5, 5, 7, 9 };
	REQUIRE(Approx(2) == StandardDeviation(s));
}


TEST_CASE("Norm", "[Statistics]") {
	TimeSignal<float> s = { 3, 4 };
	REQUIRE(Approx(5) == Norm(s));
}


TEST_CASE("Max", "[Statistics]") {
	TimeSignal<float> s = { 1, 3, 2 };
	REQUIRE(Approx(3) == Max(s));
}


TEST_CASE("Min", "[Statistics]") {
	TimeSignal<float> s = { 1, 3, 2 };
	REQUIRE(Approx(1) == Min(s));
}
