#include <dspbb/Math/Statistics.hpp>
#include <dspbb/Primitives/Signal.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>
#include <random>


using namespace dspbb;
using namespace std::complex_literals;
using Catch::Approx;


TEST_CASE("CentralMoment #0 and #1", "[Statistics]") {
	Signal<float> s = { 2, 4, 4, 4, 5, 5, 7, 9 };
	REQUIRE(0 == CentralMoment(s, 0));
	REQUIRE(0 == CentralMoment(s, 1));
}

TEST_CASE("CentralMoment #2", "[Statistics]") {
	Signal<float> s = { 2, 4, 4, 4, 5, 5, 7, 9 };
	REQUIRE(Approx(4) == CentralMoment(s, 2));
}

TEST_CASE("CentralMoment #3", "[Statistics]") {
	Signal<float> s = { 2, 4, 4, 4, 5, 5, 7, 9 };
	REQUIRE(Approx(5.25) == CentralMoment(s, 3));
}

TEST_CASE("CentralMoment #4", "[Statistics]") {
	Signal<float> s = { 2, 4, 4, 4, 5, 5, 7, 9 };
	REQUIRE(Approx(44.5) == CentralMoment(s, 4));
}

TEST_CASE("CentralMoment #5", "[Statistics]") {
	Signal<float> s = { 2, 4, 4, 4, 5, 5, 7, 9 };
	REQUIRE(Approx(101.25) == CentralMoment(s, 5));
}

TEST_CASE("StandardizedMoment #1", "[Statistics]") {
	Signal<float> s = { 2, 4, 4, 4, 5, 5, 7, 9 };
	REQUIRE(Approx(0) == StandardizedMoment(s, 1));
}

TEST_CASE("StandardizedMoment #2", "[Statistics]") {
	Signal<float> s = { 2, 4, 4, 4, 5, 5, 7, 9 };
	REQUIRE(Approx(1) == StandardizedMoment(s, 2));
}

TEST_CASE("Standard deviation", "[Statistics]") {
	Signal<float> s = { 2, 4, 4, 4, 5, 5, 7, 9 };
	REQUIRE(Approx(2) == StandardDeviation(s));
}

TEST_CASE("Variance", "[Statistics]") {
	Signal<float> s = { 2, 4, 4, 4, 5, 5, 7, 9 };
	REQUIRE(Approx(4) == Variance(s));
}

TEST_CASE("Skewness", "[Statistics]") {
	Signal<float> s = { 2, 4, 4, 4, 5, 5, 7, 9 };
	REQUIRE(Approx(5.25f / 8.f) == Skewness(s));
}

TEST_CASE("Kurtosis", "[Statistics]") {
	Signal<float> s = { 2, 4, 4, 4, 5, 5, 7, 9 };
	REQUIRE(Approx(44.5f / 16.f) == Kurtosis(s));
}



TEST_CASE("Corrected standard deviation", "[Statistics]") {
	Signal<float> s = { 2, 4, 4, 4, 5, 5, 7, 9 };
	REQUIRE(Approx(2.138089935) == CorrectedStandardDeviation(s));
}

TEST_CASE("Corrected variance", "[Statistics]") {
	Signal<float> s = { 2, 4, 4, 4, 5, 5, 7, 9 };
	REQUIRE(Approx(4.571428571) == CorrectedVariance(s));
}

TEST_CASE("Corrected skewness", "[Statistics]") {
	Signal<float> s = { 2, 4, 4, 4, 5, 5, 7, 9 };
	REQUIRE(Approx(0.818487553) == CorrectedSkewness(s));
}

TEST_CASE("Corrected kurtosis", "[Statistics]") {
	Signal<float> s(1000000);
	std::mt19937 rne(762375);
	std::normal_distribution<float> rng;
	for (auto& v : s) {
		v = rng(rne);
	}
	// Kurtosis of normal distribution should be 3.
	REQUIRE(Approx(3.0f).epsilon(0.01f) == CorrectedKurtosis(s));
}



TEST_CASE("Sum", "[Statistics]") {
	Signal<float> s = { 1, 3, 2, 4, 5, 6, 7, 8, 9, 10 };
	REQUIRE(Approx(55) == Sum(s));
}


TEST_CASE("Mean", "[Statistics]") {
	Signal<float> s = { 1, 3, 2, 4, 5, 6, 7, 8, 9, 10 };
	REQUIRE(Approx(5.5f) == Mean(s));
}


TEST_CASE("SumSquare", "[Statistics]") {
	Signal<float> s = { 1, 3, 2, 4, 5, 6, 7, 8, 9, 10 };
	REQUIRE(Approx(385) == SumSquare(s));
}


TEST_CASE("MeanSquare", "[Statistics]") {
	Signal<float> s = { 1, 3, 2, 4, 5, 6, 7, 8, 9, 10 };
	REQUIRE(Approx(38.5f) == MeanSquare(s));
}


TEST_CASE("RootMeanSquare", "[Statistics]") {
	Signal<float> s = { 1, 3, 2, 4, 5, 6, 7, 8, 9, 10 };
	REQUIRE(Approx(std::sqrt(38.5f)) == RootMeanSquare(s));
}

TEST_CASE("Norm", "[Statistics]") {
	Signal<float> s = { 1, 3, 2, 4, 5, 6, 7, 8, 9, 10 };
	REQUIRE(Approx(std::sqrt(385.f)) == Norm(s));
}


TEST_CASE("Max", "[Statistics]") {
	Signal<float> s = { 1, 3, 2, 4, 8, 9, 10, 5, 6, 7 };
	REQUIRE(Approx(10) == Max(s));
}


TEST_CASE("Min", "[Statistics]") {
	Signal<float> s = { 1, 3, 2, 4, 8, 9, 10, 5, 6, 7 };
	REQUIRE(Approx(1) == Min(s));
}


TEST_CASE("Covariance self", "[Statistics]") {
	Signal<float> s = { 1, 3, 2, 4, 8, 9, 10, 5, 6, 7 };
	Signal<float> t = { 1, 3, 2, 4, 8, 9, 10, 5, 6, 7 };
	REQUIRE(Approx(Variance(s)) == Covariance(s, t));
}

TEST_CASE("Covariance anti", "[Statistics]") {
	Signal<float> base = { 1, 3, 2, 4, 8, 9, 10, 5, 6, 7 };
	Signal<float> s = base - Mean(base);
	Signal<float> t = Mean(base) - base;
	REQUIRE(Approx(-Variance(s)) == Covariance(s, t));
}

TEST_CASE("Covariance middle", "[Statistics]") {
	Signal<float> s = { 1, 3, 2, 4, 8, 9, 10, 5, 6, 7 };
	Signal<float> t = { 3, 4, 5, 6, 3, 7, 3, 7, 4, 5 };
	REQUIRE(Approx(0.15f) == Covariance(s, t));
}

TEST_CASE("Corrected covariance self", "[Statistics]") {
	Signal<float> s = { 1, 3, 2, 4, 8, 9, 10, 5, 6, 7 };
	Signal<float> t = { 1, 3, 2, 4, 8, 9, 10, 5, 6, 7 };
	REQUIRE(Approx(CorrectedVariance(s)) == CorrectedCovariance(s, t));
}

TEST_CASE("Correlation self", "[Statistics]") {
	Signal<float> s = { 1, 3, 2, 4, 8, 9, 10, 5, 6, 7 };
	Signal<float> t = { 1, 3, 2, 4, 8, 9, 10, 5, 6, 7 };
	REQUIRE(Approx(1) == Correlation(s, t));
}

TEST_CASE("Correlation anti", "[Statistics]") {
	Signal<float> base = { 1, 3, 2, 4, 8, 9, 10, 5, 6, 7 };
	Signal<float> s = base - Mean(base);
	Signal<float> t = Mean(base) - base;
	REQUIRE(Approx(-1) == Correlation(s, t));
}

TEST_CASE("Correlation middle", "[Statistics]") {
	Signal<float> s = { 1, 3, 2, 4, 8, 9, 10, 5, 6, 7 };
	Signal<float> t = { 3, 4, 5, 6, 3, 7, 3, 7, 4, 5 };
	REQUIRE(Approx(0.15f / 4.27f) == Correlation(s, t));
}
