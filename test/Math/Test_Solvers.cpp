#include "dspbb/Utility/Numbers.hpp"
#include <dspbb/Math/Solvers.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>


using namespace dspbb;
using Catch::Approx;


TEST_CASE("Bisect exponential", "[Solvers]") {
	const auto f = [](auto x) { return std::exp(x) - 2; };
	const auto x0 = Bisect(f, -1.0, 1.0);
	REQUIRE(x0 == Approx(std::log(2)));
}

TEST_CASE("Bisect sine", "[Solvers]") {
	const auto f = [](auto x) { return std::sin(x); };
	const auto x0 = Bisect(f, 3.0, 3.5);
	REQUIRE(x0 == Approx(pi_v<double>));
}

TEST_CASE("Bisect reversed arguments", "[Solvers]") {
	const auto f = [](auto x) { return std::exp(x) - 2; };
	const auto x0 = Bisect(f, 1.0, -1.0);
	REQUIRE(x0 == Approx(std::log(2)));
}

TEST_CASE("Bisect not containing root termination", "[Solvers]") {
	const auto f = [](auto x) { return std::exp(x) - 2; };
	const auto x0 = Bisect(f, -2.0, -1.0);
	REQUIRE((x0 == Approx(-2.0) || x0 == Approx(-1.0)));
}


TEST_CASE("Newton exponential", "[Solvers]") {
	const auto f = [](auto x) { return std::exp(x) - 2; };
	const auto df = [](auto x) { return std::exp(x); };
	const auto x0 = NewtonRaphson(f, df, 0.7);
	REQUIRE(x0 == Approx(std::log(2)));
}

TEST_CASE("Newton sine", "[Solvers]") {
	const auto f = [](auto x) { return std::sin(x); };
	const auto df = [](auto x) { return std::cos(x); };
	const auto x0 = NewtonRaphson(f, df, 3.0);
	REQUIRE(x0 == Approx(pi_v<double>));
}