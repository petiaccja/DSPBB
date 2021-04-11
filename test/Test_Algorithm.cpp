#include <catch2/catch.hpp>
#include <complex>
#include <dspbb/Utility/Algorithm.hpp>

using namespace dspbb;
using namespace std::complex_literals;

TEST_CASE("Apply", "[Algorithm]") {
	Signal<float, TIME_DOMAIN> signal{ 1.f, 4.f, 9.f };
	auto result = Apply(signal, static_cast<float (*)(float)>(std::sqrt));
	REQUIRE((std::is_same<decltype(result)::value_type, float>::value));
	REQUIRE(result[0] == Approx(1.0f));
	REQUIRE(result[1] == Approx(2.0f));
	REQUIRE(result[2] == Approx(3.0f));
}


TEST_CASE("Apply convert", "[Algorithm]") {
	Signal<std::complex<float>, TIME_DOMAIN> signal{ 1.f, -4.f, 9.if };
	auto result = Apply(signal, static_cast<float (*)(const std::complex<float>&)>(std::abs));
	REQUIRE((std::is_same<decltype(result)::value_type, float>::value));
	REQUIRE(result[0] == Approx(1.0f));
	REQUIRE(result[1] == Approx(4.0f));
	REQUIRE(result[2] == Approx(9.0f));
}
