#include <catch2/catch.hpp>
#include <dspbb/Filtering/OverlapAdd.hpp>
#include <dspbb/Math/Functions.hpp>
#include <dspbb/Math/Statistics.hpp>
#include "../TestUtils.hpp"

using namespace dspbb;
using namespace std::complex_literals;


TEST_CASE("OLA real-real central", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(3);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(7);
	const auto ola = OverlapAdd(signal, filter, 4, 12, CONV_CENTRAL);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real central long", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(63);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(7);
	const auto ola = OverlapAdd(signal, filter, 4, 12, CONV_CENTRAL);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real central big chunk", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(63);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(9);
	const auto ola = OverlapAdd(signal, filter, 17, 8, CONV_CENTRAL);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real full", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(3);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(7);
	const auto ola = OverlapAdd(signal, filter, 4, 12, CONV_FULL);
	const auto conv = Convolution(signal, filter, CONV_FULL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real full long", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(63);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(7);
	const auto ola = OverlapAdd(signal, filter, 4, 12, CONV_FULL);
	const auto conv = Convolution(signal, filter, CONV_FULL);
	REQUIRE(ola.Length() == conv.Length());
	const auto diff = ola - conv;
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real full big chunk", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(63);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(9);
	const auto ola = OverlapAdd(signal, filter, 17, 8, CONV_FULL);
	const auto conv = Convolution(signal, filter, CONV_FULL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-complex", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(107);
	const auto filter = RandomSignal<std::complex<float>, TIME_DOMAIN>(16);
	const auto ola = OverlapAdd(signal, filter, 31, 15, CONV_CENTRAL);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA complex-real", "[OverlapAdd]") {
	const auto signal = RandomSignal<std::complex<float>, TIME_DOMAIN>(107);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(16);
	const auto ola = OverlapAdd(signal, filter, 31, 15, CONV_CENTRAL);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA complex-complex", "[OverlapAdd]") {
	const auto signal = RandomSignal<std::complex<float>, TIME_DOMAIN>(107);
	const auto filter = RandomSignal<std::complex<float>, TIME_DOMAIN>(16);
	const auto ola = OverlapAdd(signal, filter, 31, 15, CONV_CENTRAL);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}