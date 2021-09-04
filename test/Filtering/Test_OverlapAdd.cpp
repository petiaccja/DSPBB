#include <catch2/catch.hpp>
#include <dspbb/Filtering/OverlapAdd.hpp>
#include <dspbb/Math/Functions.hpp>
#include <dspbb/Math/Statistics.hpp>
#include <random>

using namespace dspbb;
using namespace std::complex_literals;

static std::mt19937_64 rne(723574);
static std::uniform_real_distribution<float> rng;

TimeSignal<float> RandomSignal(size_t length) {
	TimeSignal<float> s(length);
	for (auto& v : s) {
		v = rng(rne);
	}
	return s;
}

TimeSignal<std::complex<float>> RandomComplexSignal(size_t length) {
	TimeSignal<std::complex<float>> s(length);
	for (auto& v : s) {
		v.real(rng(rne));
		v.imag(rng(rne));
	}
	return s;
}

TEST_CASE("OLA real-real central", "[OverlapAdd]") {
	const auto signal = RandomSignal(3);
	const auto filter = RandomSignal(7);
	const auto ola = OverlapAdd(signal, filter, 4, 12, convolution::central);
	const auto conv = Convolution(signal, filter, convolution::central);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real central long", "[OverlapAdd]") {
	const auto signal = RandomSignal(63);
	const auto filter = RandomSignal(7);
	const auto ola = OverlapAdd(signal, filter, 4, 12, convolution::central);
	const auto conv = Convolution(signal, filter, convolution::central);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real central big chunk", "[OverlapAdd]") {
	const auto signal = RandomSignal(63);
	const auto filter = RandomSignal(9);
	const auto ola = OverlapAdd(signal, filter, 17, 8, convolution::central);
	const auto conv = Convolution(signal, filter, convolution::central);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real full", "[OverlapAdd]") {
	const auto signal = RandomSignal(3);
	const auto filter = RandomSignal(7);
	const auto ola = OverlapAdd(signal, filter, 4, 12, convolution::full);
	const auto conv = Convolution(signal, filter, convolution::full);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real full long", "[OverlapAdd]") {
	const auto signal = RandomSignal(63);
	const auto filter = RandomSignal(7);
	const auto ola = OverlapAdd(signal, filter, 4, 12, convolution::full);
	const auto conv = Convolution(signal, filter, convolution::full);
	REQUIRE(ola.Length() == conv.Length());
	const auto diff = ola - conv;
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real full big chunk", "[OverlapAdd]") {
	const auto signal = RandomSignal(63);
	const auto filter = RandomSignal(9);
	const auto ola = OverlapAdd(signal, filter, 17, 8, convolution::full);
	const auto conv = Convolution(signal, filter, convolution::full);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-complex", "[OverlapAdd]") {
	const auto signal = RandomSignal(107);
	const auto filter = RandomComplexSignal(16);
	const auto ola = OverlapAdd(signal, filter, 31, 15, convolution::central);
	const auto conv = Convolution(signal, filter, convolution::central);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA complex-real", "[OverlapAdd]") {
	const auto signal = RandomComplexSignal(107);
	const auto filter = RandomSignal(16);
	const auto ola = OverlapAdd(signal, filter, 31, 15, convolution::central);
	const auto conv = Convolution(signal, filter, convolution::central);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA complex-complex", "[OverlapAdd]") {
	const auto signal = RandomComplexSignal(107);
	const auto filter = RandomComplexSignal(16);
	const auto ola = OverlapAdd(signal, filter, 31, 15, convolution::central);
	const auto conv = Convolution(signal, filter, convolution::central);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}