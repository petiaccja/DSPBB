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
	const auto ola = OverlapAdd(signal, filter, CONV_CENTRAL, 16);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real central long", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(63);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(7);
	const auto ola = OverlapAdd(signal, filter, CONV_CENTRAL, 16);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real central big chunk", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(63);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(9);
	const auto ola = OverlapAdd(signal, filter, CONV_CENTRAL, 25);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real central small chunk", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(63);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(9);
	const auto ola = OverlapAdd(signal, filter, CONV_CENTRAL, 17);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real full", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(3);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(7);
	const auto ola = OverlapAdd(signal, filter, CONV_FULL, 16);
	const auto conv = Convolution(signal, filter, CONV_FULL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real full long", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(63);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(7);
	const auto ola = OverlapAdd(signal, filter, CONV_FULL, 16);
	const auto conv = Convolution(signal, filter, CONV_FULL);
	REQUIRE(ola.Length() == conv.Length());
	const auto diff = ola - conv;
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real full big chunk", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(63);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(9);
	const auto ola = OverlapAdd(signal, filter, CONV_FULL, 25);
	const auto conv = Convolution(signal, filter, CONV_FULL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real full small chunk", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(63);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(9);
	const auto ola = OverlapAdd(signal, filter, CONV_FULL, 17);
	const auto conv = Convolution(signal, filter, CONV_FULL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-complex", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(107);
	const auto filter = RandomSignal<std::complex<float>, TIME_DOMAIN>(16);
	const auto ola = OverlapAdd(signal, filter, CONV_CENTRAL, 46);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA complex-real", "[OverlapAdd]") {
	const auto signal = RandomSignal<std::complex<float>, TIME_DOMAIN>(107);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(16);
	const auto ola = OverlapAdd(signal, filter, CONV_CENTRAL, 46);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA complex-complex", "[OverlapAdd]") {
	const auto signal = RandomSignal<std::complex<float>, TIME_DOMAIN>(107);
	const auto filter = RandomSignal<std::complex<float>, TIME_DOMAIN>(16);
	const auto ola = OverlapAdd(signal, filter, CONV_CENTRAL, 46);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.Length() == conv.Length());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}




TEST_CASE("OLA Arbitrary offset middle", "[OverlapAdd]") {
	const auto signal = RandomSignal<std::complex<float>, TIME_DOMAIN>(107);
	const auto filter = RandomSignal<std::complex<float>, TIME_DOMAIN>(16);
	const auto ola = OverlapAdd(signal, filter, 24, 7, 33);
	const auto conv = Convolution(signal, filter, 24, 7);

	REQUIRE(ola.Length() == conv.Length());
	for (size_t i = 0; i < conv.Length(); ++i) {
		REQUIRE(ola[i] == ApproxComplex(conv[i]).margin(1e-4f));
	}
}

TEST_CASE("OLA Arbitrary offset start", "[OverlapAdd]") {
	const auto signal = RandomSignal<std::complex<float>, TIME_DOMAIN>(107);
	const auto filter = RandomSignal<std::complex<float>, TIME_DOMAIN>(16);
	const auto ola = OverlapAdd(signal, filter, 0, 7, 31);
	const auto conv = Convolution(signal, filter, 0, 7);

	REQUIRE(ola.Length() == conv.Length());
	for (size_t i = 0; i < conv.Length(); ++i) {
		REQUIRE(ola[i] == ApproxComplex(conv[i]).margin(1e-4f));
	}
}

TEST_CASE("OLA Arbitrary offset end", "[OverlapAdd]") {
	const auto signal = RandomSignal<std::complex<float>, TIME_DOMAIN>(107);
	const auto filter = RandomSignal<std::complex<float>, TIME_DOMAIN>(16);
	const auto ola = OverlapAdd(signal, filter, 100, 7, 33);
	const auto conv = Convolution(signal, filter, 100, 7);

	REQUIRE(ola.Length() == conv.Length());
	for (size_t i = 0; i < conv.Length(); ++i) {
		REQUIRE(ola[i] == ApproxComplex(conv[i]).margin(1e-4f));
	}
}

TEST_CASE("OLA Arbitrary offset out of bounds", "[OverlapAdd]") {
	const auto signal = RandomSignal<std::complex<float>, TIME_DOMAIN>(107);
	const auto filter = RandomSignal<std::complex<float>, TIME_DOMAIN>(16);

	REQUIRE_THROWS(OverlapAdd(signal, filter, 95, 30, 33));
	REQUIRE_THROWS(OverlapAdd(signal, filter, 0, 190, 33));
}

TEST_CASE("OLA 3-operand full & central", "[OverlapAdd]") {
	const auto u = RandomSignal<std::complex<float>, TIME_DOMAIN>(107);
	const auto v = RandomSignal<std::complex<float>, TIME_DOMAIN>(16);
	const auto fullExpected = Convolution(v, u, CONV_FULL);
	const auto centralExpected = Convolution(v, u, CONV_CENTRAL);
	std::decay_t<decltype(fullExpected)> fullOut(fullExpected.Size());
	std::decay_t<decltype(centralExpected)> centralOut(centralExpected.Size());

	OverlapAdd(fullOut, u, v, CONV_FULL, 33);
	OverlapAdd(centralOut, u, v, CONV_CENTRAL, 33);

	for (size_t i = 0; i < fullOut.Length(); ++i) {
		REQUIRE(fullOut[i] == ApproxComplex(fullExpected[i]).margin(1e-4f));
	}
	for (size_t i = 0; i < centralOut.Length(); ++i) {
		REQUIRE(centralOut[i] == ApproxComplex(centralExpected[i]).margin(1e-4f));
	}
}

TEST_CASE("OLA too small chunk size", "[OverlapAdd]") {
	const auto u = RandomSignal<std::complex<float>, TIME_DOMAIN>(107);
	const auto v = RandomSignal<std::complex<float>, TIME_DOMAIN>(16);
	REQUIRE_THROWS(OverlapAdd(u, v, CONV_FULL, 30));
	REQUIRE_NOTHROW(OverlapAdd(u, v, CONV_CENTRAL, 31));
}