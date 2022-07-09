#include "../TestUtils.hpp"

#include <dspbb/Math/Functions.hpp>
#include <dspbb/Math/OverlapAdd.hpp>
#include <dspbb/Math/Statistics.hpp>

#include <catch2/catch.hpp>

using namespace dspbb;
using namespace std::complex_literals;


TEST_CASE("OLA real-real central", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(3);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(7);
	const auto ola = OverlapAdd(signal, filter, CONV_CENTRAL, 16);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.size() == conv.size());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real central long", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(63);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(7);
	const auto ola = OverlapAdd(signal, filter, CONV_CENTRAL, 16);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.size() == conv.size());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real central big chunk", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(63);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(9);
	const auto ola = OverlapAdd(signal, filter, CONV_CENTRAL, 25);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.size() == conv.size());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real central small chunk", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(63);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(9);
	const auto ola = OverlapAdd(signal, filter, CONV_CENTRAL, 17);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.size() == conv.size());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real full", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(3);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(7);
	const auto ola = OverlapAdd(signal, filter, CONV_FULL, 16);
	const auto conv = Convolution(signal, filter, CONV_FULL);
	REQUIRE(ola.size() == conv.size());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real full long", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(63);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(7);
	const auto ola = OverlapAdd(signal, filter, CONV_FULL, 16);
	const auto conv = Convolution(signal, filter, CONV_FULL);
	REQUIRE(ola.size() == conv.size());
	const auto diff = ola - conv;
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real full big chunk", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(63);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(9);
	const auto ola = OverlapAdd(signal, filter, CONV_FULL, 25);
	const auto conv = Convolution(signal, filter, CONV_FULL);
	REQUIRE(ola.size() == conv.size());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-real full small chunk", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(63);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(9);
	const auto ola = OverlapAdd(signal, filter, CONV_FULL, 17);
	const auto conv = Convolution(signal, filter, CONV_FULL);
	REQUIRE(ola.size() == conv.size());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA real-complex", "[OverlapAdd]") {
	const auto signal = RandomSignal<float, TIME_DOMAIN>(107);
	const auto filter = RandomSignal<std::complex<float>, TIME_DOMAIN>(16);
	const auto ola = OverlapAdd(signal, filter, CONV_CENTRAL, 46);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.size() == conv.size());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA complex-real", "[OverlapAdd]") {
	const auto signal = RandomSignal<std::complex<float>, TIME_DOMAIN>(107);
	const auto filter = RandomSignal<float, TIME_DOMAIN>(16);
	const auto ola = OverlapAdd(signal, filter, CONV_CENTRAL, 46);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.size() == conv.size());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}

TEST_CASE("OLA complex-complex", "[OverlapAdd]") {
	const auto signal = RandomSignal<std::complex<float>, TIME_DOMAIN>(107);
	const auto filter = RandomSignal<std::complex<float>, TIME_DOMAIN>(16);
	const auto ola = OverlapAdd(signal, filter, CONV_CENTRAL, 46);
	const auto conv = Convolution(signal, filter, CONV_CENTRAL);
	REQUIRE(ola.size() == conv.size());
	REQUIRE(Max(Abs(ola - conv)) == Approx(0).margin(0.001f));
}



TEST_CASE("OLA Arbitrary offset middle", "[OverlapAdd]") {
	const auto signal = RandomSignal<std::complex<float>, TIME_DOMAIN>(107);
	const auto filter = RandomSignal<std::complex<float>, TIME_DOMAIN>(16);
	const auto ola = OverlapAdd(signal, filter, 24, 7, 33);
	const auto conv = Convolution(signal, filter, 24, 7);

	REQUIRE(ola.size() == conv.size());
	for (size_t i = 0; i < conv.size(); ++i) {
		REQUIRE(ola[i] == ApproxComplex(conv[i]).margin(1e-4f));
	}
}

TEST_CASE("OLA Arbitrary offset start", "[OverlapAdd]") {
	const auto signal = RandomSignal<std::complex<float>, TIME_DOMAIN>(107);
	const auto filter = RandomSignal<std::complex<float>, TIME_DOMAIN>(16);
	const auto ola = OverlapAdd(signal, filter, 0, 7, 31);
	const auto conv = Convolution(signal, filter, 0, 7);

	REQUIRE(ola.size() == conv.size());
	for (size_t i = 0; i < conv.size(); ++i) {
		REQUIRE(ola[i] == ApproxComplex(conv[i]).margin(1e-4f));
	}
}

TEST_CASE("OLA Arbitrary offset end", "[OverlapAdd]") {
	const auto signal = RandomSignal<std::complex<float>, TIME_DOMAIN>(107);
	const auto filter = RandomSignal<std::complex<float>, TIME_DOMAIN>(16);
	const auto ola = OverlapAdd(signal, filter, 100, 7, 33);
	const auto conv = Convolution(signal, filter, 100, 7);

	REQUIRE(ola.size() == conv.size());
	for (size_t i = 0; i < conv.size(); ++i) {
		REQUIRE(ola[i] == ApproxComplex(conv[i]).margin(1e-4f));
	}
}

TEST_CASE("OLA 3-operand full & central", "[OverlapAdd]") {
	const auto u = RandomSignal<std::complex<float>, TIME_DOMAIN>(107);
	const auto v = RandomSignal<std::complex<float>, TIME_DOMAIN>(16);
	const auto fullExpected = Convolution(v, u, CONV_FULL);
	const auto centralExpected = Convolution(v, u, CONV_CENTRAL);
	std::decay_t<decltype(fullExpected)> fullOut(fullExpected.size());
	std::decay_t<decltype(centralExpected)> centralOut(centralExpected.size());

	OverlapAdd(fullOut, u, v, CONV_FULL, 33);
	OverlapAdd(centralOut, u, v, CONV_CENTRAL, 33);

	for (size_t i = 0; i < fullOut.size(); ++i) {
		REQUIRE(fullOut[i] == ApproxComplex(fullExpected[i]).margin(1e-4f));
	}
	for (size_t i = 0; i < centralOut.size(); ++i) {
		REQUIRE(centralOut[i] == ApproxComplex(centralExpected[i]).margin(1e-4f));
	}
}

TEST_CASE("OLA optimal theoretical FFT size", "[OverlapAdd]") {
	const double s1 = impl::ola::OptimalTheoreticalSize(12, 6, 1, 2);
	REQUIRE(s1 == Approx(65.114).margin(0.001f));

	const double s2 = impl::ola::OptimalTheoreticalSize(30, 6, 1, 2);
	REQUIRE(s2 == Approx(195.815).margin(0.001f));

	const double s3 = impl::ola::OptimalTheoreticalSize(1024, 6, 1, 2);
	REQUIRE(s3 == Approx(10789.169).margin(0.001f));

	const double s4 = impl::ola::OptimalTheoreticalSize(6144, 6, 1, 2);
	REQUIRE(s4 == Approx(76793.054).margin(0.001f));
}

TEST_CASE("OLA optimal practical FFT size", "[OverlapAdd]") {
	const size_t s1 = impl::ola::OptimalPracticalSize(55000, 12, 6, 1, 2);
	REQUIRE(s1 == 128);

	const size_t s2 = impl::ola::OptimalPracticalSize(55000, 30, 6, 1, 2);
	REQUIRE(s2 == 256);

	const size_t s3 = impl::ola::OptimalPracticalSize(55000, 1024, 6, 1, 2);
	REQUIRE(s3 == 16384);

	const size_t s4 = impl::ola::OptimalPracticalSize(550000, 6144, 6, 1, 2);
	REQUIRE(s4 == 131072);
}

TEST_CASE("OLA optimal practical FFT size short signal", "[OverlapAdd]") {
	const size_t s1 = impl::ola::OptimalPracticalSize(49, 12, 6, 1, 2);
	REQUIRE(s1 == 60);

	const size_t s2 = impl::ola::OptimalPracticalSize(84, 12, 6, 1, 2);
	REQUIRE(s2 == 95);

	const size_t s3 = impl::ola::OptimalPracticalSize(86, 12, 6, 1, 2);
	REQUIRE(s3 == 128);
}