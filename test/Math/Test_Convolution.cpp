#include <dspbb/Math/Convolution.hpp>

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <complex>


using namespace dspbb;
using namespace std::complex_literals;


constexpr std::array<float, 20> ur = { 1, 3, 7, 2, 9, 2, 5, 3, 7, 2, 4, 7, 3, 6, 3, 9, 3, 5, 3, 5 };
constexpr std::array<float, 12> vr = { 4, 3, 5, 2, 6, 3, 2, 7, 8, 5, 3, 3 };
constexpr std::array<float, 9> urvr_central = { 227, 244, 238, 207, 270, 219, 242, 223, 259 };
constexpr std::array<float, 31> urvr_full = {
	4, 15, 42, 46, 89, 80, 128, 101, 169, 175, 205, 227, 244, 238, 207, 270,
	219, 242, 223, 259, 210, 205, 196, 184, 152, 122, 120, 79, 49, 24, 15
};

const std::array<std::complex<float>, 20> uc = {
	8.f + 5.if, 8.f + 5.if, 4.f + 7.if, 7.f + 8.if, 2.f + 8.if, 8.f + 3.if, 1.f + 7.if,
	3.f + 7.if, 1.f + 2.if, 1.f + 2.if, 9.f + 5.if, 7.f + 10.if, 4.f + 4.if,
	10.f + 6.if, 1.f + 3.if, 5.f + 8.if, 4.f + 3.if, 8.f + 6.if, 8.f + 7.if, 2.f + 9.if
};
const std::array<std::complex<float>, 12> vc = {
	10.f + 3.if, 6.f + 7.if, 2.f + 5.if, 2.f + 4.if, 3.f + 9.if, 9.f + 6.if,
	3.f + 6.if, 9.f + 10.if, 3.f + 3.if, 10.f + 8.if, 4.f + 8.if, 2.f + 4.if
};
const std::array<std::complex<float>, 9> ucvc_central{
	-129.f + 770.if, -96.f + 722.if, -157.f + 641.if, -128.f + 650.if, -123.f + 646.if,
	-124.f + 642.if, -74.f + 663.if, -11.f + 688.if, -79.f + 721.if
};

const std::array<std::complex<float>, 31> ucvc_full{
	65.f + 74.if, 78.f + 160.if, 23.f + 218.if, 8.f + 263.if, -70.f + 346.if, 2.f + 377.if,
	-53.f + 430.if, -102.f + 560.if, -106.f + 508.if, -104.f + 576.if, -5.f + 645.if,
	-129.f + 770.if, -96.f + 722.if, -157.f + 641.if, -128.f + 650.if, -123.f + 646.if,
	-124.f + 642.if, -74.f + 663.if, -11.f + 688.if, -79.f + 721.if, -101.f + 762.if,
	-153.f + 568.if, -72.f + 624.if, -90.f + 494.if, -77.f + 509.if, -78.f + 368.if,
	-67.f + 350.if, -17.f + 277.if, -84.f + 242.if, -76.f + 98.if, -32.f + 26.if
};


TEST_CASE("Real central", "[Convolution]") {
	Signal<float> u{ ur.begin(), ur.end() };
	Signal<float> v{ vr.begin(), vr.end() };
	Signal<float> expected = { urvr_central.begin(), urvr_central.end() };

	auto result = Convolution(u, v, CONV_CENTRAL);

	REQUIRE(result.size() == expected.size());
	for (size_t i = 0; i < expected.size(); ++i) {
		REQUIRE(result[i] == expected[i]);
	}
}

TEST_CASE("Real full", "[Convolution]") {
	Signal<float> u{ ur.begin(), ur.end() };
	Signal<float> v{ vr.begin(), vr.end() };
	Signal<float> expected = { urvr_full.begin(), urvr_full.end() };

	auto result = Convolution(u, v, CONV_FULL);

	REQUIRE(result.size() == expected.size());
	for (size_t i = 0; i < expected.size(); ++i) {
		REQUIRE(result[i] == expected[i]);
	}
}


TEST_CASE("Real-complex central BR", "[Convolution]") {
	Signal<float> u{ ur.begin(), ur.end() };
	Signal<std::complex<float>> v{ vr.begin(), vr.end() };
	Signal<std::complex<float>> expected = { urvr_central.begin(), urvr_central.end() };

	auto result = Convolution(u, v, CONV_CENTRAL);

	REQUIRE(result.size() == expected.size());
	for (size_t i = 0; i < expected.size(); ++i) {
		REQUIRE(result[i] == expected[i]);
	}
}

TEST_CASE("Complex-real central RB", "[Convolution]") {
	Signal<std::complex<float>> u{ ur.begin(), ur.end() };
	Signal<float> v{ vr.begin(), vr.end() };
	Signal<std::complex<float>> expected = { urvr_central.begin(), urvr_central.end() };

	auto result = Convolution(u, v, CONV_CENTRAL);

	REQUIRE(result.size() == expected.size());
	for (size_t i = 0; i < expected.size(); ++i) {
		REQUIRE(result[i] == expected[i]);
	}
}

TEST_CASE("Complex-complex central", "[Convolution]") {
	Signal<std::complex<float>> u{ uc.begin(), uc.end() };
	Signal<std::complex<float>> v{ vc.begin(), vc.end() };
	Signal<std::complex<float>> expected = { ucvc_central.begin(), ucvc_central.end() };

	auto result = Convolution(u, v, CONV_CENTRAL);

	REQUIRE(result.size() == expected.size());
	for (size_t i = 0; i < expected.size(); ++i) {
		REQUIRE(result[i] == expected[i]);
	}
}

TEST_CASE("Real-world signal", "[Convolution]") {
	Signal<float> u;
	Signal<float> v;

	u.resize(1000, 0.0f);
	for (int i = 0; i < 20; ++i) {
		u[200 + i] = float(i);
	}

	v.resize(51, 0.0f);
	v[0] = 1;
	v[50] = 1;

	auto result = Convolution(u, v, CONV_CENTRAL);

	REQUIRE(result.size() == 950);
	REQUIRE(result[145] == 0);
	REQUIRE(result[151] == 1);
	REQUIRE(result[169] == 19);
	REQUIRE(result[185] == 0);
	REQUIRE(result[201] == 1);
	REQUIRE(result[219] == 19);
	REQUIRE(result[225] == 0);
}



TEST_CASE("Arbitrary offset middle", "[Convolution]") {
	Signal<float> u{ ur.begin(), ur.end() };
	Signal<float> v{ vr.begin(), vr.end() };
	Signal<float> expected = { urvr_full.begin(), urvr_full.end() };

	const auto result = Convolution(u, v, 4, 6);

	REQUIRE(result.size() == 6);
	for (size_t i = 0; i < result.size(); ++i) {
		REQUIRE(result[i] == expected[i + 4]);
	}
}

TEST_CASE("Arbitrary offset start", "[Convolution]") {
	Signal<float> u{ ur.begin(), ur.end() };
	Signal<float> v{ vr.begin(), vr.end() };
	Signal<float> expected = { urvr_full.begin(), urvr_full.end() };

	const auto result = Convolution(u, v, 0, 6);

	REQUIRE(result.size() == 6);
	for (size_t i = 0; i < result.size(); ++i) {
		REQUIRE(result[i] == expected[i]);
	}
}

TEST_CASE("Arbitrary offset end", "[Convolution]") {
	Signal<float> u{ ur.begin(), ur.end() };
	Signal<float> v{ vr.begin(), vr.end() };
	Signal<float> expected = { urvr_full.begin(), urvr_full.end() };

	const auto result = Convolution(u, v, 25, 6);

	REQUIRE(result.size() == 6);
	for (size_t i = 0; i < result.size(); ++i) {
		REQUIRE(result[i] == expected[i + 25]);
	}
}

TEST_CASE("3-operand full & central", "[Convolution]") {
	Signal<float> u{ ur.begin(), ur.end() };
	Signal<float> v{ vr.begin(), vr.end() };
	Signal<float> fullExpected = { urvr_full.begin(), urvr_full.end() };
	Signal<float> centralExpected = { urvr_central.begin(), urvr_central.end() };
	Signal<float> fullOut(fullExpected.size());
	Signal<float> centralOut(centralExpected.size());

	Convolution(fullOut, u, v, CONV_FULL);
	Convolution(centralOut, u, v, CONV_CENTRAL);

	for (size_t i = 0; i < fullOut.size(); ++i) {
		REQUIRE(fullOut[i] == fullExpected[i]);
	}
	for (size_t i = 0; i < centralOut.size(); ++i) {
		REQUIRE(centralOut[i] == centralExpected[i]);
	}
}
