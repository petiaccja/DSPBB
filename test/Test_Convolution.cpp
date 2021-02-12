#include <Catch2/catch.hpp>
#include <array>
#include <complex>
#include <dspbb/Filtering/Convolution.hpp>

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
	8. + 5.i, 8. + 5.i, 4. + 7.i, 7. + 8.i, 2. + 8.i, 8. + 3.i, 1. + 7.i,
	3. + 7.i, 1. + 2.i, 1. + 2.i, 9. + 5.i, 7. + 10.i, 4. + 4.i,
	10. + 6.i, 1. + 3.i, 5. + 8.i, 4. + 3.i, 8. + 6.i, 8. + 7.i, 2. + 9.i
};
const std::array<std::complex<float>, 12> vc = {
	10. + 3.i, 6. + 7.i, 2. + 5.i, 2. + 4.i, 3. + 9.i, 9. + 6.i,
	3. + 6.i, 9. + 10.i, 3. + 3.i, 10. + 8.i, 4. + 8.i, 2. + 4.i
};
const std::array<std::complex<float>, 9> ucvc_central{
	-129. + 770.i, -96. + 722.i, -157. + 641.i, -128. + 650.i, -123. + 646.i,
	-124. + 642.i, -74. + 663.i, -11. + 688.i, -79. + 721.i
};

const std::array<std::complex<float>, 31> ucvc_full{
	65. + 74.i, 78. + 160.i, 23. + 218.i, 8.0 + 263.i, -70. + 346.i, 2.0 + 377.i,
	-53. + 430.i, -102. + 560.i, -106. + 508.i, -104. + 576.i, -5.0 + 645.i,
	-129. + 770.i, -96. + 722.i, -157. + 641.i, -128. + 650.i, -123. + 646.i,
	-124. + 642.i, -74. + 663.i, -11. + 688.i, -79. + 721.i, -101. + 762.i,
	-153. + 568.i, -72. + 624.i, -90. + 494.i, -77. + 509.i, -78. + 368.i,
	-67. + 350.i, -17. + 277.i, -84. + 242.i, -76. + 98.i, -32. + 26.i
};


TEST_CASE("Real central", "[AudioFramework:Convolution]") {
	TimeSignal<float> u{ ur.begin(), ur.end() };
	TimeSignal<float> v{ vr.begin(), vr.end() };
	TimeSignal<float> expected = { urvr_central.begin(), urvr_central.end() };

	auto result = Convolution(u, v, convolution::central);

	REQUIRE(result.Length() == expected.Length());
	for (size_t i = 0; i < expected.Length(); ++i) {
		REQUIRE(result[i] == expected[i]);
	}	
}

TEST_CASE("Real full", "[AudioFramework:Convolution]") {
	TimeSignal<float> u{ ur.begin(), ur.end() };
	TimeSignal<float> v{ vr.begin(), vr.end() };
	TimeSignal<float> expected = { urvr_full.begin(), urvr_full.end() };

	auto result = Convolution(u, v, convolution::full);

	REQUIRE(result.Length() == expected.Length());
	for (size_t i = 0; i < expected.Length(); ++i) {
		REQUIRE(result[i] == expected[i]);
	}
}


TEST_CASE("Real-complex central BR", "[AudioFramework:Convolution]") {
	TimeSignal<float> u{ ur.begin(), ur.end() };
	TimeSignal<std::complex<float>> v{ vr.begin(), vr.end() };
	TimeSignal<std::complex<float>> expected = { urvr_central.begin(), urvr_central.end() };

	auto result = Convolution(u, v, convolution::central);

	REQUIRE(result.Length() == expected.Length());
	for (size_t i = 0; i < expected.Length(); ++i) {
		REQUIRE(result[i] == expected[i]);
	}
}

TEST_CASE("Complex-real central RB", "[AudioFramework:Convolution]") {
	TimeSignal<std::complex<float>> u{ ur.begin(), ur.end() };
	TimeSignal<float> v{ vr.begin(), vr.end() };
	TimeSignal<std::complex<float>> expected = { urvr_central.begin(), urvr_central.end() };

	auto result = Convolution(u, v, convolution::central);

	REQUIRE(result.Length() == expected.Length());
	for (size_t i = 0; i < expected.Length(); ++i) {
		REQUIRE(result[i] == expected[i]);
	}
}

TEST_CASE("Complex-complex central", "[AudioFramework:Convolution]") {
	TimeSignal<std::complex<float>> u{ uc.begin(), uc.end() };
	TimeSignal<std::complex<float>> v{ vc.begin(), vc.end() };
	TimeSignal<std::complex<float>> expected = { ucvc_central.begin(), ucvc_central.end() };

	auto result = Convolution(u, v, convolution::central);

	REQUIRE(result.Length() == expected.Length());
	for (size_t i = 0; i < expected.Length(); ++i) {
		REQUIRE(result[i] == expected[i]);
	}
}


TEST_CASE("Different types", "[AudioFramework:Convolution]") {
	TimeSignal<float> u{ ur.begin(), ur.end() };
	TimeSignal<std::complex<double>> v{ vr.begin(), vr.end() };
	TimeSignal<std::complex<double>> expected = { urvr_central.begin(), urvr_central.end() };

	auto result = Convolution(u, v, convolution::central);

	REQUIRE(result.Length() == expected.Length());
	for (size_t i = 0; i < expected.Length(); ++i) {
		REQUIRE(result[i] == expected[i]);
	}
}


TEST_CASE("Real-world signal", "[AudioFramework:Convolution]") {
	TimeSignal<float> u;
	TimeSignal<float> v;

	u.Resize(1000, 0.0f);
	for (int i = 0; i < 20; ++i) {
		u[200 + i] = float(i);
	}

	v.Resize(51, 0.0f);
	v[0] = 1;
	v[50] = 1;

	auto result = Convolution(u, v, convolution::central);

	REQUIRE(result.Length() == 950);
	REQUIRE(result[145] == 0);
	REQUIRE(result[151] == 1);
	REQUIRE(result[169] == 19);
	REQUIRE(result[185] == 0);
	REQUIRE(result[201] == 1);
	REQUIRE(result[219] == 19);
	REQUIRE(result[225] == 0);
}