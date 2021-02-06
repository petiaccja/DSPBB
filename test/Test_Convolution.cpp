#include <dspbb/Filtering/Convolution.hpp>

#include <Catch2/catch.hpp>

using namespace dspbb;
using namespace std::complex_literals;

TEST_CASE("Real central", "[AudioFramework:Convolution]") {
	TimeSignal<float> u = { 1, 2, 3 };
	TimeSignal<float> v = { 1, 1 };
	TimeSignal<float> expected = { 3, 5 };

	auto result = ConvolutionFast(u, v, convolution::central);

	REQUIRE(result.Length() == 2);
	for (size_t i = 0; i < 2; ++i) {
		REQUIRE(result[i] == expected[i]);
	}
}

TEST_CASE("Real full", "[AudioFramework:Convolution]") {
	TimeSignal<float> u = { 1, 2, 3 };
	TimeSignal<float> v = { 1, 1 };
	TimeSignal<float> expected = { 1, 3, 5, 3 };

	auto result = ConvolutionFast(u, v, convolution::full);

	REQUIRE(result.Length() == 4);
	for (size_t i = 0; i < 4; ++i) {
		REQUIRE(result[i] == expected[i]);
	}
}


TEST_CASE("Real-complex central BR", "[AudioFramework:Convolution]") {
	TimeSignal<float> u = { 1, 2, 3 };
	TimeSignal<std::complex<float>> v = { 1.f + 2.if, 1.f + 3.if };
	TimeSignal<std::complex<float>> expected = { 3.f + 8.if, 5.f + 13.if };

	auto result = ConvolutionFast(u, v, convolution::central);

	REQUIRE(result.Length() == 2);
	for (size_t i = 0; i < 2; ++i) {
		REQUIRE(result[i] == expected[i]);
	}
}

TEST_CASE("Complex-real central RB", "[AudioFramework:Convolution]") {
	TimeSignal<std::complex<float>> u = { 1.f + 2.if, 1.f + 3.if };
	TimeSignal<float> v = { 1, 2, 3 };
	TimeSignal<std::complex<float>> expected = { 3.f + 8.if, 5.f + 13.if };

	auto result = ConvolutionFast(u, v, convolution::central);

	REQUIRE(result.Length() == 2);
	for (size_t i = 0; i < 2; ++i) {
		REQUIRE(result[i] == expected[i]);
	}
}

TEST_CASE("Complex-complex central", "[AudioFramework:Convolution]") {
	TimeSignal<std::complex<float>> u = { 1.f + 3.if, 2.f + 2.if, 3.f + 1.if };
	TimeSignal<std::complex<float>> v = { 1.f + 2.if, 1.f + 2.if };
	TimeSignal<std::complex<float>> expected = { -7.f + 11.if, -1.f + 13.if };

	auto result = ConvolutionFast(u, v, convolution::central);

	REQUIRE(result.Length() == 2);
	for (size_t i = 0; i < 2; ++i) {
		REQUIRE(result[i] == expected[i]);
	}
}


TEST_CASE("Different types", "[AudioFramework:Convolution]") {
	TimeSignal<std::complex<float>> u = { 1.f + 3.if, 2.f + 2.if, 3.f + 1.if };
	TimeSignal<std::complex<double>> v = { 1. + 2.i, 1. + 2.i };
	TimeSignal<std::complex<double>> expected = { -7. + 11.i, -1. + 13.i };

	auto result = ConvolutionFast(u, v, convolution::central);

	REQUIRE(result.Length() == 2);
	for (size_t i = 0; i < 2; ++i) {
		REQUIRE(result[i] == expected[i]);
	}
}


TEST_CASE("Real-world signal", "[AudioFramework:Convolution]") {
	TimeSignal<float> u;
	TimeSignal<float> v;

	u.Resize(1000, 0.0f);
	for (int i=0; i<20; ++i) {
		u[200+i] = float(i);
	}

	v.Resize(51, 0.0f);
	v[0] = 1;
	v[50] = 1;

	auto result = ConvolutionFast(u, v, convolution::central);

	REQUIRE(result.Length() == 950);
	REQUIRE(result[145] == 0);
	REQUIRE(result[151] == 1);
	REQUIRE(result[169] == 19);
	REQUIRE(result[185] == 0);
	REQUIRE(result[201] == 1);
	REQUIRE(result[219] == 19);
	REQUIRE(result[225] == 0);	
}