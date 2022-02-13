#include <dspbb/Kernels/Convolution.hpp>

#include <catch2/catch.hpp>

using namespace dspbb;


constexpr std::array<float, 20> ur = { 1, 3, 7, 2, 9, 2, 5, 3, 7, 2, 4, 7, 3, 6, 3, 9, 3, 5, 3, 5 };
constexpr std::array<float, 12> vr = { 4, 3, 5, 2, 6, 3, 2, 7, 8, 5, 3, 3 };
constexpr std::array<float, 9> urvr_central = { 227, 244, 238, 207, 270, 219, 242, 223, 259 };
constexpr std::array<float, 31> urvr_full = {
	4, 15, 42, 46, 89, 80, 128, 101, 169, 175, 205, 227, 244, 238, 207, 270,
	219, 242, 223, 259, 210, 205, 196, 184, 152, 122, 120, 79, 49, 24, 15
};


TEST_CASE("Convolution naive central", "[Kernels - Convolution]") {
	std::array<float, 9> out;
	kernels::ConvolutionNaive(ur.begin(), ur.end(), vr.begin(), vr.end(), out.begin(), out.end(), 11);
	REQUIRE(out == urvr_central);
}

TEST_CASE("Convolution naive full", "[Kernels - Convolution]") {
	std::array<float, 31> out;
	kernels::ConvolutionNaive(ur.begin(), ur.end(), vr.begin(), vr.end(), out.begin(), out.end(), 0);
	REQUIRE(out == urvr_full);
}


TEST_CASE("Convolution slide central", "[Kernels - Convolution]") {
	std::array<float, 9> out;
	kernels::ConvolutionSlide(ur.begin(), ur.end(), vr.begin(), vr.end(), out.begin(), out.end(), 11);
	REQUIRE(out == urvr_central);
}

TEST_CASE("Convolution slide full", "[Kernels - Convolution]") {
	std::array<float, 31> out;
	kernels::ConvolutionSlide(ur.begin(), ur.end(), vr.begin(), vr.end(), out.begin(), out.end(), 0);
	REQUIRE(out == urvr_full);
}

TEST_CASE("Convolution accumulate central", "[Kernels - Convolution]") {
	std::array<float, 9> out;
	kernels::ConvolutionReduce(ur.begin(), ur.end(), vr.begin(), vr.end(), out.begin(), out.end(), 11);
	REQUIRE(out == urvr_central);
}

TEST_CASE("Convolution accumulate full", "[Kernels - Convolution]") {
	std::array<float, 31> out;
	kernels::ConvolutionReduce(ur.begin(), ur.end(), vr.begin(), vr.end(), out.begin(), out.end(), 0);
	REQUIRE(out == urvr_full);
}

TEST_CASE("Convolution accumulate small filter", "[Kernels - Convolution]") {
	constexpr std::array<float, 13> u = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 8, 7, 6, 5 };
	constexpr std::array<float, 3> v = { 0.125f, 0.5f, 0.25f };
	std::array<float, 15> ref;
	std::array<float, 15> out;
	kernels::ConvolutionNaive(u.begin(), u.end(), v.begin(), v.end(), ref.begin(), ref.end(), 0);
	kernels::ConvolutionReduce(u.begin(), u.end(), v.begin(), v.end(), out.begin(), out.end(), 0);
	REQUIRE(out == ref);
}


TEST_CASE("Convolution acc_vec central", "[Kernels - Convolution]") {
	std::array<float, 9> out;
	kernels::ConvolutionReduceVec(ur.begin(), ur.end(), vr.begin(), vr.end(), out.begin(), out.end(), 11);
	REQUIRE(out == urvr_central);
}

TEST_CASE("Convolution acc_vec full", "[Kernels - Convolution]") {
	std::array<float, 31> out;
	kernels::ConvolutionReduceVec(ur.begin(), ur.end(), vr.begin(), vr.end(), out.begin(), out.end(), 0);
	REQUIRE(out == urvr_full);
}

TEST_CASE("Convolution acc_vec small filter", "[Kernels - Convolution]") {
	constexpr std::array<float, 13> u = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 8, 7, 6, 5 };
	constexpr std::array<float, 3> v = { 0.125f, 0.5f, 0.25f };
	std::array<float, 15> ref;
	std::array<float, 15> out;
	kernels::ConvolutionNaive(u.begin(), u.end(), v.begin(), v.end(), ref.begin(), ref.end(), 0);
	kernels::ConvolutionReduceVec(u.begin(), u.end(), v.begin(), v.end(), out.begin(), out.end(), 0);
	REQUIRE(out == ref);
}