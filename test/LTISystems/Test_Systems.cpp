#include "../TestUtils.hpp"

#include <dspbb/LTISystems/Systems.hpp>

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>


using namespace dspbb;
using namespace std::complex_literals;
using Catch::Approx;


template <class T>
bool ContainsNumeratorConstants(const CascadedBiquad<T>& biquad, std::vector<T> expected) {
	return std::all_of(expected.begin(), expected.end(), [&biquad](const auto& e) {
		const bool secondOrder = std::any_of(biquad.sections.begin(), biquad.sections.end(), [&e](const auto& s) {
			return s.numerator[0] == Approx(e);
		});
		const bool firstOrder = std::any_of(biquad.sections.begin(), biquad.sections.end(), [&e](const auto& s) {
			return s.numerator[0] == T(0) && s.numerator[1] == Approx(e);
		});
		return firstOrder || secondOrder;
	});
}


TEST_CASE("Biquad cascade conversion real pairing even-even", "[Biquad cascade]") {
	const DiscreteZeroPoleGain<float> sys{
		1.0f,
		{
			-2.0f,
			-4.0f,
			-1.0f,
			-6.0f,
			3.0f,
			5.0f,
			9.0f,
			1.0f,
		}
	};
	CascadedBiquad cascade{ sys };
	REQUIRE(cascade.sections.size() == 4);
	REQUIRE(ContainsNumeratorConstants(cascade, { 6, 8, 9, 15 }));
}

TEST_CASE("Biquad cascade conversion real pairing odd-even", "[Biquad cascade]") {
	const DiscreteZeroPoleGain<float> sys{
		1.0f,
		{
			-2.0f,
			-4.0f,
			-6.0f,
			3.0f,
			5.0f,
			9.0f,
			1.0f,
		}
	};
	CascadedBiquad cascade{ sys };
	REQUIRE(cascade.sections.size() == 4);
	REQUIRE(ContainsNumeratorConstants(cascade, { 4.0f, 9.0f, 12.0f, 15.0f }));
}

TEST_CASE("Biquad cascade conversion real pairing even-odd", "[Biquad cascade]") {
	const DiscreteZeroPoleGain<float> sys{
		1.0f,
		{
			-2.0f,
			-4.0f,
			-1.0f,
			-6.0f,
			3.0f,
			5.0f,
			1.0f,
		}
	};
	CascadedBiquad cascade{ sys };
	REQUIRE(cascade.sections.size() == 4);
	REQUIRE(ContainsNumeratorConstants(cascade, { -3.0f, 5.0f, 6.0f, 8.0f }));
}

TEST_CASE("Biquad cascade conversion real pairing odd-odd", "[Biquad cascade]") {
	const DiscreteZeroPoleGain<float> sys{
		1.0f,
		{
			-2.0f,
			-4.0f,
			-1.0f,
			3.0f,
			9.0f,
			1.0f,
		}
	};
	CascadedBiquad cascade{ sys };
	REQUIRE(cascade.sections.size() == 3);
	REQUIRE(ContainsNumeratorConstants(cascade, { 4, -6, 9 }));
}


TEST_CASE("Biquad cascade conversion real pairing empty-odd", "[Biquad cascade]") {
	const DiscreteZeroPoleGain<float> sys{
		1.0f,
		{ 3, 9, 1 }
	};
	CascadedBiquad cascade{ sys };
	REQUIRE(cascade.sections.size() == 2);
	REQUIRE(ContainsNumeratorConstants(cascade, { 9, -3 }));
}

TEST_CASE("Biquad cascade conversion real pairing even-empty", "[Biquad cascade]") {
	const DiscreteZeroPoleGain<float> sys{
		1.0f,
		{ -2, -5, -2, -8 }
	};
	CascadedBiquad cascade{ sys };
	REQUIRE(cascade.sections.size() == 2);
	REQUIRE(ContainsNumeratorConstants(cascade, { 16, 10 }));
}


TEST_CASE("Biquad cascade conversion complex", "[Biquad cascade]") {
	const DiscreteZeroPoleGain<float> sys{
		1.0f,
		{ 1.0f + 2.0if, 1.0f - 2.0if, 3.0f + 4.0if, 3.0f - 4.0if }
	};
	CascadedBiquad cascade{ sys };
	REQUIRE(cascade.sections.size() == 2);
	REQUIRE(ContainsNumeratorConstants(cascade, { 5.0f, 25.0f }));
}

TEST_CASE("Biquad cascade gain", "[Biquad cascade]") {
	const DiscreteZeroPoleGain<float> sys{
		2.718f,
		{ 1.0f + 2.0if, 1.0f - 2.0if, 3.0f + 4.0if, 3.0f - 4.0if },
		{ 4.0f + 2.0if, 4.0f - 2.0if, 2.0f + 4.0if, 2.0f - 4.0if }
	};
	CascadedBiquad cascade{ sys };
	const float gain = std::transform_reduce(cascade.sections.begin(), cascade.sections.end(), 1.0f, std::multiplies{}, [](const auto& section) {
		return section.numerator[2];
	});
	REQUIRE(gain == Approx(2.718f));
}

const std::array complexPoints = { 1.345f + 0.928if, 0.7823f + 2.3778if };
const std::array realPoints = { 1.345f, 0.7823f };

TEST_CASE("Biquad cascade equation evaluation 2x first order", "[Biquad cascade]") {
	const DiscreteZeroPoleGain<float> sys{
		6.67f,
		{ 3.0f, 1.0f + 2.0if, 1.0f - 2.0if, 3.0f + 4.0if, 3.0f - 4.0if },
		{ 2.0f, 4.0f + 2.0if, 4.0f - 2.0if, 2.0f + 4.0if, 2.0f - 4.0if }
	};
	CascadedBiquad cascade{ sys };
	for (auto& ri : realPoints) {
		REQUIRE(sys(ri) == Approx(cascade(ri)));
	}
	for (auto& ci : complexPoints) {
		REQUIRE(sys(ci) == ApproxComplex(cascade(ci)));
	}
}

TEST_CASE("Biquad cascade equation evaluation 1x num first order", "[Biquad cascade]") {
	const DiscreteZeroPoleGain<float> sys{
		6.67f,
		{ 3.0f, 1.0f + 2.0if, 1.0f - 2.0if, 3.0f + 4.0if, 3.0f - 4.0if },
		{ 4.0f + 2.0if, 4.0f - 2.0if, 2.0f + 4.0if, 2.0f - 4.0if }
	};
	CascadedBiquad cascade{ sys };
	for (auto& ri : realPoints) {
		REQUIRE(sys(ri) == Approx(cascade(ri)));
	}
	for (auto& ci : complexPoints) {
		REQUIRE(sys(ci) == ApproxComplex(cascade(ci)));
	}
}

TEST_CASE("Biquad cascade equation evaluation 1x den first order", "[Biquad cascade]") {
	const DiscreteZeroPoleGain<float> sys{
		6.67f,
		{ 5.0f, 2.0f, 3.0f, 1.0f + 2.0if, 1.0f - 2.0if, 3.0f + 4.0if, 3.0f - 4.0if },
		{ 4.0f + 2.0if, 4.0f - 2.0if, 2.0f + 4.0if, 2.0f - 4.0if }
	};
	CascadedBiquad cascade{ sys };
	for (auto& ri : realPoints) {
		REQUIRE(sys(ri) == Approx(cascade(ri)));
	}
	for (auto& ci : complexPoints) {
		REQUIRE(sys(ci) == ApproxComplex(cascade(ci)));
	}
}

TEST_CASE("Biquad cascade equation evaluation no first order", "[Biquad cascade]") {
	const DiscreteZeroPoleGain<float> sys{
		6.67f,
		{ 1.0f + 2.0if, 1.0f - 2.0if, 3.0f + 4.0if, 3.0f - 4.0if },
		{ 4.0f + 2.0if, 4.0f - 2.0if, 2.0f + 4.0if, 2.0f - 4.0if }
	};
	CascadedBiquad cascade{ sys };
	for (auto& ri : realPoints) {
		REQUIRE(sys(ri) == Approx(cascade(ri)));
	}
	for (auto& ci : complexPoints) {
		REQUIRE(sys(ci) == ApproxComplex(cascade(ci)));
	}
}
