#include "dspbb/Generators/Waveforms.hpp"
#include <dspbb/Filtering/Polyphase.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>


using namespace dspbb;
using Catch::Approx;


TEST_CASE("Polyphase view filter non-uniform", "[Polyphase]") {
	const Signal<float> filter = { 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2 };
	const std::array<size_t, 4> filterSizes = { 3, 3, 3, 2 };
	const auto view = PolyphaseDecompose(filter, 4);
	REQUIRE(view.num_phases() == 4);
	for (size_t i = 0; i < 4; ++i) {
		REQUIRE(std::all_of(view[i].begin(), view[i].end(), [i](float c) { return c == float(4 * i); }));
		REQUIRE(view[i].size() == filterSizes[i]);
	}
}

TEST_CASE("Polyphase view filter uniform", "[Polyphase]") {
	const Signal<float> filter = { 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3 };
	const std::array<size_t, 4> filterSizes = { 3, 3, 3, 3 };
	const auto view = PolyphaseDecompose(filter, 4);
	REQUIRE(view.num_phases() == 4);
	for (size_t i = 0; i < 4; ++i) {
		REQUIRE(std::all_of(view[i].begin(), view[i].end(), [i](float c) { return c == float(4 * i); }));
		REQUIRE(view[i].size() == filterSizes[i]);
	}
}

TEST_CASE("Polyphase normalize", "[Polyphase]") {
	const Signal<float> filter = { 1, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3 };
	const auto view = PolyphaseNormalized(PolyphaseDecompose(filter, 4));
	for (size_t i = 0; i < 4; ++i) {
		REQUIRE(Sum(view[i]) == Approx(1.0f));
	}
}


TEST_CASE("Polyphase reverse", "[Polyphase]") {
	const Signal<float> filter = { 0, 1, 2, 3 };
	const auto view = PolyphaseDecompose(filter, 2);
	REQUIRE(view[0][0] == 2 * 2);
	REQUIRE(view[0][1] == 2 * 0);
	REQUIRE(view[1][0] == 2 * 3);
	REQUIRE(view[1][1] == 2 * 1);
}
