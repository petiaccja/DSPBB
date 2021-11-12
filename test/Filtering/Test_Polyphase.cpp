#include "dspbb/Generators/Waveforms.hpp"

#include <catch2/catch.hpp>
#include <dspbb/Filtering/Polyphase.hpp>

using namespace dspbb;


TEST_CASE("Polyphase view filter non-uniform", "[Polyphase]") {
	const TimeSignal<float> filter = { 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2 };
	const std::array<size_t, 4> filterSizes = { 3, 3, 3, 2 };
	TimeSignal<float> scratch(filter.Size(), std::numeric_limits<float>::quiet_NaN());
	const auto view = PolyphaseDecompose(scratch, filter, 4);
	REQUIRE(view.numFilters == 4);
	for (size_t i = 0; i < 4; ++i) {
		REQUIRE(std::all_of(view[i].begin(), view[i].end(), [i](float c) { return c == float(4 * i); }));
		REQUIRE(view[i].Size() == filterSizes[i]);
	}
}

TEST_CASE("Polyphase view filter uniform", "[Polyphase]") {
	const TimeSignal<float> filter = { 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3 };
	const std::array<size_t, 4> filterSizes = { 3, 3, 3, 3 };
	TimeSignal<float> scratch(filter.Size(), std::numeric_limits<float>::quiet_NaN());
	const auto view = PolyphaseDecompose(scratch, filter, 4);
	REQUIRE(view.numFilters == 4);
	for (size_t i = 0; i < 4; ++i) {
		REQUIRE(std::all_of(view[i].begin(), view[i].end(), [i](float c) { return c == float(4 * i); }));
		REQUIRE(view[i].Size() == filterSizes[i]);
	}
}

TEST_CASE("Polyphase normalize", "[Polyphase]") {
	const TimeSignal<float> filter = { 1, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3 };
	TimeSignal<float> scratch(filter.Size(), std::numeric_limits<float>::quiet_NaN());
	const auto view = PolyphaseNormalized(PolyphaseDecompose(scratch, filter, 4));
	for (size_t i = 0; i < 4; ++i) {
		REQUIRE(Sum(view[i]) == Approx(1.0f));
	}
}


TEST_CASE("Polyphase reverse", "[Polyphase]") {
	const TimeSignal<float> filter = { 0, 1, 2, 3 };
	TimeSignal<float> scratch(filter.Size(), std::numeric_limits<float>::quiet_NaN());
	const auto view = PolyphaseDecompose(scratch, filter, 2);
	REQUIRE(view[0][0] == 2 * 2);
	REQUIRE(view[0][1] == 2 * 0);
	REQUIRE(view[1][0] == 2 * 3);
	REQUIRE(view[1][1] == 2 * 1);
}
