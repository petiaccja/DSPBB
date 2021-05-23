#include <catch2/catch.hpp>
#include <cmath>
#include <dspbb/Filtering/Interpolation.hpp>

using namespace dspbb;

auto MakeRamp(size_t size) {
	TimeSignal<float> signal;
	for (size_t i = 0; i < 150; ++i) {
		signal.PushBack(float(i));
	}
	return signal;
}


TEST_CASE("Polyphase interpolation replicate convolution full", "[Interpolation]") {
	constexpr int numFilters = 4;
	const TimeSignal<float> filter = FirLowPassWindowed<float>(0.5f / numFilters, 1, 31, windows::hamming);
	TimeSignal<float> scratch(filter.Size());
	const auto polyphase = PolyphaseDecompose(scratch, filter, numFilters);

	const TimeSignal<float> signal = MakeRamp(150);
	TimeSignal<float> padded(signal.Size() * numFilters);
	InterpolateZeroFill(AsView(padded), AsConstView(signal), numFilters);

	TimeSignal<float> outputConv = Convolution(padded, filter * numFilters, convolution::full);
	TimeSignal<float> output(outputConv.Size());

	Interpolate(AsView(output), AsConstView(signal), polyphase, { 1, numFilters });

	REQUIRE(Max(Abs(output - outputConv)) < 0.0001f);
}


TEST_CASE("Polyphase interpolation upsample constant", "[Interpolation]") {
	constexpr int numFilters = 4;
	const TimeSignal<float> filter = FirLowPassWindowed<float>(0.5f / numFilters, 1, 31, windows::hamming);
	TimeSignal<float> scratch(filter.Size());
	const auto polyphase = PolyphaseNormalized(PolyphaseDecompose(scratch, filter, numFilters));

	TimeSignal<float> signal(150, 1.0f);
	TimeSignal<float> output(signal.Size() - polyphase[0].Size() - 1);

	Interpolate(AsView(output), AsConstView(signal), polyphase, { 7, 11 }, { filter.Size() * 100 + 62 * numFilters, numFilters * 100 });

	REQUIRE(Min(output) == Approx(1.f));
	REQUIRE(Max(output) == Approx(1.f));
}


TEST_CASE("Polyphase interpolation upsample ramp", "[Interpolation]") {
	constexpr int numFilters = 4;
	const TimeSignal<float> filter = FirLowPassWindowed<float>(0.5f / numFilters, 1, 31, windows::hamming);
	TimeSignal<float> scratch(filter.Size());
	const auto polyphase = PolyphaseNormalized(PolyphaseDecompose(scratch, filter, numFilters));

	const TimeSignal<float> signal = MakeRamp(150);
	TimeSignal<float> output(signal.Size() - polyphase[0].Size() - 1);

	Interpolate(AsView(output), AsConstView(signal), polyphase, { 7, 11 }, { filter.Size(), numFilters });

	const TimeSignalView<const float> frontView{ output.begin(), output.Size() - 1 };
	const TimeSignalView<const float> backView{ output.begin() + 1, output.Size() - 1 };
	const auto diff = backView - frontView;
	REQUIRE(Mean(diff) == Approx(7.f / 11.f).epsilon(0.001f));
	REQUIRE(Min(diff) == Approx(7.f / 11.f).epsilon(0.02f));
	REQUIRE(Max(diff) == Approx(7.f / 11.f).epsilon(0.02f));
}


TEST_CASE("Polyphase interpolation downsample ramp mild", "[Interpolation]") {
	constexpr int numFilters = 4;
	const std::pair<uint64_t, uint64_t> ratio = { 11, 7 };
	const float ratioReal = float(ratio.first) / float(ratio.second);
	const TimeSignal<float> filter = FirLowPassWindowed<float>(0.5f / ratioReal / numFilters, 1, 31, windows::hamming);

	TimeSignal<float> scratch(filter.Size());
	const auto polyphase = PolyphaseNormalized(PolyphaseDecompose(scratch, filter, numFilters));

	const TimeSignal<float> signal = MakeRamp(150);
	TimeSignal<float> output(signal.Size() * ratio.second / ratio.first - polyphase[0].Size() - 1);

	Interpolate(AsView(output), AsConstView(signal), polyphase, ratio, { filter.Size(), numFilters });

	const TimeSignalView<const float> frontView{ output.begin(), output.Size() - 1 };
	const TimeSignalView<const float> backView{ output.begin() + 1, output.Size() - 1 };
	const auto diff = backView - frontView;
	REQUIRE(Mean(diff) == Approx(ratioReal).epsilon(0.001f));
	REQUIRE(Min(diff) == Approx(ratioReal).epsilon(0.02f));
	REQUIRE(Max(diff) == Approx(ratioReal).epsilon(0.02f));
}


TEST_CASE("Polyphase interpolation downsample ramp strong", "[Interpolation]") {
	constexpr int numFilters = 4;
	const std::pair<uint64_t, uint64_t> ratio = { 39, 7 };
	const float ratioReal = float(ratio.first) / float(ratio.second);
	const TimeSignal<float> filter = FirLowPassWindowed<float>(0.5f / ratioReal / numFilters, 1, 31, windows::hamming);

	TimeSignal<float> scratch(filter.Size());
	const auto polyphase = PolyphaseNormalized(PolyphaseDecompose(scratch, filter, numFilters));

	const TimeSignal<float> signal = MakeRamp(150);
	TimeSignal<float> output(signal.Size() * ratio.second / ratio.first - polyphase[0].Size() - 1);

	Interpolate(AsView(output), AsConstView(signal), polyphase, ratio, { filter.Size(), numFilters });

	const TimeSignalView<const float> frontView{ output.begin(), output.Size() - 1 };
	const TimeSignalView<const float> backView{ output.begin() + 1, output.Size() - 1 };
	const auto diff = backView - frontView;
	REQUIRE(Mean(diff) == Approx(ratioReal).epsilon(0.001f));
	REQUIRE(Min(diff) == Approx(ratioReal).epsilon(0.02f));
	REQUIRE(Max(diff) == Approx(ratioReal).epsilon(0.02f));
}


TEST_CASE("Polyphase interpolation shift ramp", "[Interpolation]") {
	constexpr int numFilters = 2;
	const TimeSignal<float> filter = FirLowPassWindowed<float>(0.5f / numFilters, 1, 63, windows::hamming);
	TimeSignal<float> scratch(filter.Size());
	const auto polyphase = PolyphaseNormalized(PolyphaseDecompose(scratch, filter, numFilters));

	const TimeSignal<float> signal = MakeRamp(150);
	TimeSignal<float> output(signal.Size() - polyphase[0].Size() - 1);

	const std::pair<uint64_t, uint64_t> offset = { filter.Size() * 100 + 42 * numFilters, numFilters * 100 };
	const float offsetReal = float(offset.first) / float(offset.second) - float(filter.Size() / numFilters) / 2;
	Interpolate(AsView(output), AsConstView(signal), polyphase, { 1, 1 }, offset);

	REQUIRE(output[0] == Approx(offsetReal).epsilon(0.02f));
}


TEST_CASE("Polyphase interpolation returned offset", "[Interpolation]") {
	constexpr int numFilters = 5;
	const TimeSignal<float> filter = FirLowPassWindowed<float>(0.5f / numFilters, 1, 63, windows::hamming);
	TimeSignal<float> scratch(filter.Size());
	const auto polyphase = PolyphaseNormalized(PolyphaseDecompose(scratch, filter, numFilters));

	const TimeSignal<float> signal = MakeRamp(150);
	TimeSignal<float> output(17);

	const std::pair<int64_t, uint64_t> offset = { 173, 982 };
	const std::pair<int64_t, uint64_t> ratio = { 7743, 9235 };
	const auto last = Interpolate(AsView(output), AsConstView(signal), polyphase, ratio, offset);

	REQUIRE(last.first / last.second == 17 * ratio.first / ratio.second);
	REQUIRE((last.first % last.second) / double(last.second) == Approx(0.4296632));
}
