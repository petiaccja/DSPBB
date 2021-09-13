#include <catch2/catch.hpp>
#include <cmath>
#include <dspbb/Filtering/Interpolation.hpp>

using namespace dspbb;

auto MakeRamp(size_t size) {
	TimeSignal<float> signal;
	for (size_t i = 0; i < size; ++i) {
		signal.PushBack(float(i));
	}
	return signal;
}

TEST_CASE("Decimate", "[Interpolation]") {
	const TimeSignal<float> s = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
	const TimeSignal<float> d = Decimate(s, 3);
	REQUIRE(d.Size() == 4);
	REQUIRE(d[0] == 1);
	REQUIRE(d[1] == 4);
	REQUIRE(d[2] == 7);
	REQUIRE(d[3] == 10);
}


TEST_CASE("Expand", "[Interpolation]") {
	const TimeSignal<float> s = { 1, 2, 3 };
	const TimeSignal<float> e = Expand(s, 3);
	const TimeSignal<float> exp = { 1, 0, 0, 2, 0, 0, 3, 0, 0 };

	REQUIRE(e.Size() == 9);
	REQUIRE(Max(Abs(e - exp)) == Approx(0.0f));
}


TEST_CASE("Polyphase interpolation", "[Interpolation]") {
	constexpr int numFilters = 3;
	const TimeSignal<float> filter = FirFilter<float, TIME_DOMAIN>(31, Lowpass(WINDOWED).Cutoff(1.0f / numFilters));
	TimeSignal<float> scratch(filter.Size());
	const auto polyphase = PolyphaseNormalized(PolyphaseDecompose(scratch, filter, numFilters));

	const TimeSignal<float> signal = MakeRamp(150);
	TimeSignal<float> output(numFilters * ConvolutionLength(signal.Size(), polyphase.Size(), convolution::full));
	Interpolate(output, signal, polyphase, 0);

	REQUIRE(output[0] == Approx(0).margin(0.001f));

	const TimeSignalView<const float> frontView{ output.begin() + 100, output.begin() + 350 };
	const TimeSignalView<const float> backView{ output.begin() + 101, output.begin() + 351 };
	const auto diff = backView - frontView;
	REQUIRE(Mean(diff) == Approx(1.0f / numFilters).epsilon(0.001f));
	REQUIRE(Min(diff) == Approx(1.0f / numFilters).epsilon(0.02f));
	REQUIRE(Max(diff) == Approx(1.0f / numFilters).epsilon(0.02f));
}


TEST_CASE("Polyphase resampling replicate convolution full", "[Interpolation]") {
	constexpr int numFilters = 4;
	const TimeSignal<float> filter = FirFilter<float, TIME_DOMAIN>(31, Lowpass(WINDOWED).Cutoff(1.0f / numFilters));
	TimeSignal<float> scratch(filter.Size());
	const auto polyphase = PolyphaseDecompose(scratch, filter, numFilters);

	const TimeSignal<float> signal = MakeRamp(150);
	TimeSignal<float> padded(signal.Size() * numFilters);
	Expand(padded, signal, numFilters);

	TimeSignal<float> outputConv = Convolution(padded, filter * numFilters, convolution::full);
	TimeSignal<float> output(outputConv.Size());

	Resample(output, signal, polyphase, { 1, numFilters });

	REQUIRE(Max(Abs(output - outputConv)) < 0.0001f);
}


TEST_CASE("Polyphase resampling upsample constant", "[Interpolation]") {
	constexpr int numFilters = 4;
	const TimeSignal<float> filter = FirFilter<float, TIME_DOMAIN>(31, Lowpass(WINDOWED).Cutoff(1.0f / numFilters));
	TimeSignal<float> scratch(filter.Size());
	const auto polyphase = PolyphaseNormalized(PolyphaseDecompose(scratch, filter, numFilters));

	TimeSignal<float> signal(150, 1.0f);
	TimeSignal<float> output(signal.Size() - polyphase[0].Size() - 1);

	Resample(output, signal, polyphase, { 7, 11 }, { filter.Size() * 100 + 62 * numFilters, numFilters * 100 });

	REQUIRE(Min(output) == Approx(1.f));
	REQUIRE(Max(output) == Approx(1.f));
}


TEST_CASE("Polyphase resampling upsample ramp", "[Interpolation]") {
	constexpr int numFilters = 4;
	const TimeSignal<float> filter = FirFilter<float, TIME_DOMAIN>(31, Lowpass(WINDOWED).Cutoff(1.0f / numFilters));
	TimeSignal<float> scratch(filter.Size());
	const auto polyphase = PolyphaseNormalized(PolyphaseDecompose(scratch, filter, numFilters));

	const TimeSignal<float> signal = MakeRamp(150);
	TimeSignal<float> output(signal.Size() - polyphase[0].Size() - 1);

	Resample(output, signal, polyphase, { 7, 11 }, { filter.Size(), numFilters });

	const TimeSignalView<const float> frontView{ output.begin(), output.Size() - 1 };
	const TimeSignalView<const float> backView{ output.begin() + 1, output.Size() - 1 };
	const auto diff = backView - frontView;
	REQUIRE(Mean(diff) == Approx(7.f / 11.f).epsilon(0.001f));
	REQUIRE(Min(diff) == Approx(7.f / 11.f).epsilon(0.02f));
	REQUIRE(Max(diff) == Approx(7.f / 11.f).epsilon(0.02f));
}


TEST_CASE("Polyphase resampling downsample ramp mild", "[Interpolation]") {
	constexpr int numFilters = 4;
	const std::pair<uint64_t, uint64_t> ratio = { 11, 7 };
	const float ratioReal = float(ratio.first) / float(ratio.second);
	const TimeSignal<float> filter = FirFilter<float, TIME_DOMAIN>(31, Lowpass(WINDOWED).Cutoff(1.0f / ratioReal / numFilters));

	TimeSignal<float> scratch(filter.Size());
	const auto polyphase = PolyphaseNormalized(PolyphaseDecompose(scratch, filter, numFilters));

	const TimeSignal<float> signal = MakeRamp(150);
	TimeSignal<float> output(signal.Size() * ratio.second / ratio.first - polyphase[0].Size() - 1);

	Resample(output, signal, polyphase, ratio, { filter.Size(), numFilters });

	const TimeSignalView<const float> frontView{ output.begin(), output.Size() - 1 };
	const TimeSignalView<const float> backView{ output.begin() + 1, output.Size() - 1 };
	const auto diff = backView - frontView;
	REQUIRE(Mean(diff) == Approx(ratioReal).epsilon(0.001f));
	REQUIRE(Min(diff) == Approx(ratioReal).epsilon(0.02f));
	REQUIRE(Max(diff) == Approx(ratioReal).epsilon(0.02f));
}


TEST_CASE("Polyphase resampling downsample ramp strong", "[Interpolation]") {
	constexpr int numFilters = 4;
	const std::pair<uint64_t, uint64_t> ratio = { 39, 7 };
	const float ratioReal = float(ratio.first) / float(ratio.second);
	const TimeSignal<float> filter = FirFilter<float, TIME_DOMAIN>(31, Lowpass(WINDOWED).Cutoff(1.0f / ratioReal / numFilters));

	TimeSignal<float> scratch(filter.Size());
	const auto polyphase = PolyphaseNormalized(PolyphaseDecompose(scratch, filter, numFilters));

	const TimeSignal<float> signal = MakeRamp(150);
	TimeSignal<float> output(signal.Size() * ratio.second / ratio.first - polyphase[0].Size() - 1);

	Resample(output, signal, polyphase, ratio, { filter.Size(), numFilters });

	const TimeSignalView<const float> frontView{ output.begin(), output.Size() - 1 };
	const TimeSignalView<const float> backView{ output.begin() + 1, output.Size() - 1 };
	const auto diff = backView - frontView;
	REQUIRE(Mean(diff) == Approx(ratioReal).epsilon(0.001f));
	REQUIRE(Min(diff) == Approx(ratioReal).epsilon(0.02f));
	REQUIRE(Max(diff) == Approx(ratioReal).epsilon(0.02f));
}


TEST_CASE("Polyphase resampling shift ramp", "[Interpolation]") {
	constexpr int numFilters = 2;
	const TimeSignal<float> filter = FirFilter<float, TIME_DOMAIN>(63, Lowpass(WINDOWED).Cutoff(1.0f / numFilters));
	TimeSignal<float> scratch(filter.Size());
	const auto polyphase = PolyphaseNormalized(PolyphaseDecompose(scratch, filter, numFilters));

	const TimeSignal<float> signal = MakeRamp(150);
	TimeSignal<float> output(signal.Size() - polyphase[0].Size() - 1);

	const std::pair<uint64_t, uint64_t> offset = { filter.Size() * 100 + 42 * numFilters, numFilters * 100 };
	const float offsetReal = float(offset.first) / float(offset.second) - float(filter.Size() / numFilters) / 2;
	Resample(output, signal, polyphase, { 1, 1 }, offset);

	REQUIRE(output[0] == Approx(offsetReal).epsilon(0.02f));
}


TEST_CASE("Polyphase resampling returned offset", "[Interpolation]") {
	constexpr int numFilters = 5;
	const TimeSignal<float> filter = FirFilter<float, TIME_DOMAIN>(63, Lowpass(WINDOWED).Cutoff(1.0f / numFilters));
	TimeSignal<float> scratch(filter.Size());
	const auto polyphase = PolyphaseNormalized(PolyphaseDecompose(scratch, filter, numFilters));

	const TimeSignal<float> signal = MakeRamp(150);
	TimeSignal<float> output(17);

	const std::pair<int64_t, uint64_t> offset = { 173, 982 };
	const std::pair<int64_t, uint64_t> ratio = { 7743, 9235 };
	const auto last = Resample(output, signal, polyphase, ratio, offset);

	REQUIRE(last.first / last.second == 17 * ratio.first / ratio.second);
	REQUIRE((last.first % last.second) / double(last.second) == Approx(0.4296632));
}
