#include "../TestUtils.hpp"

#include <dspbb/Filtering/Resample.hpp>
#include <dspbb/Math/Convolution.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>


using namespace dspbb;
using Catch::Approx;


template <signal_or_view SignalT>
auto InterpolateRefImpl(const SignalT& signal, const SignalT& filter, size_t rate, size_t offset, size_t length) {
	return Convolution(Expand(signal, rate), filter, offset, length) * rate;
}


template <class T>
auto SmoothTestSignal(size_t length, T center, T scale) {
	const auto func = [&](size_t idx) {
		const auto signedIdx = double(idx) - double(center);
		const auto x = double(signedIdx) * 6.0 / scale;
		return static_cast<T>(2.0 * x * std::exp(-x * x));
	};
	Signal<T> signal(length);
	for (size_t i = 0; i < signal.size(); ++i) {
		signal[i] = func(i);
	}
	return signal / Max(signal);
}


template <signal_or_view SignalT>
double SmoothTestSignalCenter(const SignalT& signal) {
	const auto minIt = std::ranges::min_element(signal);
	const auto maxIt = std::ranges::max_element(signal);

	const auto roi = SignalView<const scalar_type_t<SignalT>>(minIt, maxIt);
	const auto it = std::ranges::adjacent_find(roi, [](const auto& lhs, const auto& rhs) {
		return lhs * rhs <= 0;
	});

	if (it == roi.end()) {
		return -1.0;
	}

	const auto lhs = double(*it);
	const auto rhs = double(*(it + 1));
	const auto frac = std::abs(lhs) / (std::abs(lhs) + std::abs(rhs));
	const auto whole = double(it - roi.begin() + minIt - signal.begin());
	return whole + frac;
}


//------------------------------------------------------------------------------
// Public utilities
//------------------------------------------------------------------------------

TEST_CASE("Interpolation filter cutoff", "[Interpolation]") {
	REQUIRE(InterpolFilterCutoff(4) == Approx(0.25));
	REQUIRE(InterpolFilterCutoff(9) == Approx(0.1111111111));
}


TEST_CASE("Resampling filter cutoff", "[Interpolation]") {
	REQUIRE(ResampleFilterCutoff({ 4, 6 }, 5) == Approx(0.2));
	REQUIRE(ResampleFilterCutoff({ 6, 4 }, 5) == Approx(0.1333333333));
	REQUIRE(ResampleFilterCutoff({ 4, 71 }, 12) == Approx(0.0833333333));
	REQUIRE(ResampleFilterCutoff({ 40, 6 }, 12) == Approx(0.0125));
}


TEST_CASE("Interpolation length", "[Interpolation]") {
	SECTION("Full") {
		constexpr size_t signalSize = 2000;
		constexpr size_t filterSize = 1001;
		constexpr size_t factor = 5;

		const auto size = InterpolLength(signalSize, filterSize, factor, CONV_FULL);
		REQUIRE(size == 11000);
	}
	SECTION("Central") {
		constexpr size_t signalSize = 2000;
		constexpr size_t filterSize = 1001;
		constexpr size_t factor = 5;

		const auto size = InterpolLength(signalSize, filterSize, factor, CONV_CENTRAL);
		REQUIRE(size == 9000);
	}
}


TEST_CASE("Resampling length", "[Interpolation]") {
	SECTION("Full") {
		SECTION("Upsample exact") {
			constexpr Rational<int64_t> sampleRates = { 2, 3 };
			constexpr size_t signalSize = 2000;
			constexpr size_t filterSize = 1001;
			constexpr size_t numPhases = 5;

			constexpr auto size = ResampleLength(signalSize, filterSize, numPhases, sampleRates, CONV_FULL);
			REQUIRE(double(size) == Approx(16500.0 / 5).margin(0.01));
		}
		SECTION("Upsample inexact") {
			constexpr Rational<int64_t> sampleRates = { 3, 5 };
			constexpr size_t signalSize = 2000;
			constexpr size_t filterSize = 1001;
			constexpr size_t numPhases = 5;

			constexpr auto size = ResampleLength(signalSize, filterSize, numPhases, sampleRates, CONV_FULL);
			REQUIRE(double(size) == Approx(18333.333 / 5).margin(0.01));
		}
		SECTION("Downsample exact") {
			constexpr Rational<int64_t> sampleRates = { 11000, 3500 };
			constexpr size_t signalSize = 2000;
			constexpr size_t filterSize = 1001;
			constexpr size_t numPhases = 5;

			const auto size = ResampleLength(signalSize, filterSize, numPhases, sampleRates, CONV_FULL);
			REQUIRE(double(size) == Approx(3500.0 / 5).margin(0.01));
		}
		SECTION("Downsample inexact") {
			constexpr Rational<int64_t> sampleRates = { 22000, 7001 };
			constexpr size_t signalSize = 2000;
			constexpr size_t filterSize = 1001;
			constexpr size_t numPhases = 5;

			constexpr auto size = ResampleLength(signalSize, filterSize, numPhases, sampleRates, CONV_FULL);
			REQUIRE(double(size) == Approx(3500.5 / 5).margin(0.01));
		}
	}
	SECTION("Central") {
		SECTION("Upsample exact") {
			constexpr Rational<int64_t> sampleRates = { 9000, 14000 };
			constexpr size_t signalSize = 2000;
			constexpr size_t filterSize = 1001;
			constexpr size_t numPhases = 5;

			constexpr auto size = ResampleLength(signalSize, filterSize, numPhases, sampleRates, CONV_CENTRAL);
			REQUIRE(double(size) == Approx(14000.0 / 5).margin(0.01));
		}
		SECTION("Upsample inexact") {
			constexpr Rational<int64_t> sampleRates = { 27000, 14000 };
			constexpr size_t signalSize = 2000;
			constexpr size_t filterSize = 1001;
			constexpr size_t numPhases = 5;

			constexpr auto size = ResampleLength(signalSize, filterSize, numPhases, sampleRates, CONV_CENTRAL);
			REQUIRE(double(size) == Approx(4666.667 / 5).margin(0.01));
		}
		SECTION("Downsample exact") {
			constexpr Rational<int64_t> sampleRates = { 9000, 3500 };
			constexpr size_t signalSize = 2000;
			constexpr size_t filterSize = 1001;
			constexpr size_t numPhases = 5;

			const auto size = ResampleLength(signalSize, filterSize, numPhases, sampleRates, CONV_CENTRAL);
			REQUIRE(double(size) == Approx(3500.0 / 5).margin(0.01));
		}
		SECTION("Downsample inexact") {
			constexpr Rational<int64_t> sampleRates = { 18000, 7001 };
			constexpr size_t signalSize = 2000;
			constexpr size_t filterSize = 1001;
			constexpr size_t numPhases = 5;

			constexpr auto size = ResampleLength(signalSize, filterSize, numPhases, sampleRates, CONV_CENTRAL);
			REQUIRE(double(size) == Approx(3500.5 / 5).margin(0.01));
		}
	}
}


TEST_CASE("Interpolation offset", "[Interpolation]") {
	SECTION("Full") {
		constexpr size_t signalSize = 2000;
		constexpr size_t filterSize = 1001;
		constexpr size_t factor = 5;

		const auto offset = InterpolOffset(signalSize, filterSize, factor, CONV_FULL);
		REQUIRE(offset == 0);
	}
	SECTION("Central") {
		constexpr size_t signalSize = 2000;
		constexpr size_t filterSize = 1001;
		constexpr size_t factor = 5;

		const auto offset = InterpolOffset(signalSize, filterSize, factor, CONV_CENTRAL);
		REQUIRE(offset == 1000);
	}
}



TEST_CASE("Resampling offset", "[Interpolation]") {
	SECTION("Full") {
		SECTION("Upsample exact") {
			constexpr Rational<int64_t> sampleRates = { 2, 3 };
			constexpr size_t signalSize = 2000;
			constexpr size_t filterSize = 1001;
			constexpr size_t numPhases = 5;

			constexpr auto offset = ResampleOffset(signalSize, filterSize, numPhases, sampleRates, CONV_FULL);
			REQUIRE(double(offset) == 0.0);
		}
		SECTION("Upsample inexact") {
			constexpr Rational<int64_t> sampleRates = { 3, 5 };
			constexpr size_t signalSize = 2000;
			constexpr size_t filterSize = 1001;
			constexpr size_t numPhases = 5;

			constexpr auto offset = ResampleOffset(signalSize, filterSize, numPhases, sampleRates, CONV_FULL);
			REQUIRE(double(offset) == 0.0);
		}
		SECTION("Downsample exact") {
			constexpr Rational<int64_t> sampleRates = { 11000, 3500 };
			constexpr size_t signalSize = 2000;
			constexpr size_t filterSize = 1001;
			constexpr size_t numPhases = 5;

			const auto offset = ResampleOffset(signalSize, filterSize, numPhases, sampleRates, CONV_FULL);
			REQUIRE(double(offset) == 0.0);
		}
		SECTION("Downsample inexact") {
			constexpr Rational<int64_t> sampleRates = { 22000, 7001 };
			constexpr size_t signalSize = 2000;
			constexpr size_t filterSize = 1001;
			constexpr size_t numPhases = 5;

			constexpr auto offset = ResampleOffset(signalSize, filterSize, numPhases, sampleRates, CONV_FULL);
			REQUIRE(double(offset) == 0.0);
		}
	}
	SECTION("Central") {
		SECTION("Upsample exact") {
			constexpr Rational<int64_t> sampleRates = { 9000, 14000 };
			constexpr size_t signalSize = 2000;
			constexpr size_t filterSize = 1001;
			constexpr size_t numPhases = 5;

			constexpr auto offset = ResampleOffset(signalSize, filterSize, numPhases, sampleRates, CONV_CENTRAL);
			REQUIRE(double(offset) == Approx(14000.0 / (9000.0 * 5) * 1000.0).margin(0.01));
		}
		SECTION("Upsample inexact") {
			constexpr Rational<int64_t> sampleRates = { 27000, 14000 };
			constexpr size_t signalSize = 2000;
			constexpr size_t filterSize = 1001;
			constexpr size_t numPhases = 5;

			constexpr auto offset = ResampleOffset(signalSize, filterSize, numPhases, sampleRates, CONV_CENTRAL);
			REQUIRE(double(offset) == Approx(14000.0 / (27000.0 * 5) * 1000.0).margin(0.01));
		}
		SECTION("Downsample exact") {
			constexpr Rational<int64_t> sampleRates = { 9000, 3500 };
			constexpr size_t signalSize = 2000;
			constexpr size_t filterSize = 1001;
			constexpr size_t numPhases = 5;

			const auto offset = ResampleOffset(signalSize, filterSize, numPhases, sampleRates, CONV_CENTRAL);
			REQUIRE(double(offset) == Approx(3500.0 / (9000 * 5) * 1000.0).margin(0.01));
		}
		SECTION("Downsample inexact") {
			constexpr Rational<int64_t> sampleRates = { 18000, 7001 };
			constexpr size_t signalSize = 2000;
			constexpr size_t filterSize = 1001;
			constexpr size_t numPhases = 5;

			constexpr auto offset = ResampleOffset(signalSize, filterSize, numPhases, sampleRates, CONV_CENTRAL);
			REQUIRE(double(offset) == Approx(7001.0 / (18000.0 * 5) * 1000.0).margin(0.01));
		}
	}
}


//------------------------------------------------------------------------------
// Internal utilities
//------------------------------------------------------------------------------

TEST_CASE("Find interpolation suspension point", "[Interpolation]") {
	constexpr size_t numPhases = 6;
	constexpr size_t filterSize = 31;

	SECTION("Initial point") {
		constexpr size_t nextOutputSample = 0;
		const auto [inputIndex, startPoint] = impl::FindInterpolSuspensionPoint(nextOutputSample, filterSize, numPhases);

		REQUIRE(inputIndex == 0);
		REQUIRE(startPoint == 0);
	}
	SECTION("One off") {
		constexpr size_t nextOutputSample = 2;
		const auto [inputIndex, startPoint] = impl::FindInterpolSuspensionPoint(nextOutputSample, filterSize, numPhases);

		REQUIRE(inputIndex == 0);
		REQUIRE(startPoint == 2);
	}
	SECTION("Middle point") {
		constexpr size_t nextOutputSample = 36;
		const auto [inputIndex, startPoint] = impl::FindInterpolSuspensionPoint(nextOutputSample, filterSize, numPhases);

		REQUIRE(inputIndex == 1);
		REQUIRE(startPoint == 30);
	}
	SECTION("Far point") {
		constexpr size_t nextOutputSample = 158;
		const auto [inputIndex, startPoint] = impl::FindInterpolSuspensionPoint(nextOutputSample, filterSize, numPhases);

		REQUIRE(inputIndex == 21);
		REQUIRE(startPoint == 32);
	}
}


TEST_CASE("Find resampling suspension point", "[Interpolation]") {
	constexpr size_t numPhases = 6;
	constexpr size_t filterSize = 31;
	constexpr Rational<int64_t> sampleRates = { 4, 7 };

	SECTION("Initial point") {
		constexpr Rational<int64_t> nextOutputSample = { 0, 1 };
		const auto [inputIndex, startPoint] = impl::FindResampleSuspensionPoint(nextOutputSample, filterSize, numPhases, sampleRates);

		REQUIRE(inputIndex == 0);
		REQUIRE(double(startPoint) == Approx(0));
	}
	SECTION("One off") {
		constexpr Rational<int64_t> nextOutputSample = { 7, 7 };
		const auto [inputIndex, startPoint] = impl::FindResampleSuspensionPoint(nextOutputSample, filterSize, numPhases, sampleRates);

		REQUIRE(inputIndex == 0);
		REQUIRE(double(startPoint) == Approx(1));
	}
	SECTION("Middle point") {
		constexpr Rational<int64_t> nextOutputSample = { 6 * 7, 4 };
		const auto [inputIndex, startPoint] = impl::FindResampleSuspensionPoint(nextOutputSample, filterSize, numPhases, sampleRates);

		REQUIRE(inputIndex == 1);
		const double expectedTotalOffset = double(nextOutputSample);
		const double actualTotalOffset = double(inputIndex) / double(sampleRates) + double(startPoint);
		REQUIRE(expectedTotalOffset == Approx(actualTotalOffset));
	}
	SECTION("Far point") {
		constexpr Rational<int64_t> nextOutputSample = { 156, 1 };
		const auto [inputIndex, startPoint] = impl::FindResampleSuspensionPoint(nextOutputSample, filterSize, numPhases, sampleRates);

		REQUIRE(inputIndex == 84);
		const double expectedTotalOffset = double(nextOutputSample);
		const double actualTotalOffset = double(inputIndex) / double(sampleRates) + double(startPoint);
		REQUIRE(expectedTotalOffset == Approx(actualTotalOffset));
	}
}



TEST_CASE("Resampling: Change sample rate", "[Interpolation]") {
	constexpr int inputRate = 7;
	constexpr int outputRate = 17;

	constexpr Rational<int64_t> originalSample = { 28, 42 };

	SECTION("Regular") {
		constexpr auto newSample = impl::ChangeSampleRate(inputRate, outputRate, originalSample);

		const double inputIndexRealExpected = double(originalSample) / double(inputRate) * double(outputRate);

		REQUIRE(double(newSample) == Approx(inputIndexRealExpected));
	}
	SECTION("Simplify") {
		constexpr auto newSample = impl::ChangeSampleRate(inputRate, outputRate, originalSample);

		const double inputIndexRealExpected = double(originalSample) / double(inputRate) * double(outputRate);

		REQUIRE(double(newSample) == Approx(inputIndexRealExpected));
	}
}


TEST_CASE("Resampling: Input index 2 samples", "[Interpolation]") {
	SECTION("Zero weight") {
		const auto [firstSample, secondSample] = impl::InputIndex2Sample({ 43, 7 }, 7);
		REQUIRE(firstSample.inputIndex == 6);
		REQUIRE(firstSample.phaseIndex == 1);
		REQUIRE(firstSample.weight == 1);

		REQUIRE(secondSample.inputIndex == 6);
		REQUIRE(secondSample.phaseIndex == 2);
		REQUIRE(secondSample.weight == 0);
	}
	SECTION("Split weight") {
		const auto [firstSample, secondSample] = impl::InputIndex2Sample({ 87, 14 }, 5);
		REQUIRE(firstSample.inputIndex == 6);
		REQUIRE(firstSample.phaseIndex == 1);
		REQUIRE(firstSample.weight == 13);

		REQUIRE(secondSample.inputIndex == 6);
		REQUIRE(secondSample.phaseIndex == 2);
		REQUIRE(secondSample.weight == 1);
	}
	SECTION("Rollover") {
		const auto [firstSample, secondSample] = impl::InputIndex2Sample({ 27, 14 }, 5);
		REQUIRE(firstSample.inputIndex == 1);
		REQUIRE(firstSample.phaseIndex == 4);
		REQUIRE(firstSample.weight == 5);

		REQUIRE(secondSample.inputIndex == 2);
		REQUIRE(secondSample.phaseIndex == 0);
		REQUIRE(secondSample.weight == 9);
	}
}


TEST_CASE("Resampling: Dot product sample", "[Interpolation]") {
	const Signal<int> signal = { 1, 2, 3, 6, 5, 7 };
	const Signal<int> filter = { -1, 3, -2 }; // Convolution: -2, 3, -1
	REQUIRE(-2 == impl::DotProductSample(signal, filter, 0));
	REQUIRE(-1 == impl::DotProductSample(signal, filter, 2));
	REQUIRE(-5 == impl::DotProductSample(signal, filter, 5));
	REQUIRE(-7 == impl::DotProductSample(signal, filter, 7));
}


//------------------------------------------------------------------------------
// Resampling functions
//------------------------------------------------------------------------------

TEST_CASE("Decimate", "[Interpolation]") {
	SECTION("Out-of-place") {
		const Signal<float> input = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
		const Signal<float> result = Decimate(input, 3);
		const Signal<float> expected = { 1, 4, 7, 10 };
		REQUIRE(result == expected);
	}
	SECTION("In-place") {
		Signal<float> data = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
		Decimate(data, data, 3);
		const Signal<float> expected = { 1, 4, 7, 10 };
		REQUIRE(SignalView(data).subsignal(0, 4) == expected);
	}
}


TEST_CASE("Expand", "[Interpolation]") {
	SECTION("Out-of-place") {
		const Signal<float> input = { 1, 2, 3 };
		const Signal<float> result = Expand(input, 3);
		const Signal<float> expected = { 1, 0, 0, 2, 0, 0, 3, 0, 0 };
		REQUIRE(result == expected);
	}
	SECTION("In-place") {
		Signal<float> data = { 1, 2, 3, 0, 0, 0, 0, 0, 0 };
		Expand(data, data, 3);
		const Signal<float> expected = { 1, 0, 0, 2, 0, 0, 3, 0, 0 };
		REQUIRE(data == expected);
	}
}


TEST_CASE("Interpolation full & central", "[Interpolation]") {
	constexpr int factor = 5;
	constexpr int inputSize = 1024;

	for (const int filterSize : { 31, 33, 2047 }) {
		const auto center = inputSize / 2;
		const auto input = SmoothTestSignal<float>(inputSize, center, 0.3 * float(inputSize));
		const auto filter = DesignFilter<float, TIME_DOMAIN>(filterSize, Fir.Lowpass.Windowed.Cutoff(1.0f / factor));
		const auto polyphase = PolyphaseNormalized(PolyphaseReorder(filter, factor));

		SECTION("Full") {
			const auto expectedLength = InterpolLength(inputSize, filterSize, factor, CONV_FULL);
			const auto expectedCenter = InterpolDelay(filterSize) + center * factor;
			const auto expected = SmoothTestSignal<float>(expectedLength, expectedCenter, 0.3 * float(inputSize * factor));
			const auto result = Interpolate(input, polyphase, CONV_FULL).first;
			REQUIRE(expected.size() == result.size());
			REQUIRE(Max(Abs(expected - result)) < 1e-4f);
		}
		SECTION("Central") {
			const auto expectedLength = InterpolLength(inputSize, filterSize, factor, CONV_CENTRAL);
			const auto expectedCenter = InterpolDelay(filterSize) + center * factor - filterSize + 1;
			const auto expected = SmoothTestSignal<float>(expectedLength, expectedCenter, 0.3 * float(inputSize * factor));
			const auto result = Interpolate(input, polyphase, CONV_CENTRAL).first;
			REQUIRE(expected.size() == result.size());
			REQUIRE(Max(Abs(expected - result)) < 1e-4f);
		}
	}
}


TEST_CASE("Resampling full & central", "[Interpolation]") {
	constexpr int superSampling = 16;
	constexpr int inputSize = 1536;
	constexpr Rational<int64_t> sampleRates(7, 11);


	for (const int filterSize : { 127, 2047 }) {
		const auto center = inputSize / 2;
		const auto input = SmoothTestSignal<float>(inputSize, center, 0.3 * float(inputSize));
		const auto cutoff = ResampleFilterCutoff(sampleRates, superSampling);
		const auto filter = DesignFilter<float, TIME_DOMAIN>(filterSize, Fir.Lowpass.Windowed.Cutoff(cutoff));
		const auto polyphase = PolyphaseNormalized(PolyphaseReorder(filter, superSampling));

		SECTION("Full") {
			const auto expectedLength = ResampleLength(inputSize, filterSize, superSampling, sampleRates, CONV_FULL);
			const auto expectedCenter = ResampleDelay(filterSize, superSampling, sampleRates) + int64_t(center) / sampleRates;
			const auto expected = SmoothTestSignal<float>(floor(expectedLength), float(expectedCenter), 0.3 * float(inputSize) / float(sampleRates));
			const auto result = Resample(input, polyphase, sampleRates, CONV_FULL).first;
			REQUIRE(expected.size() == result.size());
			REQUIRE(Max(Abs(expected - result)) < 1e-4f);
		}
		SECTION("Central") {
			const auto expectedLength = ResampleLength(inputSize, filterSize, superSampling, sampleRates, CONV_CENTRAL);
			const auto expectedCenter = ResampleDelay(filterSize, superSampling, sampleRates) + int64_t(center) / sampleRates - Rational<int64_t>(filterSize - 1, superSampling) / sampleRates;
			const auto expected = SmoothTestSignal<float>(floor(expectedLength), float(expectedCenter), 0.3 * float(inputSize) / float(sampleRates));
			const auto result = Resample(input, polyphase, sampleRates, CONV_CENTRAL).first;
			REQUIRE(expected.size() == result.size());
			REQUIRE(Max(Abs(expected - result)) < 1e-4f);
		}
	}
}


TEST_CASE("Resampling various sample rates", "[Interpolation]") {
	constexpr int inputSize = 3072;

	constexpr std::array scenarios = {
		std::tuple{ Rational<int64_t>(7, 11), 16, std::array{ 255, 6143 } },
		std::tuple{ Rational<int64_t>(7, 191), 32, std::array{ 511, 8143 } },
		std::tuple{ Rational<int64_t>(191, 17), 8, std::array{ 2047, 16383 } },
		std::tuple{ Rational<int64_t>(11, 7), 16, std::array{ 255, 4095 } },
	};


	for (const auto [sampleRates, superSampling, filterSizes] : scenarios) {
		const auto center = inputSize / 2;
		const auto input = SmoothTestSignal<float>(inputSize, center, 0.5 * float(inputSize));

		for (const int filterSize : filterSizes) {
			const auto cutoff = ResampleFilterCutoff(sampleRates, superSampling);
			const auto filter = DesignFilter<float, TIME_DOMAIN>(filterSize, Fir.Lowpass.Windowed.Cutoff(0.85f * cutoff).Window(windows::flattop));
			const auto polyphase = PolyphaseNormalized(PolyphaseReorder(filter, superSampling));

			INFO("Sample rates: " << sampleRates.Numerator() << "->" << sampleRates.Denominator()
								  << ", SS: " << superSampling << "x"
								  << ", filter: " << filterSize);
			const auto expectedLength = ResampleLength(inputSize, filterSize, superSampling, sampleRates, CONV_FULL);
			const auto expectedCenter = ResampleDelay(filterSize, superSampling, sampleRates) + int64_t(center) / sampleRates;
			const auto expected = SmoothTestSignal<float>(floor(expectedLength), float(expectedCenter), 0.5 * float(inputSize) / float(sampleRates));
			const auto result = Resample(input, polyphase, sampleRates, CONV_FULL).first;
			const auto error = expected - result;
			REQUIRE(expected.size() == result.size());
			REQUIRE(Max(Abs(expected - result)) < 3e-4f);
		}
	}
}


TEST_CASE("Interpolation: sequential chunks", "[Interpolation]") {
	constexpr size_t inputSize = 3072;
	constexpr size_t filterSize = 511;
	constexpr size_t factor = 16;

	const auto center = inputSize / 2;
	const auto input = SmoothTestSignal<float>(inputSize, center, 0.5 * float(inputSize));
	const auto cutoff = InterpolFilterCutoff(factor);
	const auto filter = DesignFilter<float, TIME_DOMAIN>(filterSize, Fir.Lowpass.Windowed.Cutoff(0.85f * cutoff).Window(windows::flattop));
	const auto polyphase = PolyphaseNormalized(PolyphaseReorder(filter, factor));

	const auto expectedLength = InterpolLength(inputSize, filterSize, factor, CONV_FULL);
	const auto expectedCenter = InterpolDelay(filterSize) + center * factor;
	const auto expected = SmoothTestSignal<float>(expectedLength, expectedCenter, 0.5 * float(inputSize * factor));

	Signal<float> result(floor(expectedLength));
	InterpolSuspensionPoint sp = { .firstInputSample = 0, .outputOffset = 0 };
	constexpr size_t chunkSize = 256;
	auto outputIt = result.begin();
	auto inputIt = input.begin();
	while (outputIt != result.end()) {
		const auto truncatedChunkSize = std::min(chunkSize, size_t(result.end() - outputIt));
		const auto outputChunk = SignalView<float>(outputIt, truncatedChunkSize);
		const auto inputChunk = SignalView<const float>(inputIt, input.end());
		sp = Interpolate(outputChunk, inputChunk, polyphase, sp.outputOffset);
		outputIt += truncatedChunkSize;
		inputIt += sp.firstInputSample;
	}

	const auto error = expected - result;
	REQUIRE(expected.size() == result.size());
	REQUIRE(Max(Abs(expected - result)) < 3e-4f);
}


TEST_CASE("Resampling: sequential chunks", "[Interpolation]") {
	constexpr size_t inputSize = 3072;
	constexpr size_t filterSize = 511;
	constexpr size_t superSampling = 16;
	constexpr Rational<int64_t> sampleRates = { 7, 11 };

	const auto center = inputSize / 2;
	const auto input = SmoothTestSignal<float>(inputSize, center, 0.5 * float(inputSize));
	const auto cutoff = ResampleFilterCutoff(sampleRates, superSampling);
	const auto filter = DesignFilter<float, TIME_DOMAIN>(filterSize, Fir.Lowpass.Windowed.Cutoff(0.85f * cutoff).Window(windows::flattop));
	const auto polyphase = PolyphaseNormalized(PolyphaseReorder(filter, superSampling));

	const auto expectedLength = ResampleLength(inputSize, filterSize, superSampling, sampleRates, CONV_FULL);
	const auto expectedCenter = ResampleDelay(filterSize, superSampling, sampleRates) + int64_t(center) / sampleRates;
	const auto expected = SmoothTestSignal<float>(floor(expectedLength), float(expectedCenter), 0.5 * float(inputSize) / float(sampleRates));

	Signal<float> result(floor(expectedLength));
	ResampleSuspensionPoint sp = { .firstInputSample = 0, .outputOffset{ int64_t(0) } };
	constexpr size_t chunkSize = 256;
	auto outputIt = result.begin();
	auto inputIt = input.begin();
	while (outputIt != result.end()) {
		const auto truncatedChunkSize = std::min(chunkSize, size_t(result.end() - outputIt));
		const auto outputChunk = SignalView<float>(outputIt, truncatedChunkSize);
		const auto inputChunk = SignalView<const float>(inputIt, input.end());
		sp = Resample(outputChunk, inputChunk, polyphase, sampleRates, sp.outputOffset);
		outputIt += truncatedChunkSize;
		inputIt += sp.firstInputSample;
	}

	const auto error = expected - result;
	REQUIRE(expected.size() == result.size());
	REQUIRE(Max(Abs(expected - result)) < 3e-4f);
}