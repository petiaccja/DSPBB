#include "../TestUtils.hpp"

#include <dspbb/Filtering/FIR.hpp>
#include <dspbb/Filtering/MeasureFilter.hpp>
#include <dspbb/Filtering/Resample.hpp>
#include <dspbb/Generators/Waveforms.hpp>
#include <dspbb/Math/Convolution.hpp>
#include <dspbb/Math/Statistics.hpp>

#include <catch2/catch.hpp>

using namespace dspbb;


//------------------------------------------------------------------------------
// Filter application helpers
//------------------------------------------------------------------------------

TEST_CASE("Filter state continuity", "[FIR]") {
	constexpr int taps = 7;
	constexpr int length = 80;

	const auto signal = RandomSignal<double, TIME_DOMAIN>(length);
	const auto filter = DesignFilter<double, TIME_DOMAIN>(taps, Fir.Lowpass.LeastSquares.Cutoff(0.3f, 0.33f));

	const auto expected = Convolution(signal, filter, 0, length);

	Signal<double> state(taps - 1, 0.0f);
	Signal<double> result(length);

	SECTION("Convolution large") {
		constexpr int step = 40;
		static_assert(length % step == 0);
		for (size_t i = 0; i < length; i += step) {
			Filter(AsView(result).subsignal(i, step), AsView(signal).subsignal(i, step), filter, state, FILTER_CONV);
		}
	}
	SECTION("OLA large") {
		constexpr int step = 40;
		static_assert(length % step == 0);
		for (size_t i = 0; i < length; i += step) {
			Filter(AsView(result).subsignal(i, step), AsView(signal).subsignal(i, step), filter, state, FILTER_OLA);
		}
	}
	SECTION("Convolution small") {
		constexpr int step = 4;
		static_assert(length % step == 0);
		for (size_t i = 0; i < length; i += step) {
			Filter(AsView(result).subsignal(i, step), AsView(signal).subsignal(i, step), filter, state, FILTER_CONV);
		}
	}
	SECTION("OLA small") {
		constexpr int step = 4;
		static_assert(length % step == 0);
		for (size_t i = 0; i < length; i += step) {
			Filter(AsView(result).subsignal(i, step), AsView(signal).subsignal(i, step), filter, state, FILTER_OLA);
		}
	}
	SECTION("Convolution copy") {
		constexpr int step = 4;
		static_assert(length % step == 0);
		for (size_t i = 0; i < length; i += step) {
			const auto batch = Filter(AsView(signal).subsignal(i, step), filter, state, FILTER_CONV);
			const auto outBatch = AsView(result).subsignal(i, step);
			std::copy(batch.begin(), batch.end(), outBatch.begin());
		}
	}
	SECTION("OLA copy") {
		constexpr int step = 4;
		static_assert(length % step == 0);
		for (size_t i = 0; i < length; i += step) {
			const auto batch = Filter(AsView(signal).subsignal(i, step), filter, state, FILTER_OLA);
			const auto outBatch = AsView(result).subsignal(i, step);
			std::copy(batch.begin(), batch.end(), outBatch.begin());
		}
	}

	REQUIRE(Max(Abs(result - expected)) < 1e-7);
}

TEST_CASE("Filter central", "[FIR]") {
	constexpr int taps = 7;
	constexpr int length = 80;

	const auto signal = RandomSignal<double, TIME_DOMAIN>(length);
	const auto filter = DesignFilter<double, TIME_DOMAIN>(taps, Fir.Lowpass.LeastSquares.Cutoff(0.3f, 0.33f));

	const auto expected = Convolution(signal, filter, CONV_CENTRAL);

	SECTION("Convolution") {
		const auto result = Filter(signal, filter, CONV_CENTRAL, FILTER_CONV);
		REQUIRE(Max(Abs(result - expected)) < 1e-7);
	}
	SECTION("OLA") {
		const auto result = Filter(signal, filter, CONV_CENTRAL, FILTER_OLA);
		REQUIRE(Max(Abs(result - expected)) < 1e-7);
	}
}

TEST_CASE("Filter full", "[FIR]") {
	constexpr int taps = 7;
	constexpr int length = 80;

	const auto signal = RandomSignal<double, TIME_DOMAIN>(length);
	const auto filter = DesignFilter<double, TIME_DOMAIN>(taps, Fir.Lowpass.LeastSquares.Cutoff(0.3f, 0.33f));

	const auto expected = Convolution(signal, filter, CONV_FULL);

	SECTION("Convolution") {
		const auto result = Filter(signal, filter, CONV_FULL, FILTER_CONV);
		REQUIRE(Max(Abs(result - expected)) < 1e-7);
	}

	SECTION("OLA") {
		const auto result = Filter(signal, filter, CONV_FULL, FILTER_OLA);
		REQUIRE(Max(Abs(result - expected)) < 1e-7);
	}
}

//------------------------------------------------------------------------------
// Helpers
//------------------------------------------------------------------------------

template <class SignalT>
bool IsSymmetric(const SignalT& signal) {
	auto beg = signal.begin();
	auto end = signal.rbegin();
	while (beg <= end.base()) {
		if (*beg != Approx(*end).margin(1e-7f)) {
			return false;
		}
		++beg;
		++end;
	}
	return true;
}

template <class SignalT>
bool IsAntiSymmetric(const SignalT& signal) {
	auto beg = signal.begin();
	auto end = signal.rbegin();
	while (beg <= end.base()) {
		if (*beg != Approx(-*end).margin(1e-7f)) {
			return false;
		}
		++beg;
		++end;
	}
	return true;
}

constexpr auto TestArbitraryResponse = [](float x) {
	return 2.0f * x - 1.5f * x * x - 0.5f * std::pow(x - 1.f, 3.f);
};


//------------------------------------------------------------------------------
// Window method
//------------------------------------------------------------------------------

TEST_CASE("Windowed low-pass", "[FIR]") {
	constexpr size_t numTaps = 255;
	static constexpr float cutoff = 0.3f;

	const auto impulse = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Lowpass.Windowed.Cutoff(cutoff).Window(windows::blackman));
	REQUIRE(impulse.size() == numTaps);
	REQUIRE(IsSymmetric(impulse));

	const auto [amplitude, phase] = FrequencyResponse(impulse);
	const auto params = MeasureLowpassFilter(amplitude);
	REQUIRE(params.passbandEdge < cutoff);
	REQUIRE(params.stopbandEdge > cutoff);
	REQUIRE(params.passbandRipple < 0.05f);
	REQUIRE(params.stopbandAtten < 0.05f);
}

TEST_CASE("Windowed high-pass", "[FIR]") {
	constexpr size_t numTaps = 255;
	static constexpr float cutoff = 0.3f;

	const auto impulse = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Highpass.Windowed.Cutoff(cutoff).Window(windows::blackman));
	REQUIRE(impulse.size() == numTaps);
	REQUIRE(IsSymmetric(impulse));

	const auto [amplitude, phase] = FrequencyResponse(impulse);
	const auto params = MeasureHighpassFilter(amplitude);
	REQUIRE(params.stopbandEdge < cutoff);
	REQUIRE(params.passbandEdge > cutoff);
	REQUIRE(params.stopbandAtten < 0.05f);
	REQUIRE(params.passbandRipple < 0.05f);
}

TEST_CASE("Windowed band-pass", "[FIR]") {
	constexpr size_t numTaps = 255;
	static constexpr float bandLow = 0.3f;
	static constexpr float bandHigh = 0.6f;

	const auto impulse = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Bandpass.Windowed.Band(bandLow, bandHigh).Window(windows::blackman));
	REQUIRE(impulse.size() == numTaps);
	REQUIRE(IsSymmetric(impulse));

	const auto [amplitude, phase] = FrequencyResponse(impulse);
	const auto params = MeasureBandpassFilter(amplitude);
	REQUIRE(params.lowerStopbandEdge < bandLow);
	REQUIRE(params.passbandLowerEdge > bandLow);
	REQUIRE(params.passbandUpperEdge < bandHigh);
	REQUIRE(params.upperStopbandEdge > bandHigh);
	REQUIRE(params.lowerStopbandAtten < 0.05f);
	REQUIRE(params.passbandRipple < 0.05f);
	REQUIRE(params.upperStopbandAtten < 0.05f);
}

TEST_CASE("Windowed band-stop", "[FIR]") {
	constexpr size_t numTaps = 255;
	static constexpr float bandLow = 0.3f;
	static constexpr float bandHigh = 0.6f;

	const auto impulse = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Bandstop.Windowed.Band(bandLow, bandHigh).Window(windows::blackman));
	REQUIRE(impulse.size() == numTaps);
	REQUIRE(IsSymmetric(impulse));

	const auto [amplitude, phase] = FrequencyResponse(impulse);
	const auto params = MeasureBandstopFilter(amplitude);
	REQUIRE(params.lowerPassbandEdge < bandLow);
	REQUIRE(params.stopbandLowerEdge > bandLow);
	REQUIRE(params.stopbandUpperEdge < bandHigh);
	REQUIRE(params.upperPassbandEdge > bandHigh);
	REQUIRE(params.lowerPassbandRipple < 0.05f);
	REQUIRE(params.stopbandAtten < 0.05f);
	REQUIRE(params.upperPassbandRipple < 0.05f);
}

TEST_CASE("Windowed arbitrary", "[FIR]") {
	constexpr size_t numTaps = 255;

	const auto impulse = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Arbitrary.Windowed.Response(TestArbitraryResponse).Window(windows::blackman));
	REQUIRE(impulse.size() == numTaps);
	REQUIRE(IsSymmetric(impulse));

	const auto [amplitude, phase] = FrequencyResponse(impulse);
	auto expected = LinSpace<float, FREQUENCY_DOMAIN>(0.0f, 1.0f, amplitude.size(), true);
	std::for_each(expected.begin(), expected.end(), [](auto& v) { v = TestArbitraryResponse(v); });
	REQUIRE(Max(Abs(amplitude - expected)) < 0.02f);
}

TEST_CASE("Windowed hilbert magnitude", "[FIR]") {
	const auto odd = DesignFilter<float, TIME_DOMAIN>(377, Fir.Hilbert.Windowed.Window(windows::blackman));
	const auto even = DesignFilter<float, TIME_DOMAIN>(376, Fir.Hilbert.Windowed.Window(windows::blackman));

	const auto [amplitudeOdd, phaseOdd] = FrequencyResponse(odd);
	const auto paramsOdd = MeasureBandpassFilter(amplitudeOdd);
	REQUIRE(paramsOdd.passbandLowerEdge < 0.05f);
	REQUIRE(paramsOdd.passbandUpperEdge > 0.95f);
	REQUIRE(paramsOdd.passbandRipple < 0.05f);

	const auto [amplitudeEven, phaseEven] = FrequencyResponse(even);
	const auto paramsEven = MeasureHighpassFilter(amplitudeEven);
	REQUIRE(paramsEven.passbandEdge < 0.05f);
	REQUIRE(paramsEven.passbandRipple < 0.05f);
}

TEST_CASE("Windowed methods equal", "[FIR]") {
	constexpr size_t numTaps = 127;
	constexpr float cutoff = 0.3f;
	constexpr float bandLow = 0.2f;
	constexpr float bandHigh = 0.6f;

	const auto lp1 = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Lowpass.Windowed.Cutoff(cutoff).Window(windows::blackman));
	const auto lp2 = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Lowpass.Windowed.Cutoff(cutoff).Window(windows::blackman.operator()<float, TIME_DOMAIN>(numTaps)));

	const auto hp1 = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Highpass.Windowed.Cutoff(cutoff).Window(windows::blackman));
	const auto hp2 = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Highpass.Windowed.Cutoff(cutoff).Window(windows::blackman.operator()<float, TIME_DOMAIN>(numTaps)));

	const auto bp1 = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Bandpass.Windowed.Band(bandLow, bandHigh).Window(windows::blackman));
	const auto bp2 = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Bandpass.Windowed.Band(bandLow, bandHigh).Window(windows::blackman.operator()<float, TIME_DOMAIN>(numTaps)));

	const auto bs1 = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Bandstop.Windowed.Band(bandLow, bandHigh).Window(windows::blackman));
	const auto bs2 = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Bandstop.Windowed.Band(bandLow, bandHigh).Window(windows::blackman.operator()<float, TIME_DOMAIN>(numTaps)));

	REQUIRE(Max(Abs(lp1 - lp2)) < 1e-4f);
	REQUIRE(Max(Abs(hp1 - hp2)) < 1e-4f);
	REQUIRE(Max(Abs(bp1 - bp2)) < 1e-4f);
	REQUIRE(Max(Abs(bs1 - bs2)) < 1e-4f);
}


//------------------------------------------------------------------------------
// Least squares method
//------------------------------------------------------------------------------

TEST_CASE("Least squares low-pass", "[FIR]") {
	constexpr size_t numTaps = 255;
	constexpr float cutoffBegin = 0.28f;
	constexpr float cutoffEnd = 0.32f;
	constexpr float width = cutoffEnd - cutoffBegin;

	const auto impulse = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Lowpass.LeastSquares.Cutoff(cutoffBegin, cutoffEnd));
	REQUIRE(impulse.size() == numTaps);
	REQUIRE(IsSymmetric(impulse));

	const auto [amplitude, phase] = FrequencyResponse(impulse);
	const auto params = MeasureLowpassFilter(amplitude);
	REQUIRE(params.passbandEdge > cutoffBegin - width / 2);
	REQUIRE(params.passbandEdge < cutoffBegin + width / 2);
	REQUIRE(params.stopbandEdge > cutoffEnd - width / 2);
	REQUIRE(params.stopbandEdge < cutoffEnd + width / 2);
	REQUIRE(params.passbandRipple < 0.05f);
	REQUIRE(params.stopbandAtten < 0.05f);
}

TEST_CASE("Least squares high-pass", "[FIR]") {
	constexpr size_t numTaps = 255;
	constexpr float cutoffBegin = 0.28f;
	constexpr float cutoffEnd = 0.32f;
	constexpr float width = cutoffEnd - cutoffBegin;

	const auto impulse = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Highpass.LeastSquares.Cutoff(cutoffBegin, cutoffEnd));
	REQUIRE(impulse.size() == numTaps);
	REQUIRE(IsSymmetric(impulse));

	const auto [amplitude, phase] = FrequencyResponse(impulse);
	const auto params = MeasureHighpassFilter(amplitude);
	REQUIRE(params.stopbandEdge > cutoffBegin - width / 2);
	REQUIRE(params.stopbandEdge < cutoffBegin + width / 2);
	REQUIRE(params.passbandEdge > cutoffEnd - width / 2);
	REQUIRE(params.passbandEdge < cutoffEnd + width / 2);
	REQUIRE(params.stopbandAtten < 0.05f);
	REQUIRE(params.passbandRipple < 0.05f);
}


TEST_CASE("Least squares band-pass", "[FIR]") {
	constexpr size_t numTaps = 255;
	constexpr float bandLowBegin = 0.28f;
	constexpr float bandLowEnd = 0.32f;
	constexpr float bandHighBegin = 0.58f;
	constexpr float bandHighEnd = 0.65f;
	constexpr float lowWidth = bandLowEnd - bandLowBegin;
	constexpr float highWidth = bandHighEnd - bandHighBegin;

	const auto impulse = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Bandpass.LeastSquares.Band(bandLowBegin, bandLowEnd, bandHighBegin, bandHighEnd).Weight(1.0f, 0.1f, 1.0f, 0.1f, 1.0f));
	REQUIRE(impulse.size() == numTaps);
	REQUIRE(IsSymmetric(impulse));

	const auto [amplitude, phase] = FrequencyResponse(impulse);
	const auto params = MeasureBandpassFilter(amplitude);
	REQUIRE(params.lowerStopbandEdge > bandLowBegin - lowWidth / 2);
	REQUIRE(params.lowerStopbandEdge < bandLowBegin + lowWidth / 2);
	REQUIRE(params.passbandLowerEdge > bandLowEnd - lowWidth / 2);
	REQUIRE(params.passbandLowerEdge < bandLowEnd + lowWidth / 2);
	REQUIRE(params.passbandUpperEdge > bandHighBegin - highWidth / 2);
	REQUIRE(params.passbandUpperEdge < bandHighBegin + highWidth / 2);
	REQUIRE(params.upperStopbandEdge > bandHighEnd - highWidth / 2);
	REQUIRE(params.upperStopbandEdge < bandHighEnd + highWidth / 2);
	REQUIRE(params.lowerStopbandAtten < 0.05f);
	REQUIRE(params.passbandRipple < 0.05f);
	REQUIRE(params.upperStopbandAtten < 0.05f);
}

TEST_CASE("Least squares band-stop", "[FIR]") {
	constexpr size_t numTaps = 255;
	static constexpr float bandLowBegin = 0.28f;
	static constexpr float bandLowEnd = 0.32f;
	static constexpr float bandHighBegin = 0.58f;
	static constexpr float bandHighEnd = 0.65f;
	constexpr float lowWidth = bandLowEnd - bandLowBegin;
	constexpr float highWidth = bandHighEnd - bandHighBegin;

	const auto impulse = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Bandstop.LeastSquares.Band(bandLowBegin, bandLowEnd, bandHighBegin, bandHighEnd).Weight(1.0f, 0.1f, 1.0f, 0.1f, 1.0f));
	REQUIRE(impulse.size() == numTaps);
	REQUIRE(IsSymmetric(impulse));

	const auto [amplitude, phase] = FrequencyResponse(impulse);
	const auto params = MeasureBandstopFilter(amplitude);
	REQUIRE(params.lowerPassbandEdge > bandLowBegin - lowWidth / 2);
	REQUIRE(params.lowerPassbandEdge < bandLowBegin + lowWidth / 2);
	REQUIRE(params.stopbandLowerEdge > bandLowEnd - lowWidth / 2);
	REQUIRE(params.stopbandLowerEdge < bandLowEnd + lowWidth / 2);
	REQUIRE(params.stopbandUpperEdge > bandHighBegin - highWidth / 2);
	REQUIRE(params.stopbandUpperEdge < bandHighBegin + highWidth / 2);
	REQUIRE(params.upperPassbandEdge > bandHighEnd - highWidth / 2);
	REQUIRE(params.upperPassbandEdge < bandHighEnd + highWidth / 2);
	REQUIRE(params.lowerPassbandRipple < 0.05f);
	REQUIRE(params.stopbandAtten < 0.05f);
	REQUIRE(params.upperPassbandRipple < 0.05f);
}

TEST_CASE("Least squares arbitrary", "[FIR]") {
	constexpr size_t numTaps = 255;

	const auto impulse = DesignFilter<float, TIME_DOMAIN>(numTaps, Fir.Arbitrary.LeastSquares.Response(TestArbitraryResponse));
	REQUIRE(impulse.size() == numTaps);
	REQUIRE(IsSymmetric(impulse));

	const auto [amplitude, phase] = FrequencyResponse(impulse);
	auto expected = LinSpace<float, FREQUENCY_DOMAIN>(0.0f, 1.0f, amplitude.size(), true);
	std::for_each(expected.begin(), expected.end(), [](auto& v) { v = TestArbitraryResponse(v); });
	REQUIRE(Max(Abs(amplitude - expected)) < 0.02f);
}

TEST_CASE("Least squares hilbert magnitude", "[FIR]") {
	const float transition = 0.02f;
	const auto responseDesc = Fir.Hilbert.LeastSquares.TransitionWidth(transition).TransitionWeight(0.05f);
	const auto odd = DesignFilter<float, TIME_DOMAIN>(155, responseDesc);
	const auto even = DesignFilter<float, TIME_DOMAIN>(154, responseDesc);

	const auto [amplitudeOdd, phaseOdd] = FrequencyResponse(odd);
	const auto paramsOdd = MeasureBandpassFilter(amplitudeOdd);
	REQUIRE(paramsOdd.passbandLowerEdge < 0.05f);
	REQUIRE(paramsOdd.passbandUpperEdge > 0.95f);
	REQUIRE(paramsOdd.passbandRipple < 0.05f);

	const auto [amplitudeEven, phaseEven] = FrequencyResponse(even);
	const auto paramsEven = MeasureHighpassFilter(amplitudeEven);
	REQUIRE(paramsEven.passbandEdge < 0.05f);
	REQUIRE(paramsEven.passbandRipple < 0.05f);
}

TEST_CASE("Least squares weights", "[FIR]") {
	const auto response = [](float f) {
		if (f < 0.5f) {
			return 1.0f;
		}
		return 0.0f;
	};
	const auto weightL = [](float f) {
		if (f < 0.45f) {
			return 3.0f;
		}
		if (f < 0.55f) {
			return 0.0f;
		}
		return 1.0f;
	};
	const auto weightH = [](float f) {
		if (f < 0.45f) {
			return 1.0f;
		}
		if (f < 0.55f) {
			return 0.0f;
		}
		return 3.0f;
	};

	const auto filterL = DesignFilter<float, TIME_DOMAIN>(27, Fir.Arbitrary.LeastSquares.Response(response).Weight(weightL));
	const auto filterH = DesignFilter<float, TIME_DOMAIN>(27, Fir.Arbitrary.LeastSquares.Response(response).Weight(weightH));

	const auto [amplitudeL, phaseL] = FrequencyResponse(filterL);
	const auto [amplitudeH, phaseH] = FrequencyResponse(filterH);

	const auto paramsL = MeasureLowpassFilter(amplitudeL);
	const auto paramsH = MeasureLowpassFilter(amplitudeH);

	REQUIRE(paramsL.passbandRipple < 0.5f * paramsH.passbandRipple);
	REQUIRE(0.5f * paramsL.stopbandAtten > paramsH.stopbandAtten);
}

//------------------------------------------------------------------------------
// Hilbert band transform special checks
//------------------------------------------------------------------------------

TEST_CASE("Hilbert odd form", "[FIR]") {
	const auto filter = DesignFilter<float, TIME_DOMAIN>(247, Fir.Hilbert.Windowed);
	REQUIRE(filter.size() == 247);
	REQUIRE(IsAntiSymmetric(filter));
	const auto nonZeroSamples = Decimate(filter, 2);
	const auto zeroSamples = Decimate(AsView(filter).subsignal(1), 2);
	REQUIRE(Max(zeroSamples) == 0.0f);
	REQUIRE(Min(Abs(nonZeroSamples)) > 0.0f);
	const auto firstHalf = AsView(nonZeroSamples).subsignal(0, nonZeroSamples.size() / 2);
	const auto secondHalf = AsView(nonZeroSamples).subsignal(nonZeroSamples.size() / 2);
	REQUIRE(Max(firstHalf) < 0.0f);
	REQUIRE(Min(secondHalf) > 0.0f);
}

TEST_CASE("Hilbert even form", "[FIR]") {
	const auto filter = DesignFilter<float, TIME_DOMAIN>(246, Fir.Hilbert.Windowed);
	REQUIRE(filter.size() == 246);
	REQUIRE(IsAntiSymmetric(filter));
	REQUIRE(Min(Abs(filter)) > 0.0f);
	const auto firstHalf = AsView(filter).subsignal(0, filter.size() / 2);
	const auto secondHalf = AsView(filter).subsignal(filter.size() / 2);
	REQUIRE(Max(firstHalf) < 0.0f);
	REQUIRE(Min(secondHalf) > 0.0f);
}

TEST_CASE("Hilbert odd small form", "[FIR]") {
	const auto filter = DesignFilter<float, TIME_DOMAIN>(19, Fir.Hilbert.Windowed);
	REQUIRE(filter.size() == 19);
	REQUIRE(IsAntiSymmetric(filter));
	const auto nonZeroSamples = Decimate(filter, 2);
	const auto zeroSamples = Decimate(AsView(filter).subsignal(1), 2);
	REQUIRE(Max(zeroSamples) == 0.0f);
	REQUIRE(Min(Abs(nonZeroSamples)) > 0.0f);
	const auto firstHalf = AsView(nonZeroSamples).subsignal(0, nonZeroSamples.size() / 2);
	const auto secondHalf = AsView(nonZeroSamples).subsignal(nonZeroSamples.size() / 2);
	REQUIRE(Max(firstHalf) < 0.0f);
	REQUIRE(Min(secondHalf) > 0.0f);
}

TEST_CASE("Hilbert even small form", "[FIR]") {
	const auto filter = DesignFilter<float, TIME_DOMAIN>(10, Fir.Hilbert.Windowed);
	REQUIRE(filter.size() == 10);
	REQUIRE(IsAntiSymmetric(filter));
	REQUIRE(Min(Abs(filter)) > 0.0f);
	const auto firstHalf = AsView(filter).subsignal(0, filter.size() / 2);
	const auto secondHalf = AsView(filter).subsignal(filter.size() / 2);
	REQUIRE(Max(firstHalf) < 0.0f);
	REQUIRE(Min(secondHalf) > 0.0f);
}


TEST_CASE("Hilbert odd phase shift", "[FIR]") {
	constexpr size_t testSignalSize = 4096;
	const auto filter = DesignFilter<float, TIME_DOMAIN>(377, Fir.Hilbert.Windowed);
	const auto testSignal = SineWave<float, TIME_DOMAIN>(testSignalSize, testSignalSize, 60.0) * GaussianWindow<float, TIME_DOMAIN>(testSignalSize, 0.25);
	const auto imaginarySignal = Convolution(filter, testSignal, CONV_CENTRAL);
	const auto realSignal = AsConstView(testSignal).subsignal(filter.size() / 2, imaginarySignal.size());
	REQUIRE(std::abs(DotProduct(realSignal, imaginarySignal) / testSignalSize) < 0.000001f);
	REQUIRE(Mean(realSignal) == Approx(Mean(imaginarySignal)).margin(0.001f));
}

TEST_CASE("Hilbert even phase shift", "[FIR]") {
	constexpr size_t testSignalSize = 4096;
	const auto filter = DesignFilter<float, TIME_DOMAIN>(376, Fir.Hilbert.Windowed);
	const auto testSignal = SineWave<float, TIME_DOMAIN>(testSignalSize, testSignalSize, 60.0) * GaussianWindow<float, TIME_DOMAIN>(testSignalSize, 0.25);
	const auto imaginarySignal = Convolution(filter, testSignal, CONV_CENTRAL);
	const auto realSignal = AsConstView(testSignal).subsignal(filter.size() / 2, imaginarySignal.size());
	REQUIRE(std::abs(DotProduct(realSignal, imaginarySignal) / testSignalSize) < 0.01f);
	REQUIRE(Mean(realSignal) == Approx(Mean(imaginarySignal)).margin(0.001f));
}