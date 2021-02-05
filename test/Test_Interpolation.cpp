#include <dspbb/DSP/Interpolation.hpp>

#include <Catch2/catch.hpp>
#include <cmath>

using namespace dspbb;


TEST_CASE("Polyphase resample", "[AudioFramework:Interpolation]") {
	uint64_t sampleRateIn = 16000;
	uint64_t sampleRateOut = 44100;
	PolyphaseFilter<float> polyphase{ sampleRateIn, float(std::min(sampleRateIn, sampleRateOut)) / 2.0f, 7, 31 };

	TimeSignalF signal(260);
	for (size_t i = 0; i < signal.Size(); ++i) {
		signal[i] = std::sin(2.0f * pi_v<float> * float(i) / float(sampleRateIn) * 100.f);
	}

	auto filtered = InterpolatePolyphase(SignalView<const float, TIME_DOMAIN>{ signal }, polyphase, sampleRateIn, sampleRateOut, 0, 100);
	REQUIRE(filtered.Size() > signal.Size());
}
