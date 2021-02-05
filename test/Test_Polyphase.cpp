#include <dspbb/DSP/PolyphaseFilter.hpp>

#include <Catch2/catch.hpp>
#include <cmath>

using namespace dspbb;


TEST_CASE("Polyphase upsample", "[AudioFramework:Polyphase]") {
	PolyphaseFilter<float> polyphase{ 44100, 22050.f, 4, 31 };

	TimeSignalF signal(200);
	for (size_t i = 0; i < signal.Size(); ++i) {
		signal[i] = std::sin(2.0f * pi_v<float> * float(i) / 44100.f * 100.f);
	}

	auto filtered = polyphase(SignalView<const float, TIME_DOMAIN>{ signal.begin(), signal.end() }, convolution::full);
	REQUIRE(filtered.Size() > signal.Size());
}
