#include "dspbb/Generators/Sine.hpp"

#include <dspbb/DSP/PolyphaseFilter.hpp>

#include <Catch2/catch.hpp>
#include <cmath>

using namespace dspbb;


TEST_CASE("Polyphase upsample", "[AudioFramework:Polyphase]") {
	PolyphaseFilter<float> polyphase{ 44100, 22050.f, 4, 31 };

	const auto signal = SineWave<float, TIME_DOMAIN>(200, 44100, 100.0);

	auto filtered = polyphase(SignalView<const float, TIME_DOMAIN>{ signal.begin(), signal.end() }, convolution::full);
	REQUIRE(filtered.Size() > signal.Size());
}
