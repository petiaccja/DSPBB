#include "dspbb/Generators/Sine.hpp"

#include <dspbb/Filtering/PolyphaseFilter.hpp>
#include <dspbb/Math/Statistics.hpp>

#include <catch2/catch.hpp>
#include <cmath>

using namespace dspbb;


TEST_CASE("Polyphase upsample", "[Polyphase]") {
	constexpr int factor = 4;
	const auto filter = FirLowPassWindowed(22050.f, factor * 44100, HammingWindow<float>(128));
	const PolyphaseFilter<float> polyphase{ AsConstView(filter), factor };

	const auto signal = SineWave<float, TIME_DOMAIN>(221, 44100, 100.0);
	TimeSignal<float> interspersed(signal.Size()*factor, 0.0f);
	for (size_t i = 0; i < signal.Size(); ++i) {
		interspersed[factor*i] = signal[i];
	}

	const auto filtered = polyphase(SignalView<const float, TIME_DOMAIN>{ signal.begin(), signal.end() }, convolution::full);
	const auto control = Convolution(interspersed, factor*filter, convolution::full);
	REQUIRE(filtered.Size() == control.Size() - factor + 1);
	const auto diff = Max(Abs(filtered - TimeSignalView<const float>{ control.begin(), filtered.Size() }));
	REQUIRE(diff < 0.001f);
}
