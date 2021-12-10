#pragma once

#include <functional>
#include <dspbb/Primitives/SignalView.hpp>

using PlayStereoCallback = std::function<size_t(dspbb::SignalView<float>, dspbb::SignalView<float>)>;
using PlayMonoCallback = std::function<size_t(dspbb::SignalView<float>)>;

void PlayStereo(uint64_t sampleRate, const PlayStereoCallback& callback);
void PlayMono(uint64_t sampleRate, const PlayMonoCallback& callback);

void PlayStereo(uint64_t sampleRate, dspbb::SignalView<float> samplesLeft, dspbb::SignalView<float> samplesRight);
void PlayMono(uint64_t sampleRate, dspbb::SignalView<float> samples);