#pragma once

#include <filesystem>
#include <dspbb/Primitives/Signal.hpp>

struct StereoSound {
	dspbb::Signal<float> leftChannel;
	dspbb::Signal<float> rightChannel;
	uint64_t sampleRate;
};

StereoSound LoadStereoSound(const std::filesystem::path& file);