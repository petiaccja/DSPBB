#pragma once

#include <dspbb/Primitives/Signal.hpp>

#include <filesystem>

struct StereoSound {
	dspbb::Signal<float> leftChannel;
	dspbb::Signal<float> rightChannel;
	uint64_t sampleRate;
};

StereoSound LoadStereoSound(const std::filesystem::path& file);