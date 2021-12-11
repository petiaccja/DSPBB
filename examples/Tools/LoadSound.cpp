#include "LoadSound.hpp"

#include <dspbb/Filtering/Interpolation.hpp>

#include <sndfile.hh>

using namespace dspbb;

StereoSound LoadStereoSound(const std::filesystem::path& file) {
	SndfileHandle soundFile(file.string().c_str());
	if (!soundFile) {
		throw std::runtime_error("Failed to open demo sound file.");
	}
	const size_t numFrames = soundFile.frames();
	const size_t numChannels = soundFile.channels();
	const uint64_t sampleRate = soundFile.samplerate();
	if (numChannels != 2) {
		throw std::runtime_error("You can only load stereo files at the moment.");
	}

	Signal<float> interleaved(numFrames * numChannels);
	soundFile.read(interleaved.Data(), interleaved.Size());

	Signal<float> leftChannel = Decimate(interleaved, 2);
	Signal<float> rightChannel = Decimate(AsView(interleaved).SubSignal(1), 2);

	return { std::move(leftChannel), std::move(rightChannel), sampleRate };
}
