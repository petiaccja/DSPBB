#include "PlaySound.hpp"

#include "../RtAudio/RtAudio.h"

#include <stdexcept>
#include <thread>

using namespace dspbb;

static void Play(uint64_t sampleRate, unsigned numChannels, RtAudioCallback callbackWrapper, void* callback);

int StereoCallbackWrapper(void* outputBuffer, void*, unsigned int nBufferFrames, double, RtAudioStreamStatus, void* userData) {
	auto* const callback = static_cast<const PlayStereoCallback*>(userData);
	Signal<float> left(nBufferFrames, 0.0f);
	Signal<float> right(nBufferFrames, 0.0f);
	const size_t writtenFrames = (*callback)(left, right);

	float* outputIt = static_cast<float*>(outputBuffer);
	auto leftIt = left.begin();
	auto rightIt = right.begin();
	for (; leftIt != left.end() && rightIt != right.end(); ++leftIt, ++rightIt, outputIt += 2) {
		outputIt[0] = *leftIt;
		outputIt[1] = *rightIt;
	}

	return writtenFrames == nBufferFrames ? 0 : 1;
}

int MonoCallbackWrapper(void* outputBuffer, void*, unsigned int nBufferFrames, double, RtAudioStreamStatus, void* userData) {
	auto* const callback = static_cast<PlayMonoCallback*>(userData);
	Signal<float> out(nBufferFrames, 0.0f);
	const size_t writtenFrames = (*callback)(out);
	auto* outputIt = static_cast<float*>(outputBuffer);
	for (auto& v : out) {
		*outputIt++ = v;
		*outputIt++ = v;
	}
	return writtenFrames == nBufferFrames ? 0 : 1;
}

void PlayStereo(uint64_t sampleRate, const PlayStereoCallback& callback) {
	Play(sampleRate, 2, &StereoCallbackWrapper, (void*)&callback);
}

void PlayMono(uint64_t sampleRate, const PlayMonoCallback& callback) {
	Play(sampleRate, 2, &MonoCallbackWrapper, (void*)&callback);
}

void PlayStereo(uint64_t sampleRate, SignalView<float> samplesLeft, SignalView<float> samplesRight) {
	size_t currentSample = 0;
	PlayStereo(sampleRate, [&samplesLeft, &samplesRight, &currentSample](SignalView<float> leftOut, SignalView<float> rightOut) -> size_t {
		assert(leftOut.size() == rightOut.size());
		const size_t validSize = std::min(samplesLeft.size() - currentSample, leftOut.size());
		const auto leftChunk = samplesLeft.subsignal(currentSample, validSize);
		const auto rightChunk = samplesRight.subsignal(currentSample, validSize);
		std::copy(leftChunk.begin(), leftChunk.end(), leftOut.begin());
		std::copy(rightChunk.begin(), rightChunk.end(), rightOut.begin());
		currentSample += validSize;
		return validSize;
	});
}

void PlayMono(uint64_t sampleRate, SignalView<float> samples) {
	size_t currentSample = 0;
	PlayMono(sampleRate, [&samples, &currentSample](SignalView<float> out) -> size_t {
		const size_t validSize = std::min(samples.size() - currentSample, out.size());
		const auto chunk = samples.subsignal(currentSample, validSize);
		std::copy(chunk.begin(), chunk.end(), out.begin());
		currentSample += validSize;
		return validSize;
	});
}

static void Play(uint64_t sampleRate, unsigned numChannels, RtAudioCallback callbackWrapper, void* callback) {
	RtAudio soundOutput;
	if (soundOutput.getDeviceCount() < 1) {
		throw std::runtime_error("Could not find sound output devices.");
	}
	RtAudio::StreamParameters parameters;
	parameters.deviceId = soundOutput.getDefaultOutputDevice();
	parameters.nChannels = numChannels;
	parameters.firstChannel = 0;
	unsigned numOutputFrames = static_cast<unsigned>((sampleRate + 5) / 6);

	soundOutput.openStream(&parameters, nullptr, RTAUDIO_FLOAT32, (unsigned)sampleRate, &numOutputFrames, callbackWrapper, callback);
	soundOutput.startStream();
	while (soundOutput.isStreamRunning()) {
		std::this_thread::sleep_for(std::chrono::milliseconds(16));
	}
}