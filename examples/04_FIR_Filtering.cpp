//------------------------------------------------------------------------------
// 04. FIR filtering
//
// This demo implements a simple bass/mid/treble equalizer using FIR filters.
// You can adjust the three parameters at you desire, and listen to how they
// make music sound.
//------------------------------------------------------------------------------

#include "Tools/LoadSound.hpp"
#include "Tools/PlaySound.hpp"

#include "dspbb/Filtering/FilterParameters.hpp"
#include <dspbb/Filtering/FIR.hpp>

#include <iostream>

using namespace dspbb;

// The transition bands determine the frequency regions where the bass and treble
// adjustments kick in. Bass can be considered to be below 200 Hz, while
// treble can be anything above 5000 Hz.
constexpr std::pair<float, float> transitionBass{ 160.0f, 320.f }; // Bass below 160 Hz.
constexpr std::pair<float, float> transitionTreble{ 5000.0f, 6500.f }; // Treble above 6000 Hz.

// A function that is zero for x<0 and 1 for x>1, with a smooth transition in between.
float Smoothstep(float x) {
	x = std::clamp(x, 0.0f, 1.0f);
	return 3 * x * x - 2 * x * x * x;
}

// Returns the amplification of the equalizing filter at a specific frequency, given the
// amplification levels for bass, mid and treble.
float EqualizedResponse(float frequency, float bass, float mid, float treble) {
	const float bassCurve = Smoothstep((transitionBass.second - frequency) / (transitionBass.second - transitionBass.first));
	const float midCurve = Smoothstep((transitionTreble.second - frequency) / (transitionTreble.second - transitionTreble.first))
						   + Smoothstep((frequency - transitionBass.first) / (transitionBass.second - transitionBass.first))
						   - 1.0f;
	const float trebleCurve = Smoothstep((frequency - transitionTreble.first) / (transitionTreble.second - transitionTreble.first));
	return bass * bassCurve + mid * midCurve + treble * trebleCurve;
}

class Equalizer {
public:
	Equalizer(size_t filterSize, uint64_t sampleRate)
		: m_filter(filterSize),
		  m_leftState(filterSize - 1),
		  m_rightState(filterSize - 1),
		  m_sampleRate(sampleRate) {}

	// Generates a new FIR filter to apply the equalization as given by the bass/mid/treble amplification levels.
	void SetLevels(float bass, float mid, float treble) {
		const auto normalizedResponse = [&](float nf) { return EqualizedResponse(nf * float(m_sampleRate) / 2.0f, bass, mid, treble); };
		// We will use a least squares FIR filter with no weighting and default grid size.
		FirFilter(m_filter, Arbitrary(LEAST_SQUARES).Response(normalizedResponse));
	}

	void Filter(SignalView<const float> leftIn, SignalView<const float> rightIn, SignalView<float> leftOut, SignalView<float> rightOut) {
		// The states here work the very same way as they do for the IIR filters (go check out the example).
		// We could have used plain convolution, but overlap-add seemed faster in debug mode.
		dspbb::Filter(leftOut, leftIn, m_filter, m_leftState, FILTER_OLA, 2048);
		dspbb::Filter(rightOut, rightIn, m_filter, m_rightState, FILTER_OLA, 2048);
	}

	void Reset() {
		// Filter state can be reset by filling it with zeros, just like IIR filters.
		std::fill(m_leftState.begin(), m_leftState.end(), 0.0f);
		std::fill(m_rightState.begin(), m_rightState.end(), 0.0f);
	}

	SignalView<const float> GetFilter() const { return AsConstView(m_filter); }

private:
	Signal<float> m_filter;
	Signal<float> m_leftState;
	Signal<float> m_rightState;
	const uint64_t m_sampleRate;
};

// A simple loop that let's you type the EQ parameters: bass, mid and treble.
// Keep in mind the values are not in decibels, but in ratios. You should not type
// values larger than 1, as those will just cause clipping in the output.
// Suppressing the mid and high frequencies only might also cause clipping, you may want to
// reduce bass just a bit too to avoid that.
int main() {
	try {
		// Feel free to replace this sample with whatever clip you like.
		const auto [left, right, sampleRate] = LoadStereoSound(DSPBB_EXAMPLES_DATA "sample.ogg");

		std::cout << "Welcome to the FIR filtering demo.\n"
				  << "Type three space-separated numbers for the bass/mid/treble levels.\n"
				  << "The values should be between 0 and 1.\n"
				  << "Don't forget to turn up the volume!\n"
				  << "You can quit by typing 'exit'.\n"
				  << std::endl;

		while (true) {
			std::cout << "EQ parameters: ";
			std::string userInput;
			std::getline(std::cin, userInput);
			if (userInput == "exit") {
				break;
			}

			std::stringstream ss(userInput);
			float bass = 1.0f;
			float mid = 1.0f;
			float treble = 1.0f;
			ss >> bass;
			ss >> mid;
			ss >> treble;

			Equalizer equalizer{ 513, sampleRate };
			equalizer.SetLevels(bass, mid, treble);

			size_t currentSample = 0;
			PlayStereo(sampleRate, [&left, &right, &currentSample, &equalizer](SignalView<float> leftOut, SignalView<float> rightOut) -> size_t {
				assert(leftOut.Size() == rightOut.Size());
				const size_t validSize = std::min(left.Size() - currentSample, leftOut.Size());
				equalizer.Filter(AsView(left).SubSignal(currentSample, validSize),
								 AsView(right).SubSignal(currentSample, validSize),
								 leftOut.SubSignal(0, validSize),
								 rightOut.SubSignal(0, validSize));
				currentSample += validSize;
				return validSize;
			});
		}
	}
	catch (std::exception& ex) {
		std::cout << ex.what() << std::endl;
	}
}