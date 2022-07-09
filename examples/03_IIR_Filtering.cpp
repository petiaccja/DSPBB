//------------------------------------------------------------------------------
// 03. IIR filtering
//
// In this example, we will use infinite impulse response (IIR) filters to
// detect dialed phone numbers encoded with dual-tone multi-frequency (DTMF)
// signaling. In DTMF signaling, a sum of two signals of predefined frequencies
// are used to encode digits 0-9, the letters A-D and the cahracters # and *.
// Decoding could be done digitally by an FFT, but we will use IIR filters in the
// spirit of the analog age this technique comes from.
//------------------------------------------------------------------------------

#include "Tools/PlaySound.hpp"

#include <dspbb/Filtering/IIR.hpp>
#include <dspbb/Generators/Waveforms.hpp>
#include <dspbb/Math/Statistics.hpp>

#include <array>
#include <iostream>
#include <map>

using namespace dspbb;

// We will use a sample rate of 4000 Hz for the telephone line.
// That's enough to encode both speech and the DTMF signals.
constexpr unsigned sampleRate = 4000;

// The four possible predefined frequencies of the first tone.
const std::array<float, 4> frequencies1 = {
	697.f,
	770.f,
	852.f,
	941.f,
};

// The four possible predefined frequencies of the second tone.
const std::array<float, 4> frequencies2 = {
	1209.f,
	1336.f,
	1477.f,
	1633.f,
};

// Which tones are used to encode the given character in DTMF.
const std::map<char, std::pair<int, int>> characters = {
	{ '1', { 0, 0 } }, // The sum of tones 697 Hz & 1209 Hz encode digit '1'
	{ '2', { 0, 1 } },
	{ '3', { 0, 2 } },
	{ '4', { 1, 0 } },
	{ '5', { 1, 1 } },
	{ '6', { 1, 2 } },
	{ '7', { 2, 0 } },
	{ '8', { 2, 1 } },
	{ '9', { 2, 2 } },
	{ '0', { 3, 1 } },
	{ 'A', { 0, 3 } },
	{ 'B', { 1, 3 } },
	{ 'C', { 2, 3 } },
	{ 'D', { 3, 3 } },
	{ '#', { 3, 2 } },
	{ '*', { 3, 0 } },
};

// We will generate the dialed tones for this demo, but we could also load
// an audio clip from disk.
Signal<float> DialTone(char character) {
	const auto it = characters.find(character);
	if (it == characters.end()) {
		throw std::invalid_argument("Provide a valid DTMF character: 0-9 A-D #*");
	}

	const auto& toneIndices = it->second;
	const float frequency1 = frequencies1[toneIndices.first];
	const float frequency2 = frequencies2[toneIndices.second];

	// The dialed tone is simply the sum of sine waves of the frequencies
	// that encode the requested character.
	return SineWave<float, TIME_DOMAIN>(3000, sampleRate, frequency1)
		   + SineWave<float, TIME_DOMAIN>(3000, sampleRate, frequency2);
}

// We will need a filter bank, with one narrow bandpass filter tuned to
// each of the 8 DTMF signaling frequencies.
// To do it efficiently:
//  - We use an elliptic filter because we need a sharp transition, not smooth response and clean phase
//  - We set a loose pass-band ripple because the precise magnitude of the picked-up tone is not important
//  - We set a strict stop-band ripple to heavily suppress noise outside the narrow band
constexpr int filterOrder = 6;
const auto filterDesc = Iir.Bandpass.Elliptic.PassbandRipple(0.15f).StopbandRipple(0.02f);
float Normalize(float f) { return NormalizedFrequency(f, sampleRate); }; // A little helper since we need normalized frequencies.

// DesignFilter returns the zero-pole representation of the designed filter. To apply it
// to a signal, you have to convert it to a transfer function or a cascaded biquad.
// Unless you have a good reason, use cascaded biquads for their superior stability and accuracy.
const std::array<CascadedBiquad<float>, 4> filterBank1 = {
	CascadedBiquad{ DesignFilter<float>(filterOrder, filterDesc.Band(Normalize(frequencies1[0] - 10.f), Normalize(frequencies1[0] + 10.f))) },
	CascadedBiquad{ DesignFilter<float>(filterOrder, filterDesc.Band(Normalize(frequencies1[1] - 10.f), Normalize(frequencies1[1] + 10.f))) },
	CascadedBiquad{ DesignFilter<float>(filterOrder, filterDesc.Band(Normalize(frequencies1[2] - 10.f), Normalize(frequencies1[2] + 10.f))) },
	CascadedBiquad{ DesignFilter<float>(filterOrder, filterDesc.Band(Normalize(frequencies1[3] - 10.f), Normalize(frequencies1[3] + 10.f))) },
};

const std::array<CascadedBiquad<float>, 4> filterBank2 = {
	CascadedBiquad{ DesignFilter<float>(filterOrder, filterDesc.Band(Normalize(frequencies2[0] - 10.f), Normalize(frequencies2[0] + 10.f))) },
	CascadedBiquad{ DesignFilter<float>(filterOrder, filterDesc.Band(Normalize(frequencies2[1] - 10.f), Normalize(frequencies2[1] + 10.f))) },
	CascadedBiquad{ DesignFilter<float>(filterOrder, filterDesc.Band(Normalize(frequencies2[2] - 10.f), Normalize(frequencies2[2] + 10.f))) },
	CascadedBiquad{ DesignFilter<float>(filterOrder, filterDesc.Band(Normalize(frequencies2[3] - 10.f), Normalize(frequencies2[3] + 10.f))) },
};


// The detection routine runs the signal through every filter in the two filter banks,
// and if they let through any signal, the corresponding tone is present in the DTMF signal.
// If there are exactly two filters triggered, we have a valid DTMF character.
std::optional<char> Detect(SignalView<const float> signal) {
	// If a filter's output is stronger than the threshold, it's going to be considered
	// a present DTMF tone.
	const float threshold = 0.15f * RootMeanSquare(signal);

	std::array<bool, 4> containedTones1;
	std::array<bool, 4> containedTones2;

	// Even though we will run 8 separate filters, we can reuse the same memory for output.
	Signal<float> filtered(signal.size());

	// Due to the recursive nature of the IIR filters, they need a structure to store state.
	CascadedForm<float> state{ filterOrder };

	for (int i = 0; i < 4; ++i) {
		state.Reset(); // The state is reset not to carry over garbage from the previous filtering.
		Filter(filtered, signal, filterBank1[i], state);
		const float passed = RootMeanSquare(filtered);
		containedTones1[i] = passed > threshold;
	}

	for (int i = 0; i < 4; ++i) {
		state.Reset();
		Filter(filtered, signal, filterBank2[i], state);
		const float passed = RootMeanSquare(filtered);
		containedTones2[i] = passed > threshold;
	}

	// Valid DTMF signals contain exactly one of each key frequency sets.
	const size_t containedCount1 = std::count(containedTones1.begin(), containedTones1.end(), true);
	const size_t containedCount2 = std::count(containedTones2.begin(), containedTones2.end(), true);
	if (containedCount1 != 1 || containedCount2 != 1) {
		return {};
	}

	// Let's find the the signaled character using the character->frequency map.
	const int toneIndex1 = int(std::find(containedTones1.begin(), containedTones1.end(), true) - containedTones1.begin());
	const int toneIndex2 = int(std::find(containedTones2.begin(), containedTones2.end(), true) - containedTones2.begin());
	for (const auto& [character, indices] : characters) {
		if (indices == std::make_pair(toneIndex1, toneIndex2)) {
			return { character };
		}
	}
	return {};
}

// A simple loop that asks you to type a number to dial on the phone.
// The above functions are used to encode the dialed digit as a signal
// and then to decode the generated signal. It's pointless, but how
// else am I gonna make a demonstrative example?
int main() {
	std::cout << "Welcome to the IIR filtering example program.\n"
			  << "You can type the digits 0-9, A-D, * and # to dial on the phone.\n"
			  << "If you turn up the volume, you can also hear the dial tones played -- \n"
			  << "they should sound familiar.\n"
			  << "Type 'exit' to exit.\n"
			  << std::endl;

	std::string userInput;
	do {
		std::cout << "Enter a digit to dial: ";
		std::cin >> userInput;
		if (userInput.size() == 1) {
			Signal<float> signal;
			try {
				std::cout << "   Dialing..." << std::endl;
				signal = DialTone(userInput[0]);
				PlayMono(sampleRate, signal);
			}
			catch (std::exception&) {
				std::cout << "   Invalid digit, try again." << std::endl;
				continue;
			}

			std::cout << "   Detecting..." << std::endl;

			const auto detected = Detect(signal);
			if (detected) {
				std::cout << "   You dialed: " << detected.value() << std::endl;
			}
			else {
				std::cout << "   Dial tone does not represent any digit." << std::endl;
			}
		}
	} while (userInput != "exit");
}
