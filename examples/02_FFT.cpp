//------------------------------------------------------------------------------
// 02. FFT
//
// This demo code uses the Fast Fourier Transform to measure the solar cycle.
// The number of sunspots on the surface of the Sun has been recorded since 1749,
// and is listed in a file in the data folder. The variation in the number of
// sunspots follow a cycle of ~11 years, and the periodicity should be easy to
// find by an FFT (or auto-correlation).
//------------------------------------------------------------------------------

#include <dspbb/Filtering/Windowing.hpp>
#include <dspbb/Math/FFT.hpp>
#include <dspbb/Primitives/Signal.hpp>

#include <fstream>
#include <iostream>

using namespace dspbb;

// Load the monthly average number of sunspots. There is one sample per month.
Signal<float> LoadSunspotHistory() {
	Signal<float> sunspotHistory;
	std::ifstream file(DSPBB_EXAMPLES_DATA R"(SN_m_tot_V2.0.txt)");
	while (file.good()) {
		float numSunspots = 0.0;
		file >> numSunspots;
		sunspotHistory.push_back(numSunspots);
	}
	return sunspotHistory;
}

int main() {
	// Since we have one sample per month, we assume a sampling rate of 12
	// so that we can easily work with periodicity in years.
	constexpr unsigned sampleRate = 12;

	// Load the time domain data on sunspot counts.
	const Signal<float> sunspotHistory = LoadSunspotHistory();

	// Apply the fourier transform to the time-domain data to reveal periodicity.
	const Signal<float> window = BlackmanHarrisWindow<float, TIME_DOMAIN>(sunspotHistory.size());
	const Spectrum<std::complex<float>> spectrum = Fft(sunspotHistory * window, FFT_HALF);
	const Spectrum<float> amplitude = Abs(spectrum);

	// Find the FFT bin with the highest amplitude. That will correspond to the frequency of the solar cycle.
	// Since we know the solar cycle's period is less than 100 years, we can exclude frequencies below
	// 0.01/year, thus also excluding the expected spike at DC.
	const size_t firstBin = FourierFrequency2Bin(1.0 / 100, sunspotHistory.size(), sampleRate);
	const auto maxBinIt = std::max_element(amplitude.begin() + firstBin, amplitude.end());
	const size_t maxBin = maxBinIt - amplitude.begin();
	const float solarCycleFrequency = (float)FourierBin2Frequency(maxBin, sunspotHistory.size(), sampleRate);
	const float solarCyclePeriod = 1.0f / solarCycleFrequency;

	// Should be about 11 years.
	std::cout << "Solar cycle: " << solarCyclePeriod << " years" << std::endl;
}