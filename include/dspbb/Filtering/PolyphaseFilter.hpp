#pragma once

#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "Convolution.hpp"
#include "FIR.hpp"


namespace dspbb {


template <class T>
class PolyphaseFilter {
public:
	PolyphaseFilter(uint64_t sampleRate, T cutoffFrequency, unsigned numFilters, unsigned numTaps);
	PolyphaseFilter(TimeSignalView<const T> filter, unsigned numFilters);

	size_t NumFilters() const;
	SignalView<const T, TIME_DOMAIN> Filter(size_t i) const;
	size_t NumTaps() const;

	template <class PaddingMode>
	size_t operator()(SignalView<const T, eSignalDomain::TIME> input, SignalView<T, eSignalDomain::TIME> output, PaddingMode) const;
	template <class PaddingMode>
	TimeSignal<T> operator()(SignalView<const T, eSignalDomain::TIME> input, PaddingMode) const;

	static TimeSignal<T> FirLowPassWindowed(uint64_t sampleRateIn, float cutoffFrequency, unsigned numFilters, unsigned numTaps);

private:
	static std::vector<TimeSignal<T>> Split(TimeSignalView<const T> filter, unsigned numFilters);

private:
	std::vector<TimeSignal<T>> m_filterBank;
};


template <class T>
PolyphaseFilter<T>::PolyphaseFilter(uint64_t sampleRate, T cutoffFrequency, unsigned numFilters, unsigned numTaps) {
	const auto lowPass = FirLowPassWindowed(sampleRate, cutoffFrequency, numFilters, numTaps * numFilters);
	m_filterBank = Split(AsConstView(lowPass), numFilters);
}

template <class T>
PolyphaseFilter<T>::PolyphaseFilter(TimeSignalView<const T> filter, unsigned numFilters) {
	m_filterBank = Split(filter, numFilters);
}

template <class T>
size_t PolyphaseFilter<T>::NumFilters() const {
	return m_filterBank.size();
}

template <class T>
SignalView<const T, TIME_DOMAIN> PolyphaseFilter<T>::Filter(size_t i) const {
	return { m_filterBank[i].begin(), m_filterBank[i].end() };
}

template <class T>
size_t PolyphaseFilter<T>::NumTaps() const {
	return !m_filterBank.empty() ? m_filterBank[0].Size() : 0;
}

template <class T>
TimeSignal<T> PolyphaseFilter<T>::FirLowPassWindowed(uint64_t sampleRateIn, float cutoffFrequency, unsigned numFilters, unsigned numTaps) {
	return dspbb::FirLowPassWindowed<T>(sampleRateIn * numFilters, cutoffFrequency, numTaps);
}

template <class T>
std::vector<TimeSignal<T>> PolyphaseFilter<T>::Split(TimeSignalView<const T> filter, unsigned numFilters) {
	std::vector<TimeSignal<T>> filterBank(numFilters);
	size_t i = 0;
	while (i < filter.Size()) {
		for (auto& filterPhase : filterBank) {
			filterPhase.PushBack(i < filter.Size() ? numFilters * filter[i] : T(0));
			++i;
		}
	}
	return filterBank;
}

template <class T>
template <class PaddingMode>
size_t PolyphaseFilter<T>::operator()(SignalView<const T, eSignalDomain::TIME> input, SignalView<T, eSignalDomain::TIME> output, PaddingMode) const {
	size_t offset = 0;
	for (auto& filter : m_filterBank) {
		auto filtered = Convolution(input, AsConstView(filter), PaddingMode{});
		size_t outIndex = offset;
		for (auto& v : filtered) {
			output[outIndex] = v;
			outIndex += m_filterBank.size();
		}
		++offset;
	}
	return 0;
}

template <class T>
template <class PaddingMode>
TimeSignal<T> PolyphaseFilter<T>::operator()(SignalView<const T, eSignalDomain::TIME> input, PaddingMode) const {
	TimeSignal<T> output(ConvolutionLength(input.Size(), m_filterBank[0].Size(), PaddingMode{}) * m_filterBank.size(), 0.0f);
	operator()(input, output, PaddingMode{});
	return output;
}


} // namespace dspbb