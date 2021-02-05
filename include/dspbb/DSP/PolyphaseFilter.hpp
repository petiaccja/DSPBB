#pragma once

#include "Convolution.hpp"
#include "FIR.hpp"

#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"


namespace dspbb {

template <class T>
class PolyphaseFilter {
public:
	PolyphaseFilter(uint64_t sampleRate, float cutoffFrequency, unsigned numFilters, unsigned numTaps);

	size_t NumFilters() const;
	SignalView<const T, TIME_DOMAIN> Filter(size_t i) const;
	size_t NumTaps() const;
	
	template <class PaddingMode>
	size_t operator()(SignalView<const T, eSignalDomain::TIME> input, SignalView<T, eSignalDomain::TIME> output, PaddingMode) const;
	template <class PaddingMode>
	TimeSignal<T> operator()(SignalView<const T, eSignalDomain::TIME> input, PaddingMode);

private:
	TimeSignal<T> CreateLowPass(uint64_t sampleRate, float cutoffFrequency, unsigned numFilters, unsigned numTaps);
	std::vector<TimeSignal<T>> Split(TimeSignal<T> filter, unsigned numFilters);

private:
	std::vector<TimeSignal<T>> m_filterBank;
};


template <class T>
PolyphaseFilter<T>::PolyphaseFilter(uint64_t sampleRate, float cutoffFrequency, unsigned numFilters, unsigned numTaps) {
	const auto lowPass = CreateLowPass(sampleRate, cutoffFrequency, numFilters, numTaps * numFilters);
	m_filterBank = Split(lowPass, numFilters);
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
TimeSignal<T> PolyphaseFilter<T>::CreateLowPass(uint64_t sampleRate, float cutoffFrequency, unsigned numFilters, unsigned numTaps) {
	return WindowedLowPass<T>(sampleRate * numFilters, cutoffFrequency, numTaps);
}

template <class T>
std::vector<TimeSignal<T>> PolyphaseFilter<T>::Split(TimeSignal<T> filter, unsigned numFilters) {
	std::vector<TimeSignal<T>> filterBank(numFilters);
	for (size_t i = 0; i < filter.Size(); ++i) {
		filterBank[i % numFilters].PushBack(*(filter.begin() + i) * T(numFilters));
	}
	std::reverse(filterBank.begin(), filterBank.end());
	return filterBank;
}

template <class T>
template <class PaddingMode>
size_t PolyphaseFilter<T>::operator()(SignalView<const T, eSignalDomain::TIME> input, SignalView<T, eSignalDomain::TIME> output, PaddingMode) const {
	size_t offset = 0;
	for (auto& filter : m_filterBank) {
		auto filtered = ConvolutionFast(input, SignalView<const T, TIME_DOMAIN>{ filter.begin(), filter.end() }, PaddingMode{});
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
TimeSignal<T> PolyphaseFilter<T>::operator()(SignalView<const T, eSignalDomain::TIME> input, PaddingMode) {
	TimeSignal<T> output(ConvolutionLength(input.Size(), m_filterBank[0].Size(), PaddingMode{}) * m_filterBank.size(), 0.0f);
	operator()(input, output, PaddingMode{});
	return output;
}


} // namespace dspbb