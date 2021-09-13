#pragma once

#include "../Math/Statistics.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "FIR.hpp"


namespace dspbb {


template <class T, eSignalDomain Domain>
struct PolyphaseDecomposition {
	SignalView<T, Domain> filterBank;
	size_t numFilters;
	SignalView<T, Domain> operator[](size_t index) {
		auto loc = SubSignalLocation(index);
		return filterBank.SubSignal(loc.first, loc.second);
	}
	SignalView<const T, Domain> operator[](size_t index) const {
		auto loc = SubSignalLocation(index);
		return filterBank.SubSignal(loc.first, loc.second);
	}
	size_t Size() const {
		return (filterBank.Size() + numFilters - 1) / numFilters;
	}

private:
	std::pair<size_t, size_t> SubSignalLocation(size_t index) const {
		assert(index < numFilters);
		const size_t numExtended = filterBank.Size() % numFilters;
		const size_t baseFilterSize = filterBank.Size() / numFilters;
		const size_t thisFilterSize = baseFilterSize + size_t(index < numExtended);
		const size_t offset = baseFilterSize * index + std::min(numExtended, index);
		return { offset, thisFilterSize };
	}
};

template <class T, eSignalDomain Domain>
void PolyphaseNormalize(PolyphaseDecomposition<T, Domain>& polyphase) {
	for (size_t i = 0; i < polyphase.numFilters; ++i) {
		polyphase[i] *= T(1) / Sum(polyphase[i]);
	}
}

template <class T, eSignalDomain Domain>
auto PolyphaseNormalized(PolyphaseDecomposition<T, Domain>&& polyphase) {
	PolyphaseNormalize(polyphase);
	return polyphase;
}

template <class SignalR, class SignalT, std::enable_if_t<is_same_domain_v<SignalR, SignalT> && is_mutable_signal_v<SignalR>, int> = 0>
auto PolyphaseDecompose(SignalR&& output, const SignalT& filter, size_t numFilters) {
	assert(output.Size() == filter.Size());
	assert(output.Data() != filter.Data());

	PolyphaseDecomposition<typename signal_traits<std::decay_t<SignalR>>::type, signal_traits<std::decay_t<SignalR>>::domain> view{
		AsView(output),
		numFilters,
	};

	for (size_t phaseIdx = 0; phaseIdx < numFilters; ++phaseIdx) {
		auto filterPhase = view[phaseIdx];
		size_t coeffIdx = 0;
		for (auto it = filterPhase.rbegin(); it != filterPhase.rend(); ++it) {
			*it = filter[phaseIdx + coeffIdx * numFilters];
			++coeffIdx;
		}
		filterPhase *= numFilters;
	}

	return view;
}


} // namespace dspbb