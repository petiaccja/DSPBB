#pragma once

#include "../Math/Statistics.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "FIR.hpp"


namespace dspbb {


template <class T, eSignalDomain Domain>
class PolyphaseView {
public:
	PolyphaseView() noexcept = default;
	PolyphaseView(BasicSignalView<T, Domain> data, size_t numFilters) noexcept
		: m_bufferView(data), m_filterCount(numFilters) {}

	BasicSignalView<T, Domain> operator[](size_t index) {
		auto loc = SubSignalLocation(index);
		return m_bufferView.subsignal(loc.first, loc.second);
	}
	BasicSignalView<const T, Domain> operator[](size_t index) const {
		auto loc = SubSignalLocation(index);
		return m_bufferView.subsignal(loc.first, loc.second);
	}
	size_t PhaseSize() const noexcept {
		return (m_bufferView.size() + m_filterCount - 1) / m_filterCount;
	}
	size_t OriginalSize() const noexcept {
		return m_bufferView.size();
	}
	size_t FilterCount() const noexcept {
		return m_filterCount;
	}

private:
	std::pair<size_t, size_t> SubSignalLocation(size_t index) const {
		assert(index < m_filterCount);
		const size_t numExtended = m_bufferView.size() % m_filterCount;
		const size_t baseFilterSize = m_bufferView.size() / m_filterCount;
		const size_t thisFilterSize = baseFilterSize + size_t(index < numExtended);
		const size_t offset = baseFilterSize * index + std::min(numExtended, index);
		return { offset, thisFilterSize };
	}

private:
	BasicSignalView<T, Domain> m_bufferView;
	size_t m_filterCount = 1;
};

template <class T, eSignalDomain Domain>
class PolyphaseFilter : public PolyphaseView<T, Domain> {
public:
	PolyphaseFilter() = default;
	PolyphaseFilter(size_t hrFilterSize, size_t numPhases) : m_buffer(hrFilterSize) {
		PolyphaseView<T, Domain>::operator=({ AsView(m_buffer), numPhases });
	}
	PolyphaseFilter(const PolyphaseFilter& rhs) : m_buffer(rhs.m_buffer) {
		PolyphaseView<T, Domain>::operator=({ AsView(m_buffer), rhs.FilterCount() });
	}
	PolyphaseFilter(PolyphaseFilter&& rhs) noexcept : m_buffer(std::move(rhs.m_buffer)) {
		PolyphaseView<T, Domain>::operator=(rhs);
		rhs.PolyphaseView<T, Domain>::operator=({});
	}
	PolyphaseFilter& operator=(const PolyphaseFilter& rhs) {
		m_buffer = rhs.m_buffer;
		PolyphaseView<T, Domain>::operator=({ AsView(m_buffer), rhs.FilterCount() });
		return *this;
	}
	PolyphaseFilter& operator=(PolyphaseFilter&& rhs) noexcept {
		m_buffer = std::move(rhs.m_buffer);
		PolyphaseView<T, Domain>::operator=(rhs);
		rhs.PolyphaseView<T, Domain>::operator=({});
		return *this;
	}
	~PolyphaseFilter() = default;

	BasicSignalView<T, Domain> Buffer() {
		return AsView(m_buffer);
	}
	BasicSignalView<const T, Domain> Buffer() const {
		return AsView(m_buffer);
	}

private:
	BasicSignal<T, Domain> m_buffer;
};

template <class T, eSignalDomain Domain>
void PolyphaseNormalize(PolyphaseView<T, Domain>& polyphase) {
	for (size_t i = 0; i < polyphase.FilterCount(); ++i) {
		polyphase[i] *= T(1) / Sum(polyphase[i]);
	}
}

template <class T, eSignalDomain Domain>
auto PolyphaseNormalized(PolyphaseFilter<T, Domain> polyphase) {
	PolyphaseNormalize(polyphase);
	return polyphase;
}

template <class SignalR, class SignalT, std::enable_if_t<is_same_domain_v<SignalR, SignalT> && is_mutable_signal_v<SignalR>, int> = 0>
auto PolyphaseDecompose(SignalR&& output, const SignalT& filter, size_t numFilters) {
	assert(output.size() == filter.size());
	assert(output.data() != filter.data());

	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	constexpr auto Domain = signal_traits<std::decay_t<SignalR>>::domain;
	PolyphaseView<R, Domain> view{ AsView(output), numFilters };

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

template <class SignalT>
auto PolyphaseDecompose(const SignalT& filter, size_t numFilters) {
	using R = typename signal_traits<std::decay_t<SignalT>>::type;
	constexpr auto Domain = signal_traits<std::decay_t<SignalT>>::domain;
	PolyphaseFilter<R, Domain> polyphase{ filter.size(), numFilters };

	PolyphaseDecompose(polyphase.Buffer(), filter, numFilters);

	return polyphase;
}


} // namespace dspbb