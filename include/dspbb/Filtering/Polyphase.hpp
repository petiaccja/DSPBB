#pragma once

#include "../Math/Statistics.hpp"
#include "../Signal/Signal.hpp"
#include "../Signal/SignalView.hpp"
#include "FIR.hpp"
#include "Polyphase.hpp"

#include <ranges>


namespace dspbb {


namespace impl {

	/// <summary> Return a view to the nth phase of a reordered polyphase filter. </summary>
	template <class SignalT>
		requires signal_or_view<std::remove_cvref_t<SignalT>>
	auto GetPolyphasePhase(SignalT&& data, size_t numPhases, size_t index) {
		assert(index < numPhases);
		const size_t numExtended = data.size() % numPhases;
		const size_t baseFilterSize = data.size() / numPhases;
		const size_t thisFilterSize = baseFilterSize + size_t(index < numExtended);
		const size_t offset = baseFilterSize * index + std::min(numExtended, index);
		return BasicSignalView(data).subsignal(offset, thisFilterSize);
	}

} // namespace impl


/// <summary> A polyphase-ordered filter with methods to access individual phases. </summary>
template <class T, eSignalDomain Domain>
class PolyphaseFilter {
public:
	using container = BasicSignal<T, Domain>;
	using value_type = typename container::value_type;

	PolyphaseFilter() = default;
	PolyphaseFilter(const PolyphaseFilter& rhs) noexcept = default;
	PolyphaseFilter(PolyphaseFilter&& rhs) noexcept = default;
	PolyphaseFilter& operator=(const PolyphaseFilter& rhs) noexcept = default;
	PolyphaseFilter& operator=(PolyphaseFilter&& rhs) noexcept = default;

	/// <summary> Construct a polyphase filter with uninitialized coefficients. </summary>
	/// <param name="originalSize"> The size of the original filter. </param>
	/// <param name="numPhases"> The number of phases. </param>
	PolyphaseFilter(size_t originalSize, size_t numPhases)
		: m_container(originalSize), m_numPhases(numPhases) {}

	/// <summary> Construct a polyphase filter with coefficients. </summary>
	/// <param name="coefficients"> The coefficients of the filter, must be in polyphase order. </param>
	/// <param name="numPhases"> The number of phases. </param>
	PolyphaseFilter(container coefficients, size_t numPhases) noexcept
		: m_container(std::move(coefficients)), m_numPhases(numPhases) {}


	/// <summary> Get the <paramref name="index"/>th phase of the filter. </summary>
	BasicSignalView<T, Domain> get_phase(size_t index) noexcept {
		return impl::GetPolyphasePhase(m_container, m_numPhases, index);
	}

	/// <summary> Get the <paramref name="index"/>th phase of the filter. </summary>
	BasicSignalView<T, Domain> operator[](size_t index) noexcept {
		return get_phase(index);
	}

	/// <summary> Get the <paramref name="index"/>th phase of the filter. </summary>
	BasicSignalView<const T, Domain> get_phase(size_t index) const noexcept {
		return impl::GetPolyphasePhase(m_container, m_numPhases, index);
	}

	/// <summary> Get the <paramref name="index"/>th phase of the filter. </summary>
	BasicSignalView<const T, Domain> operator[](size_t index) const noexcept {
		return get_phase(index);
	}

	/// <summary> Get the number of filter coefficients in a single phase. </summary>
	/// <remarks> Returns the size of the longest phase. </remarks>
	size_t size_per_phase() const noexcept {
		return (m_container.size() + m_numPhases - 1) / m_numPhases;
	}

	/// <summary> Get the size of the original filter, before polyphase reordering. </summary>
	size_t size_original() const noexcept {
		return m_container.size();
	}

	/// <summary> Get the number of phases. </summary>
	size_t num_phases() const noexcept {
		return m_numPhases;
	}

	/// <summary> Get the underlying container. </summary>
	/// <remarks> This view contains the coefficients in polyphase order. </remarks>
	auto get_container() noexcept {
		return BasicSignalView(m_container);
	}

	/// <summary> Get the underlying container. </summary>
	/// <remarks> This view contains the coefficients in polyphase order. </remarks>
	const container& get_container() const noexcept {
		return m_container;
	}

private:
	container m_container;
	size_t m_numPhases = 1;
};


/// <summary> A view into a polyphase-ordered filter that helps access individual phases. </summary>
template <class T, eSignalDomain Domain>
class PolyphaseFilterView {
public:
	using container = BasicSignalView<T, Domain>;
	static constexpr bool is_const = container::is_const;
	using value_type = typename container::value_type;

	PolyphaseFilterView() noexcept = default;

	/// <summary> Construct a view to a polyphase filter. </summary>
	/// <param name="data"> The filter coefficient. Must already be in polyphase order. </param>
	/// <param name="numPhases"> The number of the phases of the polyphase filter. </param>
	/// <remarks> The data is not reordered, it must have already been reordered to be in
	///		polyphase order. </remarks>
	PolyphaseFilterView(container data, size_t numPhases) noexcept
		: m_container(data),
		  m_numPhases(numPhases) {}


	/// <summary> Construct view from a polyphase filter. </summary>
	template <std::same_as<value_type> U>
		requires is_const
	PolyphaseFilterView(const PolyphaseFilter<U, Domain>& rhs)
		: m_container(rhs.get_container()),
		  m_numPhases(rhs.num_phases()) {}


	/// <summary> Construct view from a polyphase filter. </summary>
	PolyphaseFilterView(PolyphaseFilter<value_type, Domain>& rhs)
		: m_container(rhs.get_container()),
		  m_numPhases(rhs.num_phases()){};


	/// <summary> Construct view from another view. </summary>
	template <std::same_as<value_type> U>
		requires is_const
	PolyphaseFilterView(const PolyphaseFilterView<U, Domain>& rhs)
		: m_container(rhs.get_container()),
		  m_numPhases(rhs.num_phases()) {}


	/// <summary> Get the <paramref name="index"/>th phase of the filter. </summary>
	BasicSignalView<T, Domain> get_phase(size_t index) const noexcept {
		return impl::GetPolyphasePhase(m_container, m_numPhases, index);
	}

	/// <summary> Get the <paramref name="index"/>th phase of the filter. </summary>
	BasicSignalView<T, Domain> operator[](size_t index) const noexcept {
		return get_phase(index);
	}

	/// <summary> Get the number of filter coefficients in a single phase. </summary>
	/// <remarks> Returns the size of the longest phase. </remarks>
	size_t size_per_phase() const noexcept {
		return (m_container.size() + m_numPhases - 1) / m_numPhases;
	}

	/// <summary> Get the size of the original filter, before polyphase reordering. </summary>
	size_t size_original() const noexcept {
		return m_container.size();
	}

	/// <summary> Get the number of phases. </summary>
	size_t num_phases() const noexcept {
		return m_numPhases;
	}

	/// <summary> Get the underlying container. </summary>
	/// <remarks> This view contains the coefficients in polyphase order. </remarks>
	const container& get_container() const noexcept {
		return m_container;
	}

private:
	container m_container;
	size_t m_numPhases = 1;
};


template <class T, eSignalDomain Domain>
PolyphaseFilterView(const PolyphaseFilter<T, Domain>&) -> PolyphaseFilterView<const T, Domain>;


template <class T, eSignalDomain Domain>
PolyphaseFilterView(PolyphaseFilter<T, Domain>&) -> PolyphaseFilterView<T, Domain>;


template <class T, eSignalDomain Domain>
PolyphaseFilterView(const PolyphaseFilterView<T, Domain>&) -> PolyphaseFilterView<T, Domain>;


template <class T>
struct is_polyphase_or_view : std::false_type {};


template <class T, eSignalDomain Domain>
struct is_polyphase_or_view<PolyphaseFilter<T, Domain>> : std::true_type {};


template <class T, eSignalDomain Domain>
struct is_polyphase_or_view<PolyphaseFilterView<T, Domain>> : std::true_type {};


template <class T>
inline constexpr bool is_polyphase_or_view_v = is_polyphase_or_view<T>::value;


template <class T>
concept polyphase_or_view = is_polyphase_or_view_v<T>;


template <class T, eSignalDomain Domain>
struct is_mutable<PolyphaseFilter<T, Domain>> {
	static constexpr bool value = true;
};


template <class T, eSignalDomain Domain>
struct is_mutable<const PolyphaseFilter<T, Domain>> {
	static constexpr bool value = false;
};


template <class T, eSignalDomain Domain>
struct is_mutable<PolyphaseFilterView<T, Domain>> {
	static constexpr bool value = !std::is_const_v<T>;
};


template <class T>
concept mutable_polyphase_or_view = polyphase_or_view<T> && is_mutable_v<T>;


template <class T>
concept mutable_polyphase_or_view_r = polyphase_or_view<std::remove_reference_t<T>> && is_mutable_v<std::remove_reference_t<T>>;


/// <summary> Normalize each phase of a polyphase filter. </summary>
template <mutable_polyphase_or_view_r PolyphaseTy>
void PolyphaseNormalize(PolyphaseTy&& filter) noexcept {
	using T = typename std::decay_t<PolyphaseTy>::value_type;
	for (size_t i = 0; i < filter.num_phases(); ++i) {
		filter[i] *= T(1) / Sum(filter[i]);
	}
}


/// <summary> Normalize each phase of a polyphase filter. </summary>
/// <returns> The normalized polyphase filter. </returns>
template <polyphase_or_view PolyphaseTy>
auto PolyphaseNormalized(const PolyphaseTy& filter) {
	using T = typename PolyphaseTy::value_type;
	constexpr auto Domain = domain_v<typename PolyphaseTy::container>;
	PolyphaseFilter<T, Domain> normalized({ filter.get_container().begin(), filter.get_container().end() }, filter.num_phases());
	PolyphaseNormalize(normalized);
	return normalized;
}


/// <summary> Reorder the elements of a filter into polyphase order. </summary>
/// <param name="output"> The reordered output. </param>
/// <param name="filter"> The filter to reorder. </param>
/// <param name="numPhases"> The number of phases of the polyphase filter. </param>
template <mutable_polyphase_or_view_r PolyphaseTy, same_domain_as_r<typename std::decay_t<PolyphaseTy>::container> SignalT>
void PolyphaseReorder(PolyphaseTy&& output, const SignalT& filter) {
	assert(output.get_container().size() == filter.size());
	assert(!IsAliasing(output.get_container(), filter));

	using T = typename std::decay_t<PolyphaseTy>::container::value_type;

	for (size_t phaseIdx = 0; phaseIdx < output.num_phases(); ++phaseIdx) {
		const auto phase = std::ranges::reverse_view(output[phaseIdx]);
		for (size_t coeffIdx = 0; coeffIdx < phase.size(); ++coeffIdx) {
			phase[coeffIdx] = filter[coeffIdx * output.num_phases() + phaseIdx];
		}
	}
	output.get_container() *= T(output.num_phases());
}


/// <summary> Reorder the elements of a filter into polyphase order. </summary>
/// <param name="filter"> The filter to reorder. </param>
/// <param name="numPhases"> The number of phases of the polyphase filter. </param>
/// <returns> The reordered filter. </returns>
template <class SignalT>
auto PolyphaseReorder(const SignalT& filter, size_t numPhases) {
	using R = scalar_type_t<std::decay_t<SignalT>>;
	constexpr auto Domain = domain_v<std::decay_t<SignalT>>;
	PolyphaseFilter<R, Domain> polyphase{ filter.size(), numPhases };

	PolyphaseReorder(polyphase, filter);

	return polyphase;
}

} // namespace dspbb