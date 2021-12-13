#pragma once

#include "../../LTISystems/Systems.hpp"
#include "../../Math/DotProduct.hpp"
#include "../../Primitives/Signal.hpp"

#include <algorithm>


namespace dspbb {

//------------------------------------------------------------------------------
// Direct form I
//------------------------------------------------------------------------------

template <class T>
class DirectFormI {
public:
	DirectFormI() = default;
	explicit DirectFormI(size_t order);

	void Order(size_t order);
	void Reset();
	size_t Order() const;

	template <class InputT, class SystemT, std::enable_if_t<std::is_convertible_v<InputT, T> && std::is_convertible_v<SystemT, T>, int> = 0>
	T Feed(const InputT& input, const DiscreteTransferFunction<SystemT>& sys);

private:
	BasicSignal<T, eSignalDomain::DOMAINLESS> recursiveState;
	BasicSignal<T, eSignalDomain::DOMAINLESS> forwardState;
};

template <class T>
DirectFormI<T>::DirectFormI(size_t order) {
	Order(order);
}

template <class T>
void DirectFormI<T>::Order(size_t order) {
	recursiveState.Resize(order, T(0));
	forwardState.Resize(order + 1, T(0));
}

template <class T>
void DirectFormI<T>::Reset() {
	std::fill(recursiveState.begin(), recursiveState.end(), T(0));
	std::fill(forwardState.begin(), forwardState.end(), T(0));
}

template <class T>
size_t DirectFormI<T>::Order() const {
	return recursiveState.Size();
}

template <class T>
template <class InputT, class SystemT, std::enable_if_t<std::is_convertible_v<InputT, T> && std::is_convertible_v<SystemT, T>, int>>
T DirectFormI<T>::Feed(const InputT& input, const DiscreteTransferFunction<SystemT>& sys) {
	assert(!forwardState.Empty() && Order() >= sys.Order());

	const auto fwFull = AsConstView(sys.numerator.Coefficients());
	const auto recFull = AsConstView(sys.denominator.Coefficients());
	const auto recSec = recFull.SubSignal(0, recFull.Size() - 1);

	const auto normalization = *recFull.rbegin();

	forwardState[0] = static_cast<T>(input);
	std::rotate(forwardState.begin(), ++forwardState.begin(), forwardState.end());

	const auto fwSum = DotProduct(AsView(forwardState).SubSignal(forwardState.Size() - fwFull.Size()), fwFull);
	const auto recSum = DotProduct(AsView(recursiveState).SubSignal(recursiveState.Size() - recSec.Size()), recSec);
	const auto out = (fwSum - recSum) / normalization;

	if (recursiveState.Size() > 0) {
		recursiveState[0] = static_cast<T>(out);
		std::rotate(recursiveState.begin(), ++recursiveState.begin(), recursiveState.end());
	}

	return static_cast<T>(out);
}

//------------------------------------------------------------------------------
// Direct form II
//------------------------------------------------------------------------------

template <class T>
class DirectFormII {
public:
	DirectFormII() = default;
	explicit DirectFormII(size_t order);

	void Order(size_t order);
	void Reset();
	size_t Order() const;

	template <class InputT, class SystemT, std::enable_if_t<std::is_convertible_v<InputT, T> && std::is_convertible_v<SystemT, T>, int> = 0>
	T Feed(const InputT& input, const DiscreteTransferFunction<SystemT>& sys);

private:
	BasicSignal<T, eSignalDomain::DOMAINLESS> m_state;
};

template <class T>
DirectFormII<T>::DirectFormII(size_t order) {
	m_state.Resize(order + 1, T(0));
}

template <class T>
void DirectFormII<T>::Order(size_t order) {
	m_state.Resize(order + 1, T(0));
}

template <class T>
void DirectFormII<T>::Reset() {
	std::fill(m_state.begin(), m_state.end(), T(0));
}

template <class T>
size_t DirectFormII<T>::Order() const {
	return !m_state.Empty() ? m_state.Size() - 1 : 0;
}

template <class T>
template <class InputT, class SystemT, std::enable_if_t<std::is_convertible_v<InputT, T> && std::is_convertible_v<SystemT, T>, int>>
T DirectFormII<T>::Feed(const InputT& input, const DiscreteTransferFunction<SystemT>& sys) {
	assert(!m_state.Empty() && Order() >= sys.Order());

	const auto fwFull = AsConstView(sys.numerator.Coefficients());
	const auto recFull = AsConstView(sys.denominator.Coefficients());
	const auto recSec = recFull.SubSignal(0, recFull.Size() - 1);

	const auto stateFwView = AsView(m_state).SubSignal(m_state.Size() - fwFull.Size());
	const auto stateRecView = AsView(m_state).SubSignal(m_state.Size() - recSec.Size());
	const auto normalization = *recFull.rbegin();

	const auto stateNext = input / normalization - DotProduct(recSec, stateRecView);
	m_state[0] = static_cast<T>(stateNext);
	std::rotate(m_state.begin(), ++m_state.begin(), m_state.end());

	return static_cast<T>(DotProduct(fwFull, stateFwView));
}

//------------------------------------------------------------------------------
// Cascaded form
//------------------------------------------------------------------------------

template <class T>
class CascadedForm {
public:
	CascadedForm() = default;
	explicit CascadedForm(size_t order);

	void Order(size_t order);
	void Reset();
	size_t Order() const;

	template <class InputT, class SystemT, std::enable_if_t<std::is_convertible_v<InputT, T> && std::is_convertible_v<SystemT, T>, int> = 0>
	T Feed(const InputT& input, const CascadedBiquad<SystemT>& sys);

private:
	using Section = std::array<T, 3>;
	std::vector<Section> m_sections;
};


template <class T>
CascadedForm<T>::CascadedForm(size_t order) {
	const size_t numSections = 1 + (order + 1) / 2;
	m_sections.resize(numSections, { T(0), T(0), T(0) });
}

template <class T>
void CascadedForm<T>::Order(size_t order) {
	const size_t numSections = 1 + (order + 1) / 2;
	m_sections.resize(numSections, { T(0), T(0), T(0) });
}

template <class T>
void CascadedForm<T>::Reset() {
	for (auto& section : m_sections) {
		section = { T(0), T(0), T(0) };
	}
}

template <class T>
size_t CascadedForm<T>::Order() const {
	return (std::max(size_t(1), m_sections.size()) - 1) * 2;
}

template <class T>
template <class InputT, class SystemT, std::enable_if_t<std::is_convertible_v<InputT, T> && std::is_convertible_v<SystemT, T>, int>>
T CascadedForm<T>::Feed(const InputT& input, const CascadedBiquad<SystemT>& sys) {
	assert(sys.sections.size() + 1 <= m_sections.size());
	auto output = static_cast<T>(input);
	for (size_t i = 0; i < m_sections.size(); ++i) {
		const BasicSignalView<T, DOMAINLESS> stateFwView{ m_sections[i].begin(), m_sections[i].end() };
		stateFwView[0] = output;
		std::rotate(stateFwView.begin(), stateFwView.begin() + 1, stateFwView.end());

		if (i < sys.sections.size()) {
			const BasicSignalView<T, DOMAINLESS> stateRecView{ m_sections[i + 1].begin() + 1, m_sections[i + 1].end() };
			BasicSignalView<const SystemT, DOMAINLESS> fwView{ sys.sections[i].numerator.begin(), sys.sections[i].numerator.end() };
			BasicSignalView<const SystemT, DOMAINLESS> recView{ sys.sections[i].denominator.begin(), sys.sections[i].denominator.end() };

			const auto fwSum = DotProduct(stateFwView, fwView);
			const auto recSum = DotProduct(stateRecView, recView);
			output = static_cast<T>(fwSum - recSum);
		}
	}
	return output;
}

} // namespace dspbb