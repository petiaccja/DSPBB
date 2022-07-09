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

	template <class InIter, class OutIter, class SystemT, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<InIter>()), T> && std::is_convertible_v<SystemT, T>, int> = 0>
	void Feed(InIter first, InIter last, OutIter outFirst, const DiscreteTransferFunction<SystemT>& sys);

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
	recursiveState.resize(order, T(0));
	forwardState.resize(order + 1, T(0));
}

template <class T>
void DirectFormI<T>::Reset() {
	std::fill(recursiveState.begin(), recursiveState.end(), T(0));
	std::fill(forwardState.begin(), forwardState.end(), T(0));
}

template <class T>
size_t DirectFormI<T>::Order() const {
	return recursiveState.size();
}

template <class T>
template <class InputT, class SystemT, std::enable_if_t<std::is_convertible_v<InputT, T> && std::is_convertible_v<SystemT, T>, int>>
T DirectFormI<T>::Feed(const InputT& input, const DiscreteTransferFunction<SystemT>& sys) {
	assert(!forwardState.empty() && Order() >= sys.Order());

	T output;
	Feed(&input, &input + 1, &output, sys);
	return output;
}

template <class T>
template <class InIter, class OutIter, class SystemT, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<InIter>()), T> && std::is_convertible_v<SystemT, T>, int>>
void DirectFormI<T>::Feed(InIter first, InIter last, OutIter outFirst, const DiscreteTransferFunction<SystemT>& sys) {
	assert(!forwardState.empty() && Order() >= sys.Order());

	const auto fwFull = AsConstView(sys.numerator.Coefficients());
	const auto recFull = AsConstView(sys.denominator.Coefficients());
	const auto recSec = recFull.subsignal(0, recFull.size() - 1);

	const auto fwStateView = AsView(forwardState).subsignal(forwardState.size() - fwFull.size());
	const auto recStateView = AsView(recursiveState).subsignal(recursiveState.size() - recSec.size());

	const auto normalization = T(1) / static_cast<T>(*recFull.rbegin());

	while (first != last) {
		const auto input = *first++;

		std::move(++forwardState.begin(), forwardState.end(), forwardState.begin());
		*forwardState.rbegin() = static_cast<T>(input);

		const auto fwSum = std::inner_product(fwStateView.begin(), fwStateView.end(), fwFull.begin(), T(0));
		const auto recSum = std::inner_product(recStateView.begin(), recStateView.end(), recSec.begin(), T(0));
		const auto out = (fwSum - recSum) * normalization;

		if (recursiveState.size() > 0) {
			std::move(++recursiveState.begin(), recursiveState.end(), recursiveState.begin());
			*recursiveState.rbegin() = static_cast<T>(out);
		}
		*outFirst++ = static_cast<T>(out);
	}
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

	template <class InIter, class OutIter, class SystemT, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<InIter>()), T> && std::is_convertible_v<SystemT, T>, int> = 0>
	void Feed(InIter first, InIter last, OutIter outFirst, const DiscreteTransferFunction<SystemT>& sys);

private:
	BasicSignal<T, eSignalDomain::DOMAINLESS> m_state;
};

template <class T>
DirectFormII<T>::DirectFormII(size_t order) {
	m_state.resize(order + 1, T(0));
}

template <class T>
void DirectFormII<T>::Order(size_t order) {
	m_state.resize(order + 1, T(0));
}

template <class T>
void DirectFormII<T>::Reset() {
	std::fill(m_state.begin(), m_state.end(), T(0));
}

template <class T>
size_t DirectFormII<T>::Order() const {
	return !m_state.empty() ? m_state.size() - 1 : 0;
}

template <class T>
template <class InputT, class SystemT, std::enable_if_t<std::is_convertible_v<InputT, T> && std::is_convertible_v<SystemT, T>, int>>
T DirectFormII<T>::Feed(const InputT& input, const DiscreteTransferFunction<SystemT>& sys) {
	assert(!m_state.empty() && Order() >= sys.Order());

	T output;
	Feed(&input, &input + 1, &output, sys);
	return output;
}

template <class T>
template <class InIter, class OutIter, class SystemT, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<InIter>()), T> && std::is_convertible_v<SystemT, T>, int>>
void DirectFormII<T>::Feed(InIter first, InIter last, OutIter outFirst, const DiscreteTransferFunction<SystemT>& sys) {
	assert(!m_state.empty() && Order() >= sys.Order());

	const auto fwFull = AsConstView(sys.numerator.Coefficients());
	const auto recFull = AsConstView(sys.denominator.Coefficients());
	const auto recSec = recFull.subsignal(0, recFull.size() - 1);

	const auto stateFwView = AsView(m_state).subsignal(m_state.size() - fwFull.size());
	const auto stateRecView = AsView(m_state).subsignal(m_state.size() - recSec.size());
	const auto normalization = T(1) / static_cast<T>(*recFull.rbegin());

	while (first != last) {
		const auto input = *first++;
		const auto stateNext = input * normalization - std::inner_product(recSec.begin(), recSec.end(), stateRecView.begin(), T(0));
		std::move(++m_state.begin(), m_state.end(), m_state.begin());
		*m_state.rbegin() = static_cast<T>(stateNext);
		const auto output = static_cast<T>(std::inner_product(fwFull.begin(), fwFull.end(), stateFwView.begin(), T(0)));
		*outFirst++ = output;
	}
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

	template <class InIter, class OutIter, class SystemT, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<InIter>()), T> && std::is_convertible_v<SystemT, T>, int> = 0>
	void Feed(InIter first, InIter last, OutIter outFirst, const CascadedBiquad<SystemT>& sys);

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
		auto& currentSection = m_sections[i];
		currentSection[0] = currentSection[1];
		currentSection[1] = currentSection[2];
		currentSection[2] = output;

		if (i < sys.sections.size()) {
			const auto& nextSection = m_sections[i + 1];
			const auto& sysSectionNum = sys.sections[i].numerator;
			const auto& sysSectionDen = sys.sections[i].denominator;

			const auto fwSum = currentSection[0] * sysSectionNum[0]
							   + currentSection[1] * sysSectionNum[1]
							   + currentSection[2] * sysSectionNum[2];
			const auto recSum = nextSection[1] * sysSectionDen[0]
								+ nextSection[2] * sysSectionDen[1];
			output = static_cast<T>(fwSum - recSum);
		}
	}
	return output;
}

template <class T>
template <class InIter, class OutIter, class SystemT, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<InIter>()), T> && std::is_convertible_v<SystemT, T>, int>>
void CascadedForm<T>::Feed(InIter first, InIter last, OutIter outFirst, const CascadedBiquad<SystemT>& sys) {
	while (first != last) {
		*outFirst++ = Feed(*first++, sys);
	}
}

} // namespace dspbb