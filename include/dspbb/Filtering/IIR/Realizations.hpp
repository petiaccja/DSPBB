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

	T Feed(const T& input, const DiscreteTransferFunction<T>& sys);
	template <class U>
	T Feed(const U& input, const DiscreteTransferFunction<T>& sys) { return Feed(static_cast<T>(input), sys); }
	template <class SignalT>
	T Feed(const T& input, const SignalT& forward, const SignalT& recursive);

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
T DirectFormI<T>::Feed(const T& input, const DiscreteTransferFunction<T>& sys) {
	return Feed(input, sys.numerator.Coefficients(), sys.denominator.Coefficients());
}

template <class T>
template <class SignalT>
T DirectFormI<T>::Feed(const T& input, const SignalT& forward, const SignalT& recursive) {
	assert(!forwardState.Empty() && Order() + 1 >= std::max(forward.Size(), recursive.Size()));
	using U = std::decay_t<typename SignalT::value_type>;

	const auto fwView = BasicSignalView<const U, eSignalDomain::DOMAINLESS>{ forward.begin(), forward.end() };
	const auto recView = BasicSignalView<const U, eSignalDomain::DOMAINLESS>{ recursive.begin(), recursive.end() - 1 };
	const auto normalization = *recursive.rbegin();

	forwardState[0] = input;
	std::rotate(forwardState.begin(), ++forwardState.begin(), forwardState.end());

	const auto fwSum = DotProduct(AsView(forwardState).SubSignal(forwardState.Size() - forward.Size()), fwView);
	const auto recSum = DotProduct(AsView(recursiveState).SubSignal(recursiveState.Size() - recursive.Size() + 1), recView);
	const auto out = (fwSum - recSum) / normalization;

	if (recursiveState.Size() > 0) {
		recursiveState[0] = out;
		std::rotate(recursiveState.begin(), ++recursiveState.begin(), recursiveState.end());
	}

	return out;
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

	T Feed(const T& input, const DiscreteTransferFunction<T>& sys);
	template <class U>
	T Feed(const U& input, const DiscreteTransferFunction<T>& sys) { return Feed(static_cast<T>(input), sys); }
	template <class SignalT>
	T Feed(const T& input, const SignalT& forward, const SignalT& recursive);

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
T DirectFormII<T>::Feed(const T& input, const DiscreteTransferFunction<T>& sys) {
	return Feed(input, sys.numerator.Coefficients(), sys.denominator.Coefficients());
}

template <class T>
template <class SignalT>
T DirectFormII<T>::Feed(const T& input, const SignalT& forward, const SignalT& recursive) {
	assert(!m_state.Empty() && Order() + 1 >= forward.Size());
	using U = std::decay_t<typename SignalT::value_type>;

	const auto fwView = BasicSignalView<const U, eSignalDomain::DOMAINLESS>{ forward.begin(), forward.end() };
	const auto recView = BasicSignalView<const U, eSignalDomain::DOMAINLESS>{ recursive.begin(), recursive.end() - 1 };
	const auto stateFwView = AsView(m_state).SubSignal(m_state.Size() - forward.Size());
	const auto stateRecView = AsView(m_state).SubSignal(m_state.Size() - recursive.Size() + 1);
	const auto normalization = *recursive.rbegin();

	const T stateNext = input / normalization - DotProduct(recView, stateRecView);
	m_state[0] = stateNext;
	std::rotate(m_state.begin(), ++m_state.begin(), m_state.end());

	return DotProduct(fwView, stateFwView);
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

	T Feed(const T& input, const CascadedBiquad<T>& sys);
	template <class U>
	T Feed(const U& input, const CascadedBiquad<T>& sys) { return Feed(static_cast<T>(input), sys); }

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
T CascadedForm<T>::Feed(const T& input, const CascadedBiquad<T>& sys) {
	assert(sys.sections.size() + 1 <= m_sections.size());
	auto output = input;
	for (size_t i = 0; i < m_sections.size(); ++i) {
		BasicSignalView<T, DOMAINLESS> stateFwView{ m_sections[i].begin(), m_sections[i].end() };
		stateFwView[0] = output;
		std::rotate(stateFwView.begin(), stateFwView.begin() + 1, stateFwView.end());

		if (i < sys.sections.size()) {
			BasicSignalView<T, DOMAINLESS> stateRecView{ m_sections[i + 1].begin() + 1, m_sections[i + 1].end() };
			BasicSignalView<const T, DOMAINLESS> fwView{ sys.sections[i].numerator.begin(), sys.sections[i].numerator.end() };
			BasicSignalView<const T, DOMAINLESS> recView{ sys.sections[i].denominator.begin(), sys.sections[i].denominator.end() };

			const auto fwSum = DotProduct(stateFwView, fwView);
			const auto recSum = DotProduct(stateRecView, recView);
			output = fwSum - recSum;
		}
	}
	return output;
}


} // namespace dspbb