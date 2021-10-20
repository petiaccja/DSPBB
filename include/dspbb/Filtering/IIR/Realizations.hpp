#pragma once

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

	template <class SignalT>
	T Feed(const T& input, const SignalT& forward, const SignalT& recursive);

private:
	Signal<T, eSignalDomain::DOMAINLESS> recursiveState;
	Signal<T, eSignalDomain::DOMAINLESS> forwardState;
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
	return forwardState.Size();
}

template <class T>
template <class SignalT>
T DirectFormI<T>::Feed(const T& input, const SignalT& forward, const SignalT& recursive) {
	using U = std::decay_t<typename SignalT::value_type>;
	const auto normalization = *recursive.rbegin();

	forwardState[0] = input;
	std::rotate(forwardState.begin(), ++forwardState.begin(), forwardState.end());

	const auto fwView = SignalView<const U, eSignalDomain::DOMAINLESS>{ forward.begin(), forward.end() };
	const auto recView = SignalView<const U, eSignalDomain::DOMAINLESS>{ recursive.begin(), recursive.end() - 1 };
	const float fwSum = DotProduct(forwardState, fwView);
	const float recSum = DotProduct(recursiveState, recView);
	const float out = (fwSum - recSum) / normalization;

	if (recursiveState.Size() > 0) {
		recursiveState[0] = out;
		std::rotate(recursiveState.begin(), ++recursiveState.begin(), recursiveState.end());
	}

	return out;
}


//------------------------------------------------------------------------------
// Other forms
//------------------------------------------------------------------------------

template <class T>
struct DirectFormII {
	Signal<T, eSignalDomain::DOMAINLESS> state;
};

template <class T>
struct CascadedFormI {
	Signal<T, eSignalDomain::DOMAINLESS> feedback;
	Signal<T, eSignalDomain::DOMAINLESS> input;
};

template <class T>
struct CascadedFormII {
	Signal<T, eSignalDomain::DOMAINLESS> state;
};


} // namespace dspbb