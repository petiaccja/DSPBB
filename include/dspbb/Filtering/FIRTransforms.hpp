#pragma once


#include "../Primitives/SignalTraits.hpp"


namespace dspbb {


template <class SignalR, class SignalT, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT>, int> = 0>
void MirrorResponse(SignalR&& mirrored, const SignalT& filter) {
	assert(mirrored.Size() == filter.Size());
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	T sign = T(1);
	for (size_t i = 0; i < filter.Size(); ++i, sign *= T(-1)) {
		mirrored[i] = R(sign * filter[i]);
	}
}

template <class SignalR, class SignalT, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT>, int> = 0>
void ComplementaryResponse(SignalR&& complementary, const SignalT& filter) {
	assert(filter.Size() % 2 == 1);
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	Multiply(complementary, filter, T(-1));
	complementary[complementary.Size() / 2] += R(1);
}

template <class SignalR, class SignalT, class U, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT>, int> = 0>
void ShiftResponse(SignalR&& moved, const SignalT& filter, U normalizedFrequency) {
	U offset = static_cast<U>(filter.Size() / 2);
	U scale = pi_v<U> * normalizedFrequency;
	const size_t size = filter.Size();
	for (size_t i = 0; i < size / 2; ++i) {
		const U x = (U(i) - offset) * scale;
		const U c = std::cos(x);
		moved[i] = c * filter[i];
		moved[size - i - 1] = c * filter[size - i - 1];
	}
	moved *= typename signal_traits<SignalT>::type(2);
}


} // namespace dspbb