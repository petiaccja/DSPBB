#pragma once

#include <sstream>
#include <stdexcept>

namespace dspbb {

template <class T, class U>
T NormalizedFrequency(T frequency, U sampleRate) {
	return T(2) * frequency / T(sampleRate);
}

namespace impl {
	template <class T>
	void ThrowIfNotNormalized(T frequency) {
		if (!(T(0) <= frequency && frequency <= T(1))) {
			std::stringstream msg;
			msg << "The frequency must be normalized to between 0 and 1. (You gave " << frequency << ".)";
			throw std::domain_error(msg.str());
		}
	}

	template <class T, class... Ts, std::enable_if_t<std::conjunction_v<std::is_same<T, Ts>...>, int> = 0>
	void ThrowIfNotSorted(T frequency, Ts... frequencies) {
		struct {
			explicit operator bool() const {
				return result;
			}
			auto operator<(T rhs) {
				result = result && current <= rhs;
				current = rhs;
				return *this;
			}
			T current;
			bool result = true;
		} chain{ frequency };
		if (!(chain < ... < frequencies)) {
			std::stringstream msg;
			msg << "The frequencies must be in increasing order. (You gave ";
			msg << frequency << ",";
			(msg << ... << (std::to_string(frequencies) + ','));
			msg << ".)";
			throw std::invalid_argument(msg.str());
		}
	}
} // namespace impl


} // namespace dspbb