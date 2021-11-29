#pragma once

#include "../Primitives/Signal.hpp"
#include "../Utility/Numbers.hpp"
#include "../Utility/TypeTraits.hpp"

#include <cstdint>
#include <type_traits>


namespace dspbb {



namespace impl {

	template <class SignalR, class WaveFunc, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
	void GenericWave(SignalR&& output, uint64_t sampleRate, double frequency, double phase, WaveFunc waveFunc) {
		using R = typename signal_traits<std::decay_t<SignalR>>::type;
		using T = remove_complex_t<R>;
		size_t idx = 0;
		for (auto& v : output) {
			const double time = double(idx) / double(sampleRate);
			const double totalPhase = 2.0 * pi_v<double> * time * frequency + phase;
			v = R(T(waveFunc(totalPhase)));
			++idx;
		}
	}

	template <class SignalR, class WaveFunc, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
	void GenericChirp(SignalR&& output, uint64_t sampleRate, double startFrequency, double endFrequency, double phase, WaveFunc waveFunc) {
		using R = typename signal_traits<std::decay_t<SignalR>>::type;
		using T = remove_complex_t<R>;
		const double length = double(output.Size()) / double(sampleRate);
		size_t idx = 0;
		for (auto& v : output) {
			const double time = double(idx) / double(sampleRate);
			// Integrate the linear function frequency(time).
			const double totalPhase = 2.0 * pi_v<double> * (time * startFrequency + time * time / 2.0 * (endFrequency - startFrequency) / length) + phase;
			v = R(T(waveFunc(totalPhase)));
			++idx;
		}
	}

	inline double Sawtooth(double phase, double tilt) {
		constexpr auto pi2 = 2.0 * pi_v<double>;
		double intpart;
		const double unitPhase = modf(phase / pi2, &intpart);
		const double length = unitPhase > tilt ? (1.0 - tilt) : tilt + std::numeric_limits<double>::denorm_min();
		const double distance = std::abs(unitPhase - tilt);
		const double value = 1.0 - distance / length;
		return 2.0 * value - 1.0;
	}

	inline double Pwm(double phase, double fill) {
		constexpr auto pi2 = 2.0 * pi_v<double>;
		double intpart;
		const double unitPhase = modf(phase / pi2, &intpart);
		return double(unitPhase < fill || fill >= 1.0);
	}

} // namespace impl

//------------------------------------------------------------------------------
// Constant tone
//------------------------------------------------------------------------------

template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void SineWave(SignalR&& output, uint64_t sampleRate, double frequency, double phase = 0) {
	impl::GenericWave(output, sampleRate, frequency, phase, [](const auto& arg) { return std::sin(arg); });
}

template <class T, eSignalDomain Domain>
BasicSignal<T, Domain> SineWave(size_t length, uint64_t sampleRate, double frequency, double phase = 0) {
	BasicSignal<T, Domain> signal(length);
	SineWave(signal, sampleRate, frequency, phase);
	return signal;
}


template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void SawtoothWave(SignalR&& output, uint64_t sampleRate, double frequency, double phase = 0, double tilt = 1.0) {
	impl::GenericWave(output, sampleRate, frequency, phase, [tilt](const auto& arg) { return impl::Sawtooth(arg, tilt); });
}

template <class T, eSignalDomain Domain>
BasicSignal<T, Domain> SawtoothWave(size_t length, uint64_t sampleRate, double frequency, double phase = 0, double tilt = 1.0) {
	BasicSignal<T, Domain> signal(length);
	SawtoothWave(signal, sampleRate, frequency, phase, tilt);
	return signal;
}


template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void PwmWave(SignalR&& output, uint64_t sampleRate, double frequency, double phase = 0, double fill = 0.5) {
	impl::GenericWave(output, sampleRate, frequency, phase, [fill](const auto& arg) { return impl::Pwm(arg, fill); });
}

template <class T, eSignalDomain Domain>
BasicSignal<T, Domain> PwmWave(size_t length, uint64_t sampleRate, double frequency, double phase = 0, double fill = 0.5) {
	BasicSignal<T, Domain> signal(length);
	PwmWave(signal, sampleRate, frequency, phase, fill);
	return signal;
}


template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void SquareWave(SignalR&& output, uint64_t sampleRate, double frequency, double phase = 0) {
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	impl::GenericWave(output, sampleRate, frequency, phase, [](const auto& arg) { return impl::Pwm(arg, 0.5f); });
	output *= R(2.0);
	output -= R(1.0);
}

template <class T, eSignalDomain Domain>
BasicSignal<T, Domain> SquareWave(size_t length, uint64_t sampleRate, double frequency, double phase = 0) {
	BasicSignal<T, Domain> signal(length);
	SquareWave(signal, sampleRate, frequency, phase);
	return signal;
}


//------------------------------------------------------------------------------
// Chirp
//------------------------------------------------------------------------------

template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void SineChirp(SignalR&& output, uint64_t sampleRate, double startFrequency, double endFrequency, double phase = 0) {
	impl::GenericChirp(output, sampleRate, startFrequency, endFrequency, phase, [](const auto& arg) { return std::sin(arg); });
}

template <class T, eSignalDomain Domain>
BasicSignal<T, Domain> SineChirp(size_t length, uint64_t sampleRate, double startFrequency, double endFrequency, double phase = 0) {
	BasicSignal<T, Domain> signal(length);
	SineChirp(signal, sampleRate, startFrequency, endFrequency, phase);
	return signal;
}


template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void SawtoothChirp(SignalR&& output, uint64_t sampleRate, double startFrequency, double endFrequency, double phase = 0, double tilt = 1.0) {
	impl::GenericChirp(output, sampleRate, startFrequency, endFrequency, phase, [tilt](const auto& arg) { return impl::Sawtooth(arg, tilt); });
}

template <class T, eSignalDomain Domain>
BasicSignal<T, Domain> SawtoothChirp(size_t length, uint64_t sampleRate, double startFrequency, double endFrequency, double phase = 0, double tilt = 1.0) {
	BasicSignal<T, Domain> signal(length);
	SawtoothChirp(signal, sampleRate, startFrequency, endFrequency, phase, tilt);
	return signal;
}


template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void PwmChirp(SignalR&& output, uint64_t sampleRate, double startFrequency, double endFrequency, double phase = 0, double fill = 0.5) {
	impl::GenericChirp(output, sampleRate, startFrequency, endFrequency, phase, [fill](const auto& arg) { return impl::Pwm(arg, fill); });
}

template <class T, eSignalDomain Domain>
BasicSignal<T, Domain> PwmChirp(size_t length, uint64_t sampleRate, double startFrequency, double endFrequency, double phase = 0, double fill = 0.5) {
	BasicSignal<T, Domain> signal(length);
	PwmChirp(signal, sampleRate, startFrequency, endFrequency, phase, fill);
	return signal;
}


template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void SquareChirp(SignalR&& output, uint64_t sampleRate, double startFrequency, double endFrequency, double phase = 0) {
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	impl::GenericChirp(output, sampleRate, startFrequency, endFrequency, phase, [](const auto& arg) { return impl::Pwm(arg, 0.5f); });
	output *= R(2.0);
	output -= R(1.0);
}

template <class T, eSignalDomain Domain>
BasicSignal<T, Domain> SquareChirp(size_t length, uint64_t sampleRate, double startFrequency, double endFrequency, double phase = 0) {
	BasicSignal<T, Domain> signal(length);
	SquareChirp(signal, sampleRate, startFrequency, endFrequency, phase);
	return signal;
}



} // namespace dspbb
