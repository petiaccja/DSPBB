#pragma once

#include <cstdint>
#include <type_traits>
#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/Utility/Numbers.hpp>
#include <dspbb/Utility/TypeTraits.hpp>


namespace dspbb {




namespace impl {

	template <class SignalR, class WaveFunc, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
	void Waveform(SignalR&& output, uint64_t sampleRate, double frequency, double phase, WaveFunc waveFunc) {
		using R = typename signal_traits<std::decay_t<SignalR>>::type;
		size_t idx = 0;
		for (auto& v : output) {
			double time = double(idx) / double(sampleRate);
			double totalPhase = 2.0 * pi_v<double> * time * double(frequency) + double(phase);
			v = R(remove_complex_t<R>(waveFunc(totalPhase)));
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


template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void SineWave(SignalR&& output, uint64_t sampleRate, double frequency, double phase = 0) {
	impl::Waveform(output, sampleRate, frequency, phase, [](const auto& arg) { return std::sin(arg); });
}

template <class T, eSignalDomain Domain>
Signal<T, Domain> SineWave(size_t length, uint64_t sampleRate, double frequency, double phase = 0) {
	Signal<T, Domain> signal(length);
	SineWave(signal, sampleRate, frequency, phase);
	return signal;
}


template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void SawtoothWave(SignalR&& output, uint64_t sampleRate, double frequency, double phase = 0, double tilt = 1.0) {
	impl::Waveform(output, sampleRate, frequency, phase, [tilt](const auto& arg) { return impl::Sawtooth(arg, tilt); });
}

template <class T, eSignalDomain Domain>
Signal<T, Domain> SawtoothWave(size_t length, uint64_t sampleRate, double frequency, double phase = 0, double tilt = 1.0) {
	Signal<T, Domain> signal(length);
	SawtoothWave(signal, sampleRate, frequency, phase, tilt);
	return signal;
}


template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void PwmWave(SignalR&& output, uint64_t sampleRate, double frequency, double phase = 0, double fill = 0.5) {
	impl::Waveform(output, sampleRate, frequency, phase, [fill](const auto& arg) { return impl::Pwm(arg, fill); });
}

template <class T, eSignalDomain Domain>
Signal<T, Domain> PwmWave(size_t length, uint64_t sampleRate, double frequency, double phase = 0, double fill = 0.5) {
	Signal<T, Domain> signal(length);
	PwmWave(signal, sampleRate, frequency, phase, fill);
	return signal;
}


template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void SquareWave(SignalR&& output, uint64_t sampleRate, double frequency, double phase = 0) {
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	impl::Waveform(output, sampleRate, frequency, phase, [](const auto& arg) { return impl::Pwm(arg, 0.5f); });
	output *= R(2.0);
	output -= R(1.0);
}

template <class T, eSignalDomain Domain>
Signal<T, Domain> SquareWave(size_t length, uint64_t sampleRate, double frequency, double phase = 0) {
	Signal<T, Domain> signal(length);
	SquareWave(signal, sampleRate, frequency, phase);
	return signal;
}


} // namespace dspbb
