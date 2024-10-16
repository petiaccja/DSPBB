#pragma once

#include "../Signal/Signal.hpp"
#include "../Utility/TypeTraits.hpp"

#include <cstdint>
#include <numbers>
#include <type_traits>



namespace dspbb {



namespace impl {

	template <mutable_signal_or_view_r SignalOut, class WaveFunc, std::floating_point Real>
	void GenericWave(SignalOut&& output, Real sampleRate, Real frequency, Real startPhase, WaveFunc waveFunc) {
		using Scalar = typename std::decay_t<SignalOut>::value_type;
		size_t idx = 0;
		for (auto& v : output) {
			const auto period = sampleRate / frequency;
			const auto basePhase = Real(2) * std::numbers::pi_v<Real> * std::fmod(Real(idx) / period, Real(1));
			const auto phase = basePhase + startPhase;
			v = static_cast<Scalar>(waveFunc(phase));
			++idx;
		}
	}

	template <mutable_signal_or_view_r SignalOut, class WaveFunc, std::floating_point Real>
	void GenericChirp(SignalOut&& output, Real sampleRate, Real startFrequency, Real endFrequency, Real startPhase, WaveFunc waveFunc) {
		// The frequency in function of time can be written as f(t) = (l - t) / l * f0 + t / l * f1.
		// With l = length of signal (seconds), t = time, (f0, f1) = start and end frequencies.
		// The phase phi(t) can be acquired by the definite integral of 2*pi*f(t') from 0 to t.

		using Scalar = typename std::decay_t<SignalOut>::value_type;
		const auto length = Real(output.size()) / sampleRate;
		size_t idx = 0;
		for (auto& v : output) {
			const auto time = Real(idx) / Real(sampleRate);
			const auto integratedTime = time * (startFrequency + time * (endFrequency - startFrequency) / (Real(2) * length));
			const auto basePhase = Real(2) * std::numbers::pi_v<Real> * integratedTime;
			const auto phase = basePhase + startPhase;
			v = static_cast<Scalar>(waveFunc(phase));
			++idx;
		}
	}

	template <std::floating_point Real>
	inline Real Sawtooth(Real phase, Real tilt) {
		constexpr auto pi2 = Real(2) * std::numbers::pi_v<Real>;
		const auto unitPhase = std::fmod(phase / pi2, Real(1));
		const auto length = unitPhase > tilt ? (1.0 - tilt) : tilt + std::numeric_limits<Real>::denorm_min();
		const auto distance = std::abs(unitPhase - tilt);
		const auto value = 1.0 - distance / length;
		return 2.0 * value - 1.0;
	}

	template <std::floating_point Real>
	inline Real Pwm(Real phase, Real fill) {
		constexpr auto pi2 = Real(2) * std::numbers::pi_v<Real>;
		const auto unitPhase = std::fmod(phase / pi2, Real(1));
		return Real(unitPhase < fill || fill >= 1.0);
	}

} // namespace impl

//------------------------------------------------------------------------------
// Constant tone
//------------------------------------------------------------------------------

/// <summary> Generate a sine wave. </summary>
/// <param name="output"> The generated wave is written here. </param>
/// <param name="sampleRate"> The sample rate of the generated wave. </param>
/// <param name="frequency"> The frequency of the wave. </param>
/// <param name="startPhase"> The starting phase of the wave. </param>
template <mutable_signal_or_view_r SignalOut, std::floating_point Real>
void SineWave(SignalOut&& output,
			  Real sampleRate,
			  Real frequency,
			  Real startPhase = Real(0)) {
	impl::GenericWave(output,
					  sampleRate,
					  frequency,
					  startPhase,
					  [](const auto& arg) { return std::sin(arg); });
}


/// <summary> Generate a sine wave. </summary>
///	<param name="length"> The length of the generated wave in samples. </param>
/// <param name="sampleRate"> The sample rate of the generated wave. </param>
/// <param name="frequency"> The frequency of the wave. </param>
/// <param name="startPhase"> The starting phase of the wave. </param>
template <class T, eSignalDomain Domain, std::floating_point Real>
BasicSignal<T, Domain> SineWave(size_t length,
								Real sampleRate,
								Real frequency,
								Real startPhase = Real(0)) {
	BasicSignal<T, Domain> signal(length);
	SineWave(signal, sampleRate, frequency, startPhase);
	return signal;
}


/// <summary> Generate a sawtooth wave. </summary>
/// <param name="output"> The generated wave is written here. </param>
/// <param name="sampleRate"> The sample rate of the generated wave. </param>
/// <param name="frequency"> The frequency of the wave. </param>
/// <param name="startPhase"> The starting phase of the wave. </param>
///	<param name="tilt"> 0 for a vertical rising edge, 1 for a vertical falling edge. </param>
template <mutable_signal_or_view_r SignalOut, std::floating_point Real>
void SawtoothWave(SignalOut&& output,
				  Real sampleRate,
				  Real frequency,
				  Real startPhase = Real(0),
				  Real tilt = Real(1.0)) {
	impl::GenericWave(output,
					  sampleRate,
					  frequency,
					  startPhase,
					  [tilt](const auto& arg) { return impl::Sawtooth(arg, tilt); });
}


/// <summary> Generate a sawtooth wave. </summary>
///	<param name="length"> The length of the generated wave in samples. </param>
/// <param name="sampleRate"> The sample rate of the generated wave. </param>
/// <param name="frequency"> The frequency of the wave. </param>
/// <param name="startPhase"> The starting phase of the wave. </param>
///	<param name="tilt"> 0 for a vertical rising edge, 1 for a vertical falling edge. </param>
template <class T, eSignalDomain Domain, std::floating_point Real>
BasicSignal<T, Domain> SawtoothWave(size_t length,
									Real sampleRate,
									Real frequency,
									Real startPhase = Real(0),
									Real tilt = Real(1)) {
	BasicSignal<T, Domain> signal(length);
	SawtoothWave(signal, sampleRate, frequency, startPhase, tilt);
	return signal;
}


/// <summary> Generate a PWM wave. </summary>
/// <param name="output"> The generated wave is written here. </param>
/// <param name="sampleRate"> The sample rate of the generated wave. </param>
/// <param name="frequency"> The frequency of the wave. </param>
/// <param name="startPhase"> The starting phase of the wave. </param>
///	<param name="dutyCycle"> The PWM duty cycle. </param>
template <mutable_signal_or_view_r SignalOut, std::floating_point Real>
void PwmWave(SignalOut&& output,
			 Real sampleRate,
			 Real frequency,
			 Real startPhase = Real(0),
			 Real dutyCycle = Real(0.5)) {
	impl::GenericWave(output,
					  sampleRate,
					  frequency,
					  startPhase,
					  [dutyCycle](const auto& arg) { return impl::Pwm(arg, dutyCycle); });
}


/// <summary> Generate a PWM wave. </summary>
/// <param name="length"> The length of the generated wave in samples. </param>
/// <param name="sampleRate"> The sample rate of the generated wave. </param>
/// <param name="frequency"> The frequency of the wave. </param>
/// <param name="startPhase"> The starting phase of the wave. </param>
///	<param name="dutyCycle"> The PWM duty cycle. </param>
template <class T, eSignalDomain Domain, std::floating_point Real>
BasicSignal<T, Domain> PwmWave(size_t length,
							   Real sampleRate,
							   Real frequency,
							   Real startPhase = Real(0),
							   Real dutyCycle = Real(0.5)) {
	BasicSignal<T, Domain> signal(length);
	PwmWave(signal, sampleRate, frequency, startPhase, dutyCycle);
	return signal;
}


/// <summary> Generate a square wave. </summary>
/// <param name="output"> The generated wave is written here. </param>
/// <param name="sampleRate"> The sample rate of the generated wave. </param>
/// <param name="frequency"> The frequency of the wave. </param>
/// <param name="startPhase"> The starting phase of the wave. </param>
template <mutable_signal_or_view_r SignalOut, std::floating_point Real>
void SquareWave(SignalOut&& output,
				Real sampleRate,
				Real frequency,
				Real startPhase = Real(0)) {
	using R = scalar_type_t<std::decay_t<SignalOut>>;
	impl::GenericWave(output,
					  sampleRate,
					  frequency,
					  startPhase,
					  [](const auto& arg) { return impl::Pwm(arg, Real(0.5)); });
	output *= R(2.0);
	output -= R(1.0);
}


/// <summary> Generate a square wave. </summary>
/// <param name="length"> The length of the generated wave in samples. </param>
/// <param name="sampleRate"> The sample rate of the generated wave. </param>
/// <param name="frequency"> The frequency of the wave. </param>
/// <param name="startPhase"> The starting phase of the wave. </param>
template <class T, eSignalDomain Domain, std::floating_point Real>
BasicSignal<T, Domain> SquareWave(size_t length,
								  Real sampleRate,
								  Real frequency,
								  Real startPhase = Real(0)) {
	BasicSignal<T, Domain> signal(length);
	SquareWave(signal, sampleRate, frequency, startPhase);
	return signal;
}


//------------------------------------------------------------------------------
// Chirp
//------------------------------------------------------------------------------


/// <summary> Generate a sine wave with a frequency sweep. </summary>
/// <param name="output"> The generated wave is written here. </param>
/// <param name="sampleRate"> The sample rate of the generated wave. </param>
/// <param name="startFrequency"> The frequency at the start of the generated wave. </param>
/// <param name="endFrequency"> The frequency at the end of the generated wave. </param>
/// <param name="startPhase"> The starting phase of the wave. </param>
/// <remarks> The frequency is linearly interpolated between the start and end frequencies. </remarks>
template <mutable_signal_or_view_r SignalOut, std::floating_point Real>
void SineChirp(SignalOut&& output,
			   Real sampleRate,
			   Real startFrequency,
			   Real endFrequency,
			   Real startPhase = Real(0)) {
	impl::GenericChirp(output,
					   sampleRate,
					   startFrequency,
					   endFrequency,
					   startPhase,
					   [](const auto& arg) { return std::sin(arg); });
}


/// <summary> Generate a sine wave with a frequency sweep. </summary>
/// <param name="length"> The length of the generated wave in samples. </param>
/// <param name="sampleRate"> The sample rate of the generated wave. </param>
/// <param name="startFrequency"> The frequency at the start of the generated wave. </param>
/// <param name="endFrequency"> The frequency at the end of the generated wave. </param>
/// <param name="startPhse"> The starting phase of the wave. </param>
/// <remarks> The frequency is linearly interpolated between the start and end frequencies. </remarks>
template <class T, eSignalDomain Domain, std::floating_point Real>
BasicSignal<T, Domain> SineChirp(size_t length,
								 Real sampleRate,
								 Real startFrequency,
								 Real endFrequency,
								 Real startPhse = Real(0)) {
	BasicSignal<T, Domain> signal(length);
	SineChirp(signal, sampleRate, startFrequency, endFrequency, startPhse);
	return signal;
}


/// <summary> Generate a sawtooth wave with a frequency sweep. </summary>
/// <param name="output"> The generated wave is written here. </param>
/// <param name="sampleRate"> The sample rate of the generated wave. </param>
/// <param name="startFrequency"> The frequency at the start of the generated wave. </param>
/// <param name="endFrequency"> The frequency at the end of the generated wave. </param>
/// <param name="phase"> The starting phase of the wave. </param>
///	<param name="tilt"> 0 for a vertical rising edge, 1 for a vertical falling edge. </param>
/// <remarks> The frequency is linearly interpolated between the start and end frequencies. </remarks>
template <mutable_signal_or_view_r SignalOut, std::floating_point Real>
void SawtoothChirp(SignalOut&& output,
				   Real sampleRate,
				   Real startFrequency,
				   Real endFrequency,
				   Real phase = Real(0),
				   Real tilt = Real(1.0)) {
	impl::GenericChirp(output,
					   sampleRate,
					   startFrequency,
					   endFrequency,
					   phase,
					   [tilt](const auto& arg) { return impl::Sawtooth(arg, tilt); });
}


/// <summary> Generate a sawtooth wave with a frequency sweep. </summary>
/// <param name="length"> The length of the generated wave in samples. </param>
/// <param name="sampleRate"> The sample rate of the generated wave. </param>
/// <param name="startFrequency"> The frequency at the start of the generated wave. </param>
/// <param name="endFrequency"> The frequency at the end of the generated wave. </param>
/// <param name="startPhase"> The starting phase of the wave. </param>
///	<param name="tilt"> 0 for a vertical rising edge, 1 for a vertical falling edge. </param>
/// <remarks> The frequency is linearly interpolated between the start and end frequencies. </remarks>
template <class T, eSignalDomain Domain, std::floating_point Real>
BasicSignal<T, Domain> SawtoothChirp(size_t length,
									 Real sampleRate,
									 Real startFrequency,
									 Real endFrequency,
									 Real startPhase = 0,
									 Real tilt = 1.0) {
	BasicSignal<T, Domain> signal(length);
	SawtoothChirp(signal, sampleRate, startFrequency, endFrequency, startPhase, tilt);
	return signal;
}


/// <summary> Generate a PWM wave with a frequency sweep. </summary>
/// <param name="output"> The generated wave is written here. </param>
/// <param name="sampleRate"> The sample rate of the generated wave. </param>
/// <param name="startFrequency"> The frequency at the start of the generated wave. </param>
/// <param name="endFrequency"> The frequency at the end of the generated wave. </param>
/// <param name="startPhase"> The starting phase of the wave. </param>
///	<param name="dutyCycle"> The PWM duty cycle. </param>
/// <remarks> The frequency is linearly interpolated between the start and end frequencies. </remarks>
template <mutable_signal_or_view_r SignalOut, std::floating_point Real>
void PwmChirp(SignalOut&& output,
			  Real sampleRate,
			  Real startFrequency,
			  Real endFrequency,
			  Real startPhase = Real(0),
			  Real dutyCycle = Real(0.5)) {
	impl::GenericChirp(output,
					   sampleRate,
					   startFrequency,
					   endFrequency,
					   startPhase,
					   [dutyCycle](const auto& arg) { return impl::Pwm(arg, dutyCycle); });
}


/// <summary> Generate a PWM wave with a frequency sweep. </summary>
/// <param name="length"> The length of the generated wave in samples. </param>
/// <param name="sampleRate"> The sample rate of the generated wave. </param>
/// <param name="startFrequency"> The frequency at the start of the generated wave. </param>
/// <param name="endFrequency"> The frequency at the end of the generated wave. </param>
/// <param name="startPhase"> The starting phase of the wave. </param>
///	<param name="dutyCycle"> The PWM duty cycle. </param>
/// <remarks> The frequency is linearly interpolated between the start and end frequencies. </remarks>
template <class T, eSignalDomain Domain, std::floating_point Real>
BasicSignal<T, Domain> PwmChirp(size_t length,
								Real sampleRate,
								Real startFrequency,
								Real endFrequency,
								Real startPhase = Real(0),
								Real dutyCycle = Real(0.5)) {
	BasicSignal<T, Domain> signal(length);
	PwmChirp(signal, sampleRate, startFrequency, endFrequency, startPhase, dutyCycle);
	return signal;
}


/// <summary> Generate a square wave with a frequency sweep. </summary>
/// <param name="output"> The generated wave is written here. </param>
/// <param name="sampleRate"> The sample rate of the generated wave. </param>
/// <param name="startFrequency"> The frequency at the start of the generated wave. </param>
/// <param name="endFrequency"> The frequency at the end of the generated wave. </param>
/// <param name="startPhase"> The starting phase of the wave. </param>
/// <remarks> The frequency is linearly interpolated between the start and end frequencies. </remarks>
template <mutable_signal_or_view_r SignalOut, std::floating_point Real>
void SquareChirp(SignalOut&& output,
				 Real sampleRate,
				 Real startFrequency,
				 Real endFrequency,
				 Real startPhase = Real(0)) {
	using R = scalar_type_t<std::decay_t<SignalOut>>;
	impl::GenericChirp(output,
					   sampleRate,
					   startFrequency,
					   endFrequency,
					   startPhase,
					   [](const auto& arg) { return impl::Pwm(arg, Real(0.5)); });
	output *= R(2.0);
	output -= R(1.0);
}


/// <summary> Generate a square wave with a frequency sweep. </summary>
/// <param name="length"> The length of the generated wave in samples. </param>
/// <param name="sampleRate"> The sample rate of the generated wave. </param>
/// <param name="startFrequency"> The frequency at the start of the generated wave. </param>
/// <param name="endFrequency"> The frequency at the end of the generated wave. </param>
/// <param name="startPhase"> The starting phase of the wave. </param>
/// <remarks> The frequency is linearly interpolated between the start and end frequencies. </remarks>
template <class T, eSignalDomain Domain, std::floating_point Real>
BasicSignal<T, Domain> SquareChirp(size_t length,
								   Real sampleRate,
								   Real startFrequency,
								   Real endFrequency,
								   Real startPhase = 0) {
	BasicSignal<T, Domain> signal(length);
	SquareChirp(signal, sampleRate, startFrequency, endFrequency, startPhase);
	return signal;
}

} // namespace dspbb