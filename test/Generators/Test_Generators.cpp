#include <dspbb/Generators/Spaces.hpp>
#include <dspbb/Generators/Waveforms.hpp>
#include <dspbb/Math/Statistics.hpp>
#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/Primitives/SignalView.hpp>

#include <catch2/catch.hpp>


using namespace dspbb;


constexpr int sampleRate = 44100;
constexpr float frequency = 89.f;
const float cycle = sampleRate / frequency;


TEST_CASE("Linspace inclusive", "[Generators]") {
	const auto s = LinSpace<float, TIME_DOMAIN>(6.28f, 2.718f, 23, true);
	REQUIRE(s.Size() == 23);
	REQUIRE(s[0] == Approx(6.28f));
	REQUIRE(s[22] == Approx(2.718f));
	REQUIRE(s[17] - s[12] == Approx(5.0f * (s[9] - s[8])));
}

TEST_CASE("Linspace exclusive", "[Generators]") {
	const auto s = LinSpace<float, TIME_DOMAIN>(6.28f, 2.718f, 23, false);
	REQUIRE(s.Size() == 23);
	REQUIRE(s[0] == Approx(6.28f));
	REQUIRE(s[22] + (s[22] - s[21]) == Approx(2.718f));
	REQUIRE(s[17] - s[12] == Approx(5.0f * (s[9] - s[8])));
}

TEST_CASE("Logspace", "[Generators]") {
	const auto s = LogSpace<float, TIME_DOMAIN>(2.0f, 4.0f, 23, 10.f, true);
	REQUIRE(s.Size() == 23);
	REQUIRE(s[0] == Approx(100.f));
	REQUIRE(s[22] == Approx(10000.f));
	const auto quot = SignalView<const float>(s.begin(), s.end() - 1) / SignalView<const float>(s.begin() + 1, s.end());
	REQUIRE(Max(quot) == Approx(Min(quot)));
}


TEST_CASE("Sine wave", "[Generators]") {
	const auto s = SineWave<float, TIME_DOMAIN>(4410, 44100, frequency, 0.5f);
	REQUIRE(s[0] == Approx(std::sin(0.5f)));
	REQUIRE(s[size_t(cycle)] == Approx(s[0]).margin(0.02f));
	REQUIRE(s[size_t(cycle * (2.f * pi_v<float> - 0.5f) / 2.f / pi_v<float>)] == Approx(0.0f).margin(0.02f));
	REQUIRE(s[size_t(cycle * (2.5f * pi_v<float> - 0.5f) / 2.f / pi_v<float>)] == Approx(1.0f).margin(0.02f));
}


TEST_CASE("Sawtooth wave fw", "[Generators]") {
	const auto s = SawtoothWave<float, TIME_DOMAIN>(4410, 44100, frequency, 0.0f, 1.0f);
	REQUIRE(s[0] == Approx(-1.0f));
	REQUIRE(s[size_t(cycle * 0.5f)] == Approx(0.0f).margin(0.02f));
	REQUIRE(s[size_t(cycle)] == Approx(1.0f).margin(0.02f));
}

TEST_CASE("Sawtooth wave bw", "[Generators]") {
	const auto s = SawtoothWave<float, TIME_DOMAIN>(4410, 44100, frequency, 0.0f, 0.0f);
	REQUIRE(s[0] == Approx(1.0f));
	REQUIRE(s[size_t(cycle * 0.5f)] == Approx(0.0f).margin(0.02f));
	REQUIRE(s[size_t(cycle)] == Approx(-1.0f).margin(0.02f));
}

TEST_CASE("Sawtooth wave triangle", "[Generators]") {
	const auto s = SawtoothWave<float, TIME_DOMAIN>(4410, 44100, frequency, 0.0f, 0.6f);
	REQUIRE(s[0] == Approx(-1.0f));
	REQUIRE(s[size_t(cycle * 0.6f)] == Approx(1.0f).margin(0.02f));
	REQUIRE(s[size_t(cycle) + 1] == Approx(-1.0f).margin(0.02f));
}

TEST_CASE("PWM wave empty", "[Generators]") {
	const auto s = PwmWave<float, TIME_DOMAIN>(4410, 44100, frequency, 0.0f, 0.0f);
	REQUIRE(Max(s) == Approx(0));
	REQUIRE(Min(s) == Approx(0));
}

TEST_CASE("PWM wave full", "[Generators]") {
	const auto s = PwmWave<float, TIME_DOMAIN>(4410, 44100, frequency, 0.0f, 1.0f);
	REQUIRE(Max(s) == Approx(1));
	REQUIRE(Min(s) == Approx(1));
}

TEST_CASE("Sawtooth wave frac", "[Generators]") {
	const auto s = PwmWave<float, TIME_DOMAIN>(4410, 44100, frequency, 0.0f, 0.6f);
	REQUIRE(Max(s) == Approx(1));
	REQUIRE(Min(s) == Approx(0));
	REQUIRE(s[0] == Approx(1.0f));
	REQUIRE(s[size_t(cycle * 0.55f)] == Approx(1.0f));
	REQUIRE(s[size_t(cycle * 0.65f)] == Approx(0.0f));
	REQUIRE(s[size_t(cycle * 0.99f)] == Approx(0.0f));
	REQUIRE(s[size_t(cycle * 1.01f)] == Approx(1.0f));
}

TEST_CASE("Square wave", "[Generators]") {
	const auto s = SquareWave<float, TIME_DOMAIN>(4410, 44100, frequency, 0.0f);
	REQUIRE(Max(s) == Approx(1));
	REQUIRE(Min(s) == Approx(-1));

	REQUIRE(s[0] == Approx(1.0f));
	REQUIRE(s[size_t(cycle * 0.45f)] == Approx(1.0f));
	REQUIRE(s[size_t(cycle * 0.55f)] == Approx(-1.0f));
	REQUIRE(s[size_t(cycle * 0.99f)] == Approx(-1.0f));
	REQUIRE(s[size_t(cycle * 1.01f)] == Approx(1.0f));
}

// Enough to test the base chirp phase function.
TEST_CASE("Chirp phase", "[Generators]") {
	Signal<float> s(512);
	const double phase = 1.55;
	const double startFrequency = 1150.;
	const double endFrequency = 2320.;
	impl::GenericChirp(s, sampleRate, startFrequency, endFrequency, phase, [](const auto& passThrough) { return passThrough; });
	REQUIRE(s[0] == Approx(phase));
	REQUIRE(*(s.begin() + 1) - *(s.begin()) == Approx(2 * pi_v<double> * startFrequency / sampleRate).epsilon(0.01f));
	REQUIRE(*(s.end() - 1) - *(s.end() - 2) == Approx(2 * pi_v<double> * endFrequency / sampleRate).epsilon(0.01f));
	REQUIRE(Max(SignalView<float>(s.begin(), s.end() - 1) - SignalView<float>(s.begin() + 1, s.end())) < 0.0f);
}

TEST_CASE("Square chirp", "[Generators]") {
	const auto s = SquareChirp<float, TIME_DOMAIN>(4410, 44100, 2 * frequency, frequency, 0.0f);
	REQUIRE(Max(s) == Approx(1));
	REQUIRE(Min(s) == Approx(-1));
	REQUIRE(s[0] == Approx(1.0f));
}