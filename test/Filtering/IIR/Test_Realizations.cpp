#include <catch2/catch.hpp>
#include <dspbb/Filtering/IIR/Realizations.hpp>
#include <dspbb/LTISystems/Systems.hpp>
#include <dspbb/Math/FFT.hpp>
#include <dspbb/Math/Functions.hpp>
#include <dspbb/Math/Statistics.hpp>


using namespace dspbb;
using namespace std::complex_literals;
using real_t = double;


const DiscreteZeroPoleGain<real_t> sys = {
	3.2f,
	{ 0.35f, 0.77f + 0.2if, 0.77f - 0.2if, 0.4f + 0.6if, 0.4f - 0.6if },
	{ -0.2f, -0.6f, -0.7f + 0.2if, -0.7f - 0.2if, -0.35f + 0.6if, -0.35f - 0.6if }
};
const TransferFunction tf{ sys };
const CascadedBiquad cascade{ sys };

Signal<real_t, FREQUENCY_DOMAIN> re;
Signal<real_t, FREQUENCY_DOMAIN> im;
const Signal<real_t, TIME_DOMAIN> input = { 0.5f, 0.9f, 1.4f, -1.3f, -0.6f, -0.3f };
const Signal<real_t, TIME_DOMAIN> response = []() {
	TimeSignal<real_t> num{ tf.numerator.Coefficients().begin(), tf.numerator.Coefficients().end() };
	TimeSignal<real_t> den{ tf.denominator.Coefficients().begin(), tf.denominator.Coefficients().end() };
	TimeSignal<real_t> padded = input;
	std::reverse(num.begin(), num.end());
	std::reverse(den.begin(), den.end());
	num.Resize(1000);
	den.Resize(1000);
	padded.Resize(1000);
	const auto numF = FourierTransform(num, false);
	const auto denF = FourierTransform(den, false);
	const auto inputF = FourierTransform(padded, false);
	return InverseFourierTransformR(inputF * numF / denF, 1000);
}();


TEST_CASE("Direct form I", "[IIR realizations]") {
	TimeSignal<real_t> out;

	DirectFormI<real_t> state{ std::max(sys.zeros.NumRoots(), sys.poles.NumRoots()) };
	for (size_t i = 0; i < 1000; ++i) {
		const real_t u = i < input.Size() ? input[i] : 0.0f;
		out.PushBack(state.Feed(u, tf));
	}

	const real_t similarity = DotProduct(response, out) / Norm(out) / Norm(response);
	REQUIRE(similarity == Approx(1));
}

TEST_CASE("Direct form II", "[IIR realizations]") {
	TimeSignal<real_t> out;

	DirectFormII<real_t> state{ std::max(sys.zeros.NumRoots(), sys.poles.NumRoots()) };
	for (size_t i = 0; i < 1000; ++i) {
		const real_t u = i < input.Size() ? input[i] : 0.0f;
		out.PushBack(state.Feed(u, tf));
	}

	const real_t similarity = DotProduct(response, out) / Norm(out) / Norm(response);
	REQUIRE(similarity == Approx(1));
}

TEST_CASE("Cascaded biquad form", "[IIR realizations]") {
	TimeSignal<real_t> out;

	CascadedForm<real_t> state{ std::max(sys.zeros.NumRoots(), sys.poles.NumRoots()) };
	for (size_t i = 0; i < 1000; ++i) {
		const real_t u = i < input.Size() ? input[i] : 0.0f;
		out.PushBack(state.Feed(u, cascade));
	}

	const real_t similarity = DotProduct(response, out) / Norm(out) / Norm(response);
	REQUIRE(similarity == Approx(1));
}