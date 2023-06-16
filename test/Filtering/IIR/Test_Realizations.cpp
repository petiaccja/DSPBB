#include <dspbb/Filtering/IIR/Realizations.hpp>
#include <dspbb/LTISystems/Systems.hpp>
#include <dspbb/Math/FFT.hpp>
#include <dspbb/Math/Functions.hpp>
#include <dspbb/Math/Statistics.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>


using namespace dspbb;
using namespace std::complex_literals;
using real_t = double;
using Catch::Approx;


const DiscreteZeroPoleGain<real_t> sys = {
	3.2f,
	{ 0.35f, 0.77f + 0.2if, 0.77f - 0.2if, 0.4f + 0.6if, 0.4f - 0.6if },
	{ -0.2f, -0.6f, -0.7f + 0.2if, -0.7f - 0.2if, -0.35f + 0.6if, -0.35f - 0.6if }
};
const TransferFunction tf{ sys };
const CascadedBiquad cascade{ sys };

const BasicSignal<real_t, TIME_DOMAIN> input = { 0.5f, 0.9f, 1.4f, -1.3f, -0.6f, -0.3f };
const BasicSignal<real_t, TIME_DOMAIN> response = []() {
	Signal<real_t> num{ tf.numerator.coefficients().begin(), tf.numerator.coefficients().end() };
	Signal<real_t> den{ tf.denominator.coefficients().begin(), tf.denominator.coefficients().end() };
	Signal<real_t> padded = input;
	std::reverse(num.begin(), num.end());
	std::reverse(den.begin(), den.end());
	num.resize(1000);
	den.resize(1000);
	padded.resize(1000);
	const auto numF = Fft(num, FFT_HALF);
	const auto denF = Fft(den, FFT_HALF);
	const auto inputF = Fft(padded, FFT_HALF);
	return Ifft(inputF * numF / denF, FFT_HALF, true);
}();

//------------------------------------------------------------------------------
// feed
//------------------------------------------------------------------------------

TEST_CASE("Direct form I feed", "[IIR realizations]") {
	Signal<real_t> out;

	DirectFormI<real_t> state{ std::max(sys.zeros.num_roots(), sys.poles.num_roots()) };
	for (size_t i = 0; i < 1000; ++i) {
		const real_t u = i < input.size() ? input[i] : 0.0f;
		out.push_back(state.feed(u, tf));
	}

	const real_t similarity = DotProduct(response, out) / Norm(out) / Norm(response);
	REQUIRE(similarity == Approx(1));
}

TEST_CASE("Direct form II feed", "[IIR realizations]") {
	Signal<real_t> out;

	DirectFormII<real_t> state{ std::max(sys.zeros.num_roots(), sys.poles.num_roots()) };
	for (size_t i = 0; i < 1000; ++i) {
		const real_t u = i < input.size() ? input[i] : 0.0f;
		out.push_back(state.feed(u, tf));
	}

	const real_t similarity = DotProduct(response, out) / Norm(out) / Norm(response);
	REQUIRE(similarity == Approx(1));
}

TEST_CASE("Cascaded biquad form feed", "[IIR realizations]") {
	Signal<real_t> out;

	CascadedForm<real_t> state{ std::max(sys.zeros.num_roots(), sys.poles.num_roots()) };
	for (size_t i = 0; i < 1000; ++i) {
		const real_t u = i < input.size() ? input[i] : 0.0f;
		out.push_back(state.feed(u, cascade));
	}

	const real_t similarity = DotProduct(response, out) / Norm(out) / Norm(response);
	REQUIRE(similarity == Approx(1));
}

//------------------------------------------------------------------------------
// feed different input type
//------------------------------------------------------------------------------

const DiscreteZeroPoleGain<double> sysd = {
	0.0,
	{ -0.6i, 0.6i },
	{ -0.55i, 0.55i }
};
const auto tfd = TransferFunction{ sysd };
const auto cascaded = CascadedBiquad{ sysd };

const DiscreteZeroPoleGain<float> sysf = {
	0.0f,
	{ -0.6if, 0.6if },
	{ -0.55if, 0.55if }
};
const auto tff = TransferFunction{ sysf };
const auto cascadef = CascadedBiquad{ sysf };

TEST_CASE("Direct form I feed float/double", "[IIR realizations]") {
	constexpr float inputf = 1.0f;
	DirectFormI<float> df1{ sys.order() };
	DirectFormII<float> df2{ sys.order() };
	CascadedForm<float> cf{ sys.order() };

	[[maybe_unused]] const auto out1 = df1.feed(inputf, tfd);
	[[maybe_unused]] const auto out2 = df2.feed(inputf, tfd);
	[[maybe_unused]] const auto out3 = cf.feed(inputf, cascaded);
	REQUIRE(std::is_same_v<float, std::decay_t<decltype(out1)>>);
	REQUIRE(std::is_same_v<float, std::decay_t<decltype(out2)>>);
	REQUIRE(std::is_same_v<float, std::decay_t<decltype(out3)>>);
}

TEST_CASE("Direct form I feed complex<float>/float", "[IIR realizations]") {
	constexpr std::complex<float> inputcf = 1.0f;
	DirectFormI<std::complex<float>> df1{ sys.order() };
	DirectFormII<std::complex<float>> df2{ sys.order() };
	CascadedForm<std::complex<float>> cf{ sys.order() };

	[[maybe_unused]] const auto out1 = df1.feed(inputcf, tff);
	[[maybe_unused]] const auto out2 = df2.feed(inputcf, tff);
	[[maybe_unused]] const auto out3 = cf.feed(inputcf, cascadef);
	REQUIRE(std::is_same_v<std::complex<float>, std::decay_t<decltype(out1)>>);
	REQUIRE(std::is_same_v<std::complex<float>, std::decay_t<decltype(out2)>>);
	REQUIRE(std::is_same_v<std::complex<float>, std::decay_t<decltype(out3)>>);
}

//------------------------------------------------------------------------------
// Direct form I
//------------------------------------------------------------------------------

TEST_CASE("Direct form I default construct", "[IIR realizations]") {
	DirectFormI<float> state;
	REQUIRE(state.order() == 0);
}

TEST_CASE("Direct form I construct", "[IIR realizations]") {
	DirectFormI<float> state{ 12 };
	REQUIRE(state.order() == 12);
}

TEST_CASE("Direct form I order", "[IIR realizations]") {
	DirectFormI<float> state;
	state.order(12);
	REQUIRE(state.order() == 12);
}

TEST_CASE("Direct form I reset", "[IIR realizations]") {
	DirectFormI<float> state{ 2 };
	DiscreteTransferFunction<float> tf2{ Polynomial<float>{ 1, 1, 1 }, Polynomial<float>{ 1, 1, 1 } };
	for (int i = 0; i < 10; ++i) {
		REQUIRE(0.0f != state.feed(1.0f, tf2));
	}
	state.reset();
	for (int i = 0; i < 10; ++i) {
		REQUIRE(0.0f == state.feed(0.0f, tf2));
	}
}

//------------------------------------------------------------------------------
// Direct form II
//------------------------------------------------------------------------------

TEST_CASE("Direct form II default construct", "[IIR realizations]") {
	DirectFormII<float> state;
	REQUIRE(state.order() == 0);
}

TEST_CASE("Direct form II construct", "[IIR realizations]") {
	DirectFormII<float> state{ 12 };
	REQUIRE(state.order() == 12);
}

TEST_CASE("Direct form II order", "[IIR realizations]") {
	DirectFormII<float> state;
	state.order(12);
	REQUIRE(state.order() == 12);
}

TEST_CASE("Direct form II reset", "[IIR realizations]") {
	DirectFormII<float> state{ 2 };
	DiscreteTransferFunction<float> tf2{ Polynomial<float>{ 1, 1, 1 }, Polynomial<float>{ 1, 1, 1 } };
	for (int i = 0; i < 10; ++i) {
		REQUIRE(0.0f != state.feed(1.0f, tf2));
	}
	state.reset();
	for (int i = 0; i < 10; ++i) {
		REQUIRE(0.0f == state.feed(0.0f, tf2));
	}
}

//------------------------------------------------------------------------------
// Cascaded form
//------------------------------------------------------------------------------

TEST_CASE("Cascaded form default construct", "[IIR realizations]") {
	CascadedForm<float> state;
	REQUIRE(state.order() == 0);
}

TEST_CASE("Cascaded form construct", "[IIR realizations]") {
	CascadedForm<float> state{ 12 };
	REQUIRE(state.order() == 12);
}

TEST_CASE("Cascaded form construct odd", "[IIR realizations]") {
	CascadedForm<float> state{ 11 };
	REQUIRE(state.order() == 12); // Can't have even orders
}

TEST_CASE("Cascaded form order", "[IIR realizations]") {
	CascadedForm<float> state;
	state.order(12);
	REQUIRE(state.order() == 12);
}

TEST_CASE("Cascaded form order odd", "[IIR realizations]") {
	CascadedForm<float> state;
	state.order(11);
	REQUIRE(state.order() == 12); // Can't have even orders
}

TEST_CASE("Cascaded form reset", "[IIR realizations]") {
	CascadedForm<float> state{ 2 };
	const CascadedBiquad s{ DiscreteZeroPoleGain<float>{ 1.0f, { 1.0f, 2.0f }, { -1.0f, -2.0f } } };
	for (int i = 0; i < 10; ++i) {
		REQUIRE(0.0f != state.feed(1.0f, s));
	}
	state.reset();
	for (int i = 0; i < 10; ++i) {
		REQUIRE(0.0f == state.feed(0.0f, s));
	}
}