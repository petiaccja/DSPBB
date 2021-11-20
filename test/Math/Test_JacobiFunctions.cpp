#include "../TestUtils.hpp"

#include <dspbb/Generators/Spaces.hpp>
#include <dspbb/Math/Functions.hpp>
#include <dspbb/Math/JacobiFunctions.hpp>
#include <dspbb/Math/Statistics.hpp>
#include <dspbb/Primitives/SignalView.hpp>

#include <catch2/catch.hpp>

using namespace dspbb;
using namespace std::complex_literals;


//------------------------------------------------------------------------------
// Constants
//------------------------------------------------------------------------------

const std::array variants = { 1, 2, 3, 4 };


//------------------------------------------------------------------------------
// Lattice transform regressions
//------------------------------------------------------------------------------

TEST_CASE("Shift scalar", "[Jacobi functions]") {
	std::array<float, 5> inputs = { -1234567891011.0f, -23.15f, 0.86f, 3.14f, 987654321.0f };
	std::array<float, 5> remainders = { 0.0f, -0.15f, -0.14f, 0.14f, 0.0f };
	std::array<int, 5> counts = { 0, 23, -1, -3, 0 };
	for (int i = 0; i < 5; ++i) {
		const auto [remainder, count] = ShiftScalar(inputs[i]);
		REQUIRE(remainder == Approx(remainders[i]).epsilon(1e-5f));
		REQUIRE(count % 8 == counts[i] % 8);
	}
}

TEST_CASE("Shift variant", "[Jacobi functions]") {
	REQUIRE(ShiftVariant(1, 5) == 1);
	REQUIRE(ShiftVariant(1, 6) == 1);
	REQUIRE(ShiftVariant(2, 5) == 2);
	REQUIRE(ShiftVariant(2, 6) == 2);
	REQUIRE(ShiftVariant(3, 5) == 4);
	REQUIRE(ShiftVariant(3, 6) == 3);
	REQUIRE(ShiftVariant(4, 5) == 3);
	REQUIRE(ShiftVariant(4, 6) == 4);
}

TEST_CASE("Shift multiplier", "[Jacobi functions]") {
	REQUIRE(ShiftMultiplier<float>(1, 5) == ApproxComplex(std::polar(1.0f, -5.0f / 4.0f * pi_v<float>)));
	REQUIRE(ShiftMultiplier<float>(1, 14) == ApproxComplex(std::polar(1.0f, -6.0f / 4.0f * pi_v<float>)));
	REQUIRE(ShiftMultiplier<float>(2, 5) == ApproxComplex(std::polar(1.0f, -5.0f / 4.0f * pi_v<float>)));
	REQUIRE(ShiftMultiplier<float>(2, 14) == ApproxComplex(std::polar(1.0f, -6.0f / 4.0f * pi_v<float>)));
	REQUIRE(ShiftMultiplier<float>(3, 5) == ApproxComplex(std::polar(1.0f, 0.0f)));
	REQUIRE(ShiftMultiplier<float>(3, 6) == ApproxComplex(std::polar(1.0f, 0.0f)));
	REQUIRE(ShiftMultiplier<float>(4, 5) == ApproxComplex(std::polar(1.0f, 0.0f)));
	REQUIRE(ShiftMultiplier<float>(4, 6) == ApproxComplex(std::polar(1.0f, 0.0f)));
}

TEST_CASE("Shift tau negative direction var 1", "[Jacobi functions]") {
	const std::complex<float> tau = { 12.99345f, 0.1f };
	const int variant = 1;
	const auto [newVariant, _, newTau, multiplier, __] = ShiftTau(variant, {}, tau);
	REQUIRE(real(newTau) == Approx(12.99345f - 13.0f).epsilon(1e-6f));
	REQUIRE(imag(newTau) == Approx(0.1f).epsilon(1e-6f));
	REQUIRE(multiplier == ApproxComplex(std::exp(1if * pi_v<float> * 13.0f / 4.0f)).epsilon(1e-6f));
	REQUIRE(newVariant == variant);
}

TEST_CASE("Shift tau positive direction var 2", "[Jacobi functions]") {
	const std::complex<float> tau = { -12.99345f, 0.1f };
	const int variant = 2;
	const auto [newVariant, _, newTau, multiplier, __] = ShiftTau(variant, {}, tau);
	REQUIRE(real(newTau) == Approx(-12.99345f + 13.0f).epsilon(1e-6f));
	REQUIRE(imag(newTau) == Approx(0.1f).epsilon(1e-6f));
	REQUIRE(multiplier == ApproxComplex(std::exp(1if * pi_v<float> * -13.0f / 4.0f)).epsilon(1e-6f));
	REQUIRE(newVariant == variant);
}

TEST_CASE("Shift tau negative direction var 3", "[Jacobi functions]") {
	const std::complex<float> tau = { 12.99345f, 0.1f };
	const int variant = 3;
	const auto [newVariant, _, newTau, multiplier, __] = ShiftTau(variant, {}, tau);
	REQUIRE(real(newTau) == Approx(12.99345f - 13.0f).epsilon(1e-6f));
	REQUIRE(imag(newTau) == Approx(0.1f).epsilon(1e-6f));
	REQUIRE(multiplier == ApproxComplex(1.0f).epsilon(1e-6f));
	REQUIRE(newVariant == 4);
}

TEST_CASE("Shift tau positive direction var 4", "[Jacobi functions]") {
	const std::complex<float> tau = { -12.99345f, 0.1f };
	const int variant = 4;
	const auto [newVariant, _, newTau, multiplier, __] = ShiftTau(variant, {}, tau);
	REQUIRE(real(newTau) == Approx(-12.99345f + 13.0f).epsilon(1e-6f));
	REQUIRE(imag(newTau) == Approx(0.1f).epsilon(1e-6f));
	REQUIRE(multiplier == ApproxComplex(1.0f).epsilon(1e-6f));
	REQUIRE(newVariant == 3);
}

TEST_CASE("Invert variant", "[Jacobi functions]") {
	REQUIRE(InvertVariant(1) == 1);
	REQUIRE(InvertVariant(2) == 4);
	REQUIRE(InvertVariant(3) == 3);
	REQUIRE(InvertVariant(4) == 2);
}

TEST_CASE("Invert multiplier #1", "[Jacobi functions]") {
	const auto z = 2.0f + 0.0if;
	const auto tau = 0.5if;
	const auto [factor, exponent] = InvertMultiplier(1, z, tau);
	REQUIRE(factor == ApproxComplex(-i_v<float> * std::sqrt(2.0f)));
	REQUIRE(exponent == ApproxComplex(-8.0f / pi_v<float>));
}

TEST_CASE("Invert multiplier #2,3,4", "[Jacobi functions]") {
	const auto z = 2.0f + 0.0if;
	const auto tau = 0.5if;
	for (int variant = 2; variant <= 4; ++variant) {
		const auto [factor, exponent] = InvertMultiplier(variant, z, tau);
		REQUIRE(factor == ApproxComplex(std::sqrt(2.0f)));
		REQUIRE(exponent == ApproxComplex(-8.0f / pi_v<float>));
	}
}

TEST_CASE("Invert tau #1", "[Jacobi functions]") {
	const auto z = 2.0f + 0.0if;
	const auto tau = 0.5if;
	const auto [newVariant, newZ, newTau, factor, exponent] = InvertTau(1, z, tau);
	REQUIRE(newZ == ApproxComplex(4if));
	REQUIRE(newTau == ApproxComplex(2if));
	REQUIRE(factor == ApproxComplex(-i_v<float> * std::sqrt(2.0f)));
	REQUIRE(exponent == ApproxComplex(-8.0f / pi_v<float>));
	REQUIRE(newVariant == 1);
}

TEST_CASE("Invert tau #2,3,4", "[Jacobi functions]") {
	const auto z = 2.0f + 0.0if;
	const auto tau = 0.5if;
	for (int variant = 2; variant <= 4; ++variant) {
		const auto [newVariant, newZ, newTau, factor, exponent] = InvertTau(variant, z, tau);
		REQUIRE(newZ == ApproxComplex(4if));
		REQUIRE(newTau == ApproxComplex(2if));
		REQUIRE(factor == ApproxComplex(std::sqrt(2.0f)));
		REQUIRE(exponent == ApproxComplex(-8.0f / pi_v<float>));
		REQUIRE(newVariant == 2 + 4 - variant);
	}
}

//------------------------------------------------------------------------------
// Lattice transform application identities
//------------------------------------------------------------------------------


template <class T>
std::pair<std::complex<T>, std::complex<T>> LatticeApplicationIdentity(const LatticeTransform<T>& control, const LatticeTransform<T>& trial) {
	const auto rcontrol = control.multiplier * ThetaSeries(control.variant, control.z, control.tau, control.exponent, 25);
	const auto rtrial = trial.multiplier * ThetaSeries(trial.variant, trial.z, trial.tau, trial.exponent, 25);
	return { rcontrol, rtrial };
}

const std::array latticeApplicationTaus = {
	0.73 + 1.49i,
	-3.17 + 0.49i,
	0.77i,
	-0.11 + 1.03i,
	-14.11 + 1.03i,
	7.83 + 1.03i,
};

const auto latticeApplicationZ = 0.7 + 0.3i;


TEST_CASE("Lattice shift application identity", "[Jacobi functions]") {
	for (int variant : variants) {
		for (auto& tau : latticeApplicationTaus) {
			UNSCOPED_INFO("variant=" << variant << "  tau=" << tau);

			const auto control = LatticeTransform<double>{ variant, latticeApplicationZ, tau };
			const auto trial = ShiftTau(variant, latticeApplicationZ, tau);
			const auto [rcontrol, rtrial] = LatticeApplicationIdentity(control, trial);

			REQUIRE(rcontrol == ApproxComplex(rtrial).epsilon(1e-3f));
		}
	}
}

TEST_CASE("Lattice inversion application identity", "[Jacobi functions]") {
	for (int variant : variants) {
		for (auto& tau : latticeApplicationTaus) {
			UNSCOPED_INFO("variant=" << variant << "  tau=" << tau);

			const auto control = LatticeTransform<double>{ variant, latticeApplicationZ, tau };
			const auto trial = InvertTau(variant, latticeApplicationZ, tau);
			const auto [rcontrol, rtrial] = LatticeApplicationIdentity(control, trial);

			REQUIRE(rcontrol == ApproxComplex(rtrial).epsilon(1e-3f));
		}
	}
}

TEST_CASE("Lattice rotation application identity", "[Jacobi functions]") {
	for (int variant : variants) {
		for (auto& tau : latticeApplicationTaus) {
			UNSCOPED_INFO("variant=" << variant << "  tau=" << tau);

			const auto control = LatticeTransform<double>{ variant, latticeApplicationZ, tau };
			const auto trial = RotateTau(variant, latticeApplicationZ, tau);
			const auto [rcontrol, rtrial] = LatticeApplicationIdentity(control, trial);

			REQUIRE(rcontrol == ApproxComplex(rtrial).epsilon(1e-3f));
		}
	}
}


//------------------------------------------------------------------------------
// Identities
//------------------------------------------------------------------------------


const std::array<std::complex<double>, 6> identityTaus = {
	0.73 + 1.49i,
	-3.17 + 0.49i,
	0.77i,
	-0.11 + 1.03i,
	-14.11 + 1.03i,
	7.83 + 1.03i,
};

const std::array<std::complex<double>, 7> identityZs = {
	-7.97,
	3.98 + 0.12i,
	-12.55,
	-0.12 - 0.77i,
	1.84 - 0.11i,
	0.798i,
	0.234,
};

// Lattice shift identities

template <class T>
std::pair<std::complex<T>, std::complex<T>> LatticeShiftIdentity1(const std::complex<T>& z, const std::complex<T>& tau) {
	const auto lhs = Theta(1, z, tau + T(1));
	const auto rhs = std::exp(pi_v<T> * i_v<T> / T(4)) * Theta(1, z, tau);
	return { lhs, rhs };
}

template <class T>
std::pair<std::complex<T>, std::complex<T>> LatticeShiftIdentity2(const std::complex<T>& z, const std::complex<T>& tau) {
	const auto lhs = Theta(2, z, tau + T(1));
	const auto rhs = std::exp(pi_v<T> * i_v<T> / T(4)) * Theta(2, z, tau);
	return { lhs, rhs };
}

template <class T>
std::pair<std::complex<T>, std::complex<T>> LatticeShiftIdentity3(const std::complex<T>& z, const std::complex<T>& tau) {
	const auto lhs = Theta(3, z, tau + T(1));
	const auto rhs = Theta(4, z, tau);
	return { lhs, rhs };
}

template <class T>
std::pair<std::complex<T>, std::complex<T>> LatticeShiftIdentity4(const std::complex<T>& z, const std::complex<T>& tau) {
	const auto lhs = Theta(4, z, tau + T(1));
	const auto rhs = Theta(3, z, tau);
	return { lhs, rhs };
}

TEST_CASE("Lattice shift identity #1", "[Jacobi functions]") {
	for (auto& z : identityZs) {
		for (auto& tau : identityTaus) {
			UNSCOPED_INFO("z=" << z << "  tau=" << tau);
			const auto [rcontrol, rtrial] = LatticeShiftIdentity1(z, tau);
			REQUIRE(rcontrol == ApproxComplex(rtrial).epsilon(1e-3f));
		}
	}
}

TEST_CASE("Lattice shift identity #2", "[Jacobi functions]") {
	for (auto& z : identityZs) {
		for (auto& tau : identityTaus) {
			UNSCOPED_INFO("z=" << z << "  tau=" << tau);
			const auto [rcontrol, rtrial] = LatticeShiftIdentity2(z, tau);
			REQUIRE(rcontrol == ApproxComplex(rtrial).epsilon(1e-3f));
		}
	}
}

TEST_CASE("Lattice shift identity #3", "[Jacobi functions]") {
	for (auto& z : identityZs) {
		for (auto& tau : identityTaus) {
			UNSCOPED_INFO("z=" << z << "  tau=" << tau);
			const auto [rcontrol, rtrial] = LatticeShiftIdentity3(z, tau);
			REQUIRE(rcontrol == ApproxComplex(rtrial).epsilon(1e-3f));
		}
	}
}

TEST_CASE("Lattice shift identity #4", "[Jacobi functions]") {
	for (auto& z : identityZs) {
		for (auto& tau : identityTaus) {
			UNSCOPED_INFO("z=" << z << "  tau=" << tau);
			const auto [rcontrol, rtrial] = LatticeShiftIdentity4(z, tau);
			REQUIRE(rcontrol == ApproxComplex(rtrial).epsilon(1e-3f));
		}
	}
}

// Lattice inversion identities identities

template <class T>
std::complex<T> LatticeInversionIdentityFactor(const std::complex<T>& z, const std::complex<T>& tau) {
	return std::sqrt(tau / i_v<T>) * std::exp(i_v<T> * tau * z * z / pi_v<T>);
}

template <class T>
std::pair<std::complex<T>, std::complex<T>> LatticeInversionIdentity1(const std::complex<T>& z, const std::complex<T>& tau) {
	const auto lhs = Theta(1, z, -T(1) / tau);
	const auto rhs = -i_v<T> * Theta(1, tau * z, tau) * LatticeInversionIdentityFactor(z, tau);
	return { lhs, rhs };
}

template <class T>
std::pair<std::complex<T>, std::complex<T>> LatticeInversionIdentity2(const std::complex<T>& z, const std::complex<T>& tau) {
	const auto lhs = Theta(2, z, -T(1) / tau);
	const auto rhs = Theta(4, tau * z, tau) * LatticeInversionIdentityFactor(z, tau);
	return { lhs, rhs };
}

template <class T>
std::pair<std::complex<T>, std::complex<T>> LatticeInversionIdentity3(const std::complex<T>& z, const std::complex<T>& tau) {
	const auto lhs = Theta(3, z, -T(1) / tau);
	const auto rhs = Theta(3, tau * z, tau) * LatticeInversionIdentityFactor(z, tau);
	return { lhs, rhs };
}

template <class T>
std::pair<std::complex<T>, std::complex<T>> LatticeInversionIdentity4(const std::complex<T>& z, const std::complex<T>& tau) {
	const auto lhs = Theta(4, z, -T(1) / tau);
	const auto rhs = Theta(2, tau * z, tau) * LatticeInversionIdentityFactor(z, tau);
	return { lhs, rhs };
}

TEST_CASE("Lattice inversion identity #1", "[Jacobi functions]") {
	for (auto& z : identityZs) {
		for (auto& tau : identityTaus) {
			UNSCOPED_INFO("z=" << z << "  tau=" << tau);
			const auto [rcontrol, rtrial] = LatticeInversionIdentity1(z, tau);
			REQUIRE(rcontrol == ApproxComplex(rtrial).epsilon(1e-3f));
		}
	}
}

TEST_CASE("Lattice inversion identity #2", "[Jacobi functions]") {
	for (auto& z : identityZs) {
		for (auto& tau : identityTaus) {
			UNSCOPED_INFO("z=" << z << "  tau=" << tau);
			const auto [rcontrol, rtrial] = LatticeInversionIdentity2(z, tau);
			REQUIRE(rcontrol == ApproxComplex(rtrial).epsilon(1e-3f));
		}
	}
}

TEST_CASE("Lattice inversion identity #3", "[Jacobi functions]") {
	for (auto& z : identityZs) {
		for (auto& tau : identityTaus) {
			UNSCOPED_INFO("z=" << z << "  tau=" << tau);
			const auto [rcontrol, rtrial] = LatticeInversionIdentity3(z, tau);
			REQUIRE(rcontrol == ApproxComplex(rtrial).epsilon(1e-3f));
		}
	}
}

TEST_CASE("Lattice inversion identity #4", "[Jacobi functions]") {
	for (auto& z : identityZs) {
		for (auto& tau : identityTaus) {
			UNSCOPED_INFO("z=" << z << "  tau=" << tau);
			const auto [rcontrol, rtrial] = LatticeInversionIdentity4(z, tau);
			REQUIRE(rcontrol == ApproxComplex(rtrial).epsilon(1e-3f));
		}
	}
}


// Periodicity

template <class T>
std::pair<std::complex<T>, std::complex<T>> PeriodicityIdentity(int variant, const std::complex<T>& z, const std::complex<T>& tau) {
	const T period = pi_v<T> * T(1 + (variant == 1 || variant == 2));
	const auto lhs = Theta(variant, z, tau);
	const auto rhs = Theta(variant, z + period, tau);
	return { lhs, rhs };
}

TEST_CASE("Periodicity identity", "[Jacobi functions]") {
	for (auto variant : variants) {
		for (auto& z : identityZs) {
			for (auto& tau : identityTaus) {
				UNSCOPED_INFO("variant=" << variant << "z=" << z << "  tau=" << tau);
				const auto [rcontrol, rtrial] = PeriodicityIdentity(variant, z, tau);
				REQUIRE(rcontrol == ApproxComplex(rtrial).epsilon(1e-3f));
			}
		}
	}
}

// Symmetry

template <class T>
std::pair<std::complex<T>, std::complex<T>> SymmetryIdentity(int variant, const std::complex<T>& z, const std::complex<T>& tau) {
	const auto lhs = Theta(variant, -z, tau);
	const auto rhs = variant == 1 ? -Theta(variant, z, tau) : Theta(variant, z, tau);
	return { lhs, rhs };
}

TEST_CASE("Symmetry identity", "[Jacobi functions]") {
	for (auto variant : variants) {
		for (auto& z : identityZs) {
			for (auto& tau : identityTaus) {
				UNSCOPED_INFO("variant=" << variant << " z=" << z << "  tau=" << tau);
				const auto [rcontrol, rtrial] = SymmetryIdentity(variant, z, tau);
				REQUIRE(rcontrol == ApproxComplex(rtrial).epsilon(1e-3f));
			}
		}
	}
}


// Conjugate symmetry

template <class T>
std::pair<std::complex<T>, std::complex<T>> ConjugateSymmetryIdentity(int variant, const std::complex<T>& z, const std::complex<T>& tau) {
	const auto lhs = Theta(variant, std::conj(z), tau);
	const auto rhs = std::conj(Theta(variant, z, -std::conj(tau)));
	return { lhs, rhs };
}

TEST_CASE("Conjugate symmetry identity", "[Jacobi functions]") {
	for (auto variant : variants) {
		for (auto& z : identityZs) {
			for (auto& tau : identityTaus) {
				UNSCOPED_INFO("variant=" << variant << " z=" << z << "  tau=" << tau);
				const auto [rcontrol, rtrial] = ConjugateSymmetryIdentity(variant, z, tau);
				REQUIRE(rcontrol == ApproxComplex(rtrial).epsilon(1e-3f));
			}
		}
	}
}


//------------------------------------------------------------------------------
// Randomized arguments
//------------------------------------------------------------------------------

std::mt19937_64 rne(7235472357);

template <class T>
std::vector<std::complex<T>> RandomTaus(size_t count) {
	std::uniform_real_distribution<T> rngMagnitude(T(0), T(0.5));
	std::uniform_real_distribution<T> rngPhase(-pi_v<T>, pi_v<T>);

	std::vector<std::complex<T>> taus;
	for (size_t i = 0; i < count; ++i) {
		const auto magnitude = rngMagnitude(rne);
		const auto phase = rngPhase(rne);
		const auto q = std::polar(magnitude, phase);
		const auto tau = -i_v<T> * std::log(q) / pi_v<T>;
		taus.push_back(tau);
	}

	return taus;
}

template <class T>
std::vector<std::complex<T>> RandomZs(size_t count) {
	std::uniform_real_distribution<T> rngReal(-T(3) * pi_v<T>, T(3) * pi_v<T>);
	std::uniform_real_distribution<T> rngImag(-T(1), T(1));

	std::vector<std::complex<T>> taus;
	for (size_t i = 0; i < count; ++i) {
		taus.emplace_back(rngReal(rne), rngImag(rne));
	}

	return taus;
}

TEST_CASE("Performance measure", "[Jacobi functions]") {
	using real_t = float;

	constexpr size_t count = 100'000;
	auto taus = RandomTaus<real_t>(count);
	auto zs = RandomZs<real_t>(count);

	const auto start = std::chrono::high_resolution_clock::now();
	float s = 0;
	for (size_t i = 0; i < count; ++i) {
		const auto r = Theta(1, zs[i], taus[i]);
		s += r.imag();
	}
	const auto end = std::chrono::high_resolution_clock::now();
	const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	const double nsPerOp = double(elapsed.count()) / double(count);
	WARN(nsPerOp << " ns / op");
	REQUIRE(s != 0.0f);
}

template <class T>
size_t HistogramBin(T value, T low, T high, size_t count) {
	const T t = (value - low) / (high - low);
	const intptr_t index = std::clamp(intptr_t(T(count - 1) * t + T(0.5)), intptr_t(0), intptr_t(count - 1));
	return size_t(index);
}

TEST_CASE("Accuracy measure", "[Jacobi functions]") {
	constexpr size_t count = 100'000;
	auto taus = RandomTaus<float>(count);
	auto zs = RandomZs<float>(count);
	std::vector<double> biterrs(count);
	std::vector<std::complex<double>> results(count);
	std::vector<std::complex<float>> resultfs(count);

	for (size_t i = 0; i < count; ++i) {
		const auto z = std::complex<double>(zs[i]);
		const auto tau = std::complex<double>(taus[i]);
		const auto zf = std::complex<float>(zs[i]);
		const auto tauf = std::complex<float>(taus[i]);

		const auto result = Theta(1, z, tau);
		const auto resultf = Theta(1, zf, tauf);
		results[i] = result;
		resultfs[i] = resultf;

		const auto err = result - std::complex<double>(resultf);
		const auto biterr = std::abs(err) / std::abs(result) / double(std::numeric_limits<float>::epsilon());
		if (abs(result) < 1000.0) {
			biterrs[i] = biterr;
		}
	}

	const auto view = AsView<DOMAINLESS>(biterrs.begin(), biterrs.end());
	const auto maxIt = std::max_element(view.begin(), view.end());
	WARN("outlier: " << *maxIt);
	WARN("z =       " << zs[maxIt - view.begin()]);
	WARN("tau =     " << taus[maxIt - view.begin()]);
	WARN("r_f32 =   " << resultfs[maxIt - view.begin()]);
	WARN("r_f64 =   " << results[maxIt - view.begin()]);
	[[maybe_unused]] const auto max = Max(Abs(view));
	[[maybe_unused]] const auto avg = Mean(Abs(view));
	[[maybe_unused]] const auto rms = RootMeanSquare(view);
	[[maybe_unused]] const auto std = CorrectedStandardDeviation(view);

	std::vector<double> histogram(100);
	for (size_t i = 0; i < count; ++i) {
		const double magQ = std::exp(-pi_v<double> * taus[i].imag());
		size_t bin = HistogramBin(magQ, 0.0, 1.0, histogram.size());
		histogram[bin] = std::min(100.0, std::max(histogram[bin], biterrs[i]));
	}

	REQUIRE(std::abs(max) < 100);
}

TEST_CASE("Accuracy debug", "[Jacobi functions]") {
	const auto z = -0.62151 - 0.535859i;
	const auto tau = -0.889212 + 0.00102383i;

	const auto result = Theta(1, z, tau);

	bool inf = std::isinf(std::abs(result));
	if (inf) {
		throw std::runtime_error("failed to eval");
	}
}


//------------------------------------------------------------------------------
// Debug
//------------------------------------------------------------------------------

TEST_CASE("Visualization debug", "[Jacobi functions]") {
	using real_t = double;

	const std::complex<real_t> q = -0.1f + 0.3if;
	const auto tau = -i_v<real_t> * std::log(q) / pi_v<real_t>;
	const auto q2 = std::exp(i_v<real_t> * pi_v<real_t> * tau);
	REQUIRE(q == ApproxComplex(q2));

	const auto xs = LinSpace<real_t, DOMAINLESS>(-2 * pi_v<real_t>, 2 * pi_v<real_t>, 400);
	std::vector<real_t> r;
	std::vector<real_t> i;
	for (auto& x : xs) {
		const auto c = theta(2, x + 0.0i, 0.73 + 1.49i);
		r.push_back(real(c));
		i.push_back(imag(c));
	}
	REQUIRE(r.size() == xs.Size());
}