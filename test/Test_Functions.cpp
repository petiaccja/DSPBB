#include <catch2/catch.hpp>
#include <cmath>
#include <complex>
#include <dspbb/Math/Functions.hpp>

using namespace dspbb;
using namespace std::complex_literals;

namespace std {
template <class T>
auto iden(T arg) {
	return arg;
}
} // namespace std


#define TEST_CASE_FUNCTION_REAL(NAME, FUNC, STDFUNC)                    \
	TEST_CASE(NAME " real", "[Functions]") {                            \
		const TimeSignal<float> signal = { 1, 8, 2, 5, 3, 6, 3, 6, 4 }; \
		const auto applied = FUNC(signal);                              \
		for (size_t i = 0; i < signal.Size(); ++i) {                    \
			REQUIRE(Approx(applied[i]) == std::STDFUNC(signal[i]));     \
		}                                                               \
	}


#define TEST_CASE_FUNCTION_CPLX(NAME, FUNC, STDFUNC)                                   \
	TEST_CASE(NAME " complex", "[Functions]") {                                        \
		const TimeSignal<std::complex<float>> csignal = { -1.f + 0.7if, 8.f + 2.6if }; \
		const auto capplied = FUNC(csignal);                                           \
		for (size_t i = 0; i < csignal.Size(); ++i) {                                  \
			REQUIRE(std::abs(capplied[i] - std::STDFUNC(csignal[i])) < 0.0001f);       \
		}                                                                              \
	}

// Complex number functions.
TEST_CASE_FUNCTION_REAL("Abs", Abs, abs);
TEST_CASE_FUNCTION_CPLX("Abs", Abs, abs);
TEST_CASE_FUNCTION_CPLX("Arg", Arg, arg);

TEST_CASE_FUNCTION_REAL("Real", Real, iden);
TEST_CASE_FUNCTION_CPLX("Real", Real, real);
TEST_CASE_FUNCTION_CPLX("Imag", Imag, imag);


// Exponential functions
TEST_CASE_FUNCTION_REAL("Log", Log, log);
TEST_CASE_FUNCTION_REAL("Log2", Log2, log2);
TEST_CASE_FUNCTION_REAL("Log10", Log10, log10);
TEST_CASE_FUNCTION_REAL("Exp", Exp, exp);


// Polynomial functions
TEST_CASE_FUNCTION_REAL("Sqrt", Sqrt, sqrt);
TEST_CASE_FUNCTION_CPLX("Sqrt", Sqrt, sqrt);
TEST_CASE_FUNCTION_REAL("Cbrt", Cbrt, cbrt);

TEST_CASE("Pow real", "[Functions]") {
	const TimeSignal<float> signal = { 1, 8 };
	const auto applied = Pow(signal, 2.5f);
	for (size_t i = 0; i < signal.Size(); ++i) {
		REQUIRE(Approx(applied[i]) == std::pow(signal[i], 2.5f));
	}
}
TEST_CASE("Pow complex", "[Functions]") {
	const TimeSignal<std::complex<float>> csignal = { -1.f + 0.7if, 8.f + 2.6if };
	const auto capplied = Pow(csignal, 2.5f);
	for (size_t i = 0; i < csignal.Size(); ++i) {
		REQUIRE(Approx(capplied[i].real()) == std::pow(csignal[i], 2.5f).real());
		REQUIRE(Approx(capplied[i].imag()) == std::pow(csignal[i], 2.5f).imag());
	}
}

// Trigonometric functions
TEST_CASE_FUNCTION_CPLX("Sin", Sin, sin);
TEST_CASE_FUNCTION_CPLX("Cos", Cos, cos);
TEST_CASE_FUNCTION_CPLX("Tan", Tan, tan);
TEST_CASE_FUNCTION_CPLX("Asin", Asin, asin);
TEST_CASE_FUNCTION_CPLX("Acos", Acos, acos);
TEST_CASE_FUNCTION_CPLX("Atan", Atan, atan);

// Hyperbolic functions
TEST_CASE_FUNCTION_CPLX("Sinh", Sinh, sinh);
TEST_CASE_FUNCTION_CPLX("Cosh", Cosh, cosh);
TEST_CASE_FUNCTION_CPLX("Tanh", Tanh, tanh);
TEST_CASE_FUNCTION_CPLX("Asinh", Asinh, asinh);
TEST_CASE_FUNCTION_CPLX("Acosh", Acosh, acosh);
TEST_CASE_FUNCTION_CPLX("Atanh", Atanh, atanh);