#include <dspbb/Math/Functions.hpp>
#include <dspbb/Primitives/Signal.hpp>

#include <catch2/catch.hpp>
#include <cmath>
#include <complex>

using namespace dspbb;
using namespace std::complex_literals;

template <class T>
auto iden(T arg) {
	return arg;
}


#define TEST_CASE_FUNCTION_REAL(NAME, FUNC, STDFUNC)                                           \
	TEST_CASE(NAME " real", "[Functions]") {                                                   \
		using namespace std;                                                                   \
		const Signal<float> signal = { 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f }; \
		const auto applied = FUNC(signal);                                                     \
		for (size_t i = 0; i < signal.size(); ++i) {                                           \
			REQUIRE(Approx(applied[i]) == STDFUNC(signal[i]));                                 \
		}                                                                                      \
	}


#define TEST_CASE_FUNCTION_CPLX(NAME, FUNC, STDFUNC)                               \
	TEST_CASE(NAME " complex", "[Functions]") {                                    \
		using namespace std;                                                       \
		const Signal<std::complex<float>> csignal = { -1.f + 0.7if, 8.f + 2.6if }; \
		const auto capplied = FUNC(csignal);                                       \
		for (size_t i = 0; i < csignal.size(); ++i) {                              \
			REQUIRE(std::abs(capplied[i] - STDFUNC(csignal[i])) < 0.0001f);        \
		}                                                                          \
	}

// Complex number functions.
TEST_CASE_FUNCTION_REAL("Abs", Abs, abs);
TEST_CASE_FUNCTION_CPLX("Abs", Abs, abs);
TEST_CASE_FUNCTION_REAL("Arg", Arg, arg);
TEST_CASE_FUNCTION_CPLX("Arg", Arg, arg);

TEST_CASE_FUNCTION_REAL("Real", Real, iden);
TEST_CASE_FUNCTION_CPLX("Real", Real, real);
TEST_CASE_FUNCTION_REAL("Imag", Imag, imag);
TEST_CASE_FUNCTION_CPLX("Imag", Imag, imag);

TEST_CASE_FUNCTION_CPLX("Conj", Conj, conj);

// Exponential functions
TEST_CASE_FUNCTION_REAL("Log", Log, log);
TEST_CASE_FUNCTION_REAL("Log2", Log2, log2);
TEST_CASE_FUNCTION_REAL("Log10", Log10, log10);
TEST_CASE_FUNCTION_REAL("Exp", Exp, exp);

TEST_CASE_FUNCTION_CPLX("Log", Log, log);
TEST_CASE_FUNCTION_CPLX("Log10", Log10, log10);
TEST_CASE_FUNCTION_CPLX("Exp", Exp, exp);


// Polynomial functions
TEST_CASE_FUNCTION_REAL("Sqrt", Sqrt, sqrt);
TEST_CASE_FUNCTION_CPLX("Sqrt", Sqrt, sqrt);
TEST_CASE_FUNCTION_REAL("Cbrt", Cbrt, cbrt);

TEST_CASE("Pow real", "[Functions]") {
	const Signal<float> signal = { 1, 8 };
	const auto applied = Pow(signal, 2.5f);
	for (size_t i = 0; i < signal.size(); ++i) {
		REQUIRE(Approx(applied[i]) == std::pow(signal[i], 2.5f));
	}
}
TEST_CASE("Pow complex", "[Functions]") {
	const Signal<std::complex<float>> csignal = { -1.f + 0.7if, 8.f + 2.6if };
	const auto capplied = Pow(csignal, 2.5f);
	for (size_t i = 0; i < csignal.size(); ++i) {
		REQUIRE(Approx(capplied[i].real()) == std::pow(csignal[i], 2.5f).real());
		REQUIRE(Approx(capplied[i].imag()) == std::pow(csignal[i], 2.5f).imag());
	}
}

// Trigonometric functions
TEST_CASE_FUNCTION_REAL("Sin", Sin, sin);
TEST_CASE_FUNCTION_REAL("Cos", Cos, cos);
TEST_CASE_FUNCTION_REAL("Tan", Tan, tan);
TEST_CASE_FUNCTION_REAL("Asin", Asin, asin);
TEST_CASE_FUNCTION_REAL("Acos", Acos, acos);
TEST_CASE_FUNCTION_REAL("Atan", Atan, atan);

TEST_CASE_FUNCTION_CPLX("Sin", Sin, sin);
TEST_CASE_FUNCTION_CPLX("Cos", Cos, cos);
TEST_CASE_FUNCTION_CPLX("Tan", Tan, tan);
TEST_CASE_FUNCTION_CPLX("Asin", Asin, asin);
TEST_CASE_FUNCTION_CPLX("Acos", Acos, acos);
TEST_CASE_FUNCTION_CPLX("Atan", Atan, atan);

// Hyperbolic functions
TEST_CASE_FUNCTION_REAL("Sinh", Sinh, sinh);
TEST_CASE_FUNCTION_REAL("Cosh", Cosh, cosh);
TEST_CASE_FUNCTION_REAL("Tanh", Tanh, tanh);
TEST_CASE_FUNCTION_REAL("Asinh", Asinh, asinh);
TEST_CASE_FUNCTION_REAL("Atanh", Atanh, atanh);

TEST_CASE_FUNCTION_CPLX("Sinh", Sinh, sinh);
TEST_CASE_FUNCTION_CPLX("Cosh", Cosh, cosh);
TEST_CASE_FUNCTION_CPLX("Tanh", Tanh, tanh);
TEST_CASE_FUNCTION_CPLX("Asinh", Asinh, asinh);
TEST_CASE_FUNCTION_CPLX("Acosh", Acosh, acosh);
TEST_CASE_FUNCTION_CPLX("Atanh", Atanh, atanh);

// Hyperbolic functions
TEST_CASE_FUNCTION_REAL("Erf", Erf, erf);
TEST_CASE_FUNCTION_REAL("Erfc", Erfc, erfc);
TEST_CASE_FUNCTION_REAL("TGamma", TGamma, tgamma);
TEST_CASE_FUNCTION_REAL("LGamma", LGamma, lgamma);