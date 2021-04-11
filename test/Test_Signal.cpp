#include <dspbb/Primitives/Signal.hpp>

#include <Catch2/catch.hpp>
#include <complex>

using namespace dspbb;
using namespace std::complex_literals;

TEST_CASE("Signal - Default construct", "[Signal]") {
	TimeSignal<float> s;
	TimeSignal<std::complex<float>> c;
	REQUIRE(s.Empty());
	REQUIRE(c.Empty());
}


TEST_CASE("Signal - Ilist construct", "[Signal]") {
	TimeSignal<float> s = { 1, 2, 3 };
	TimeSignal<std::complex<float>> c = { 1.f + 4.if, 2.f + 5.if, 3.f + 6.if };
	REQUIRE(s.Length() == 3);
	REQUIRE(c.Length() == 3);
	REQUIRE(!s.Empty());
	REQUIRE(!c.Empty());
	for (int i = 0; i < 3; ++i) {
		REQUIRE(s.Data()[i] == i + 1);
		REQUIRE(c.Data()[i].real() == i + 1);
		REQUIRE(c.Data()[i].imag() == i + 4);
	}
}


TEST_CASE("Signal - Element access", "[Signal]") {
	TimeSignal<float> s = { 1, 2, 3 };
	TimeSignal<std::complex<float>> c = { 1.f + 4.if, 2.f + 5.if, 3.f + 6.if };
	for (int i = 0; i < 3; ++i) {
		REQUIRE(s[i] == i + 1.f);
		REQUIRE(c[i] == float(i + 1) + float(i + 4) * 1.if);
	}
}


TEST_CASE("Signal - Conversion construct", "[Signal]") {
	TimeSignal<float> s = { 1, 2, 3 };
	TimeSignal<std::complex<float>> c = { 1.f + 4.if, 2.f + 5.if, 3.f + 6.if };
	TimeSignal<double> d = s;
	TimeSignal<std::complex<double>> cd = c;
	TimeSignal<std::complex<double>> cs = s;
	for (int i = 0; i < 3; ++i) {
		REQUIRE(d[i] == float(i) + 1.f);
		REQUIRE(cs[i] == double(i) + 1.0);
		REQUIRE(cd[i] == double(i + 1) + double(i + 4) * 1.i);
	}
}


TEST_CASE("Signal - Conversion operator=", "[Signal]") {
	TimeSignal<float> s = { 1, 2, 3 };
	TimeSignal<std::complex<float>> c = { 1.f + 4.if, 2.f + 5.if, 3.f + 6.if };
	TimeSignal<double> d;
	d = s;
	TimeSignal<std::complex<double>> cd;
	cd = c;
	TimeSignal<std::complex<double>> cs;
	cs = s;
	for (int i = 0; i < 3; ++i) {
		REQUIRE(d[i] == i + 1);
		REQUIRE(cs[i] == i + 1.);
		REQUIRE(cd[i] == double(i + 1) + double(i + 4) * 1.i);
	}
}


TEST_CASE("Signal - Reserve", "[Signal]") {
	TimeSignal<float> s = { 1, 2, 3 };
	TimeSignal<std::complex<float>> c = { 1.f + 4.if, 2.f + 5.if, 3.f + 6.if };
	s.Reserve(1024);
	c.Reserve(1024);
	REQUIRE(s.Capacity() >= 1024);
	REQUIRE(s.Size() == 3);
	REQUIRE(c.Capacity() >= 1024);
	REQUIRE(c.Size() == 3);
}


TEST_CASE("Signal - Resize", "[Signal]") {
	TimeSignal<float> s = { 1, 2, 3 };
	TimeSignal<std::complex<float>> c = { 1.f + 4.if, 2.f + 5.if, 3.f + 6.if };
	s.Resize(1024);
	c.Resize(1024);
	REQUIRE(s.Capacity() >= 1024);
	REQUIRE(s.Size() == 1024);
	REQUIRE(c.Capacity() >= 1024);
	REQUIRE(c.Size() == 1024);
}


TEST_CASE("Signal - Append", "[Signal]") {
	TimeSignal<float> s1 = { 1, 2, 3 };
	TimeSignal<float> s2 = { 4, 5, 6 };
	s1.Append(s2);
	REQUIRE(s2.Size() == 3);
	REQUIRE(s1.Size() == 6);
	REQUIRE(s1[2] == 3);
	REQUIRE(s1[3] == 4);
}


TEST_CASE("Signal - Prepend", "[Signal]") {
	TimeSignal<float> s1 = { 1, 2, 3 };
	TimeSignal<float> s2 = { 4, 5, 6 };
	s1.Prepend(s2);
	REQUIRE(s2.Size() == 3);
	REQUIRE(s1.Size() == 6);
	REQUIRE(s1[2] == 6);
	REQUIRE(s1[3] == 1);
}


TEST_CASE("Signal - ExtractFront", "[Signal]") {
	TimeSignal<float> s = { 1, 2, 3, 4, 5, 6 };
	TimeSignal<float> part = s.ExtractFront(2);
	REQUIRE(s.Size() == 4);
	REQUIRE(part.Size() == 2);
	REQUIRE(part[0] == 1);
	REQUIRE(part[1] == 2);
	REQUIRE(s[0] == 3);
	REQUIRE(s[1] == 4);
	REQUIRE(s[2] == 5);
	REQUIRE(s[3] == 6);
}


TEST_CASE("Signal - ExtractBack", "[Signal]") {
	TimeSignal<float> s = { 1, 2, 3, 4, 5, 6 };
	TimeSignal<float> part = s.ExtractBack(4);
	REQUIRE(s.Size() == 2);
	REQUIRE(part.Size() == 4);
	REQUIRE(s[0] == 1);
	REQUIRE(s[1] == 2);
	REQUIRE(part[0] == 3);
	REQUIRE(part[1] == 4);
	REQUIRE(part[2] == 5);
	REQUIRE(part[3] == 6);
}


TEST_CASE("Signal - Erase", "[Signal]") {
	TimeSignal<float> s = { 1, 2, 3, 4, 5, 6 };
	s.Erase(s.begin() + 3);
	REQUIRE(s.Size() == 5);
	REQUIRE(s[2] == 3);
	REQUIRE(s[3] == 5);
}

TEST_CASE("Signal - Erase range", "[Signal]") {
	TimeSignal<float> s = { 1, 2, 3, 4, 5, 6 };
	s.Erase(s.begin() + 1, s.begin() + 5);
	REQUIRE(s.Size() == 2);
	REQUIRE(s[0] == 1);
	REQUIRE(s[1] == 6);
}


TEST_CASE("Signal - Append (complex)", "[Signal]") {
	TimeSignal<std::complex<float>> s1 = { 1, 2, 3 };
	TimeSignal<float> s2 = { 4, 5, 6 };
	s1.Append(s2);
	REQUIRE(s2.Size() == 3);
	REQUIRE(s1.Size() == 6);
	REQUIRE(s1[2] == 3.f);
	REQUIRE(s1[3] == 4.f);
}


TEST_CASE("Signal - Prepend (complex)", "[Signal]") {
	TimeSignal<std::complex<float>> s1 = { 1, 2, 3 };
	TimeSignal<float> s2 = { 4, 5, 6 };
	s1.Prepend(s2);
	REQUIRE(s2.Size() == 3);
	REQUIRE(s1.Size() == 6);
	REQUIRE(s1[2] == 6.f);
	REQUIRE(s1[3] == 1.f);
}


TEST_CASE("Signal - ExtractFront (complex)", "[Signal]") {
	TimeSignal<std::complex<float>> s = { 1.f * (1.f + 1.if), 2.f * (1.f + 1.if), 3.f * (1.f + 1.if), 4.f * (1.f + 1.if), 5.f * (1.f + 1.if), 6.f * (1.f + 1.if) };
	TimeSignal<std::complex<float>> part = s.ExtractFront(2);
	REQUIRE(s.Size() == 4);
	REQUIRE(part.Size() == 2);
	REQUIRE(part[0] == 1.f * (1.f + 1.if));
	REQUIRE(part[1] == 2.f * (1.f + 1.if));
	REQUIRE(s[0] == 3.f * (1.f + 1.if));
	REQUIRE(s[1] == 4.f * (1.f + 1.if));
	REQUIRE(s[2] == 5.f * (1.f + 1.if));
	REQUIRE(s[3] == 6.f * (1.f + 1.if));
}


TEST_CASE("Signal - ExtractBack (complex)", "[Signal]") {
	TimeSignal<std::complex<float>> s = { 1.f * (1.f + 1.if), 2.f * (1.f + 1.if), 3.f * (1.f + 1.if), 4.f * (1.f + 1.if), 5.f * (1.f + 1.if), 6.f * (1.f + 1.if) };
	TimeSignal<std::complex<float>> part = s.ExtractBack(4);
	REQUIRE(s.Size() == 2);
	REQUIRE(part.Size() == 4);
	REQUIRE(s[0] == 1.f * (1.f + 1.if));
	REQUIRE(s[1] == 2.f * (1.f + 1.if));
	REQUIRE(part[0] == 3.f * (1.f + 1.if));
	REQUIRE(part[1] == 4.f * (1.f + 1.if));
	REQUIRE(part[2] == 5.f * (1.f + 1.if));
	REQUIRE(part[3] == 6.f * (1.f + 1.if));
}


TEST_CASE("Signal - Erase (complex)", "[Signal]") {
	TimeSignal<std::complex<float>> s = { 1, 2, 3, 4, 5, 6 };
	s.Erase(s.begin() + 3);
	REQUIRE(s.Size() == 5);
	REQUIRE(s[2] == 3.f);
	REQUIRE(s[3] == 5.f);
}


TEST_CASE("Signal - Erase range (complex)", "[Signal]") {
	TimeSignal<std::complex<float>> s = { 1, 2, 3, 4, 5, 6 };
	s.Erase(s.begin() + 1, s.begin() + 5);
	REQUIRE(s.Size() == 2);
	REQUIRE(s[0] == 1.f);
	REQUIRE(s[1] == 6.f);
}


TEST_CASE("Signal - Iteration", "[Signal]") {
	TimeSignal<float> s = { 1, 2, 3, 4, 5, 6 };
	float expected = 1.0f;
	for (const auto& v : s) {
		REQUIRE(v == expected);
		expected += 1.0f;
	}
}


TEST_CASE("Signal - Iteration (complex)", "[Signal]") {
	TimeSignal<std::complex<float>> s = { 1, 2, 3, 4, 5, 6 };
	std::complex<float> expected = 1.0f;
	for (const auto& v : s) {
		REQUIRE(v == expected);
		expected += 1.0f;
	}
}


TEST_CASE("Signal - Const iteration", "[Signal]") {
	const TimeSignal<float> s = { 1, 2, 3, 4, 5, 6 };
	float expected = 1.0f;
	for (const auto& v : s) {
		REQUIRE(v == expected);
		expected += 1.0f;
	}
}


TEST_CASE("Signal - Const iteration (complex)", "[Signal]") {
	const TimeSignal<std::complex<float>> s = { 1, 2, 3, 4, 5, 6 };
	std::complex<float> expected = 1.0f;
	for (const auto& v : s) {
		REQUIRE(v == expected);
		expected += 1.0f;
	}
}


//------------------------------------------------------------------------------
// Real operators
//------------------------------------------------------------------------------
/*
TEST_CASE("Signal - += signal", "[Signal]") {
	TimeSignal<float> s = { 1, 2, 3 };
	TimeSignal<float> a = { 3, 2, 1 };
	s += a;
	for (auto&& v : s) {
		REQUIRE(v == 4.0f);
	}
}

TEST_CASE("Signal - -= signal", "[Signal]") {
	TimeSignal<float> s = { 1, 2, 3 };
	TimeSignal<float> a = { 4, 5, 6 };
	s -= a;
	for (auto&& v : s) {
		REQUIRE(v == -3.0f);
	}
}

TEST_CASE("Signal - *= signal", "[Signal]") {
	TimeSignal<float> s = { 1, 2, 4 };
	TimeSignal<float> a = { 4, 2, 1 };
	s *= a;
	for (auto&& v : s) {
		REQUIRE(v == 4.0f);
	}
}

TEST_CASE("Signal - /= signal", "[Signal]") {
	TimeSignal<float> s = { 1, 2, 3 };
	TimeSignal<float> a = { 2, 4, 6 };
	s /= a;
	for (auto&& v : s) {
		REQUIRE(v == 0.5f);
	}
}


TEST_CASE("Signal - += scalar", "[Signal]") {
	TimeSignal<float> s = { 1, 1, 1 };
	s += 5.f;
	for (auto&& v : s) {
		REQUIRE(v == 6.f);
	}
}

TEST_CASE("Signal - -= scalar", "[Signal]") {
	TimeSignal<float> s = { 6, 6, 6 };
	s -= 5.0f;
	for (auto&& v : s) {
		REQUIRE(v == 1.0f);
	}
}

TEST_CASE("Signal - *= scalar", "[Signal]") {
	TimeSignal<float> s = { 2, 2, 2 };
	s *= 3.0;
	for (auto&& v : s) {
		REQUIRE(v == 6.0f);
	}
}

TEST_CASE("Signal - /= scalar", "[Signal]") {
	TimeSignal<float> s = { 2, 2, 2 };
	s /= 4.0f;
	for (auto&& v : s) {
		REQUIRE(v == 0.5f);
	}
}



//------------------------------------------------------------------------------
// Complex-real operators
//------------------------------------------------------------------------------

TEST_CASE("Signal - += signal (complex-real)", "[Signal]") {
	TimeSignal<std::complex<float>> s = { 1.f + 1.if, 2.f + 1.if, 3.f + 1.if };
	TimeSignal<float> a = { 3, 2, 1 };
	s += a;
	for (auto&& v : s) {
		REQUIRE(v.real() == 4.0f);
		REQUIRE(v.imag() == 1.0f);
	}
}

TEST_CASE("Signal - -= signal (complex-real)", "[Signal]") {
	TimeSignal<std::complex<float>> s = { 1.f + 1.if, 2.f + 1.if, 3.f + 1.if };
	TimeSignal<float> a = { 4, 5, 6 };
	s -= a;
	for (auto&& v : s) {
		REQUIRE(v.real() == -3.0f);
		REQUIRE(v.imag() == 1.0f);
	}
}

TEST_CASE("Signal - *= signal (complex-real)", "[Signal]") {
	TimeSignal<std::complex<float>> s = { 1.f + 1.if, 2.f + 2.if, 4.f + 4.if };
	TimeSignal<float> a = { 4, 2, 1 };
	s *= a;
	for (auto&& v : s) {
		REQUIRE(v.real() == 4.0f);
		REQUIRE(v.imag() == 4.0f);
	}
}

TEST_CASE("Signal - /= signal (complex-real)", "[Signal]") {
	TimeSignal<std::complex<float>> s = { 1.f + 1.if, 2.f + 2.if, 3.f + 3.if };
	TimeSignal<float> a = { 2, 4, 6 };
	s /= a;
	for (auto&& v : s) {
		REQUIRE(v.real() == 0.5f);
		REQUIRE(v.imag() == 0.5f);
	}
}


TEST_CASE("Signal - += scalar (complex-real)", "[Signal]") {
	TimeSignal<std::complex<float>> s = { 1.f + 1.if, 1.f + 1.if, 1.f + 1.if };
	s += 5.f;
	for (auto&& v : s) {
		REQUIRE(v.real() == 6.f);
		REQUIRE(v.imag() == 1.f);
	}
}

TEST_CASE("Signal - -= scalar (complex-real)", "[Signal]") {
	TimeSignal<std::complex<float>> s = { 6.f + 1.if, 6.f + 1.if, 6.f + 1.if };
	s -= 5.0f;
	for (auto&& v : s) {
		REQUIRE(v.real() == 1.0f);
		REQUIRE(v.imag() == 1.0f);
	}
}

TEST_CASE("Signal - *= scalar (complex-real)", "[Signal]") {
	TimeSignal<std::complex<float>> s = { 2.f + 1.if, 2.f + 1.if, 2.f + 1.if };
	s *= 3.0;
	for (auto&& v : s) {
		REQUIRE(v.real() == 6.0f);
		REQUIRE(v.imag() == 3.0f);
	}
}

TEST_CASE("Signal - /= scalar (complex-real)", "[Signal]") {
	TimeSignal<std::complex<float>> s = { 2.f + 4.if, 2.f + 4.if, 2.f + 4.if };
	s /= 4.0f;
	for (auto&& v : s) {
		REQUIRE(v.real() == 0.5f);
		REQUIRE(v.imag() == 1.f);
	}
}



//------------------------------------------------------------------------------
// Complex-complex operators
//------------------------------------------------------------------------------

TEST_CASE("Signal - += signal (complex-complex)", "[Signal]") {
	const TimeSignal<std::complex<float>> a = { 1.f + 1.if, 2.f + 1.if, 3.f + 1.if };
	const TimeSignal<std::complex<float>> b = { 3.f + 3.if, 2.f + 6.if, 1.f + 1.if };
	auto copy = a;
	copy += b;

	for (int i = 0; i < 3; ++i) {
		REQUIRE(copy[i] == a[i] + b[i]);
	}
}

TEST_CASE("Signal - -= signal (complex-complex)", "[Signal]") {
	const TimeSignal<std::complex<float>> a = { 1.f + 1.if, 2.f + 1.if, 3.f + 1.if };
	const TimeSignal<std::complex<float>> b = { 3.f + 3.if, 2.f + 6.if, 1.f + 1.if };
	auto copy = a;
	copy -= b;

	for (int i = 0; i < 3; ++i) {
		REQUIRE(copy[i] == a[i] - b[i]);
	}
}

TEST_CASE("Signal - *= signal (complex-complex)", "[Signal]") {
	const TimeSignal<std::complex<float>> a = { 1.f + 1.if, 2.f + 1.if, 3.f + 1.if };
	const TimeSignal<std::complex<float>> b = { 3.f + 3.if, 2.f + 6.if, 1.f + 1.if };
	auto copy = a;
	copy *= b;

	for (int i = 0; i < 3; ++i) {
		REQUIRE(copy[i] == a[i] * b[i]);
	}
}

TEST_CASE("Signal - /= signal (complex-complex)", "[Signal]") {
	const TimeSignal<std::complex<float>> a = { 1.f + 1.if, 2.f + 1.if, 3.f + 1.if };
	const TimeSignal<std::complex<float>> b = { 3.f + 3.if, 2.f + 6.if, 1.f + 1.if };
	auto copy = a;
	copy /= b;

	for (int i = 0; i < 3; ++i) {
		REQUIRE(Approx(copy[i].real()) == (a[i] / b[i]).real());
		REQUIRE(Approx(copy[i].imag()) == (a[i] / b[i]).imag());
	}
}


TEST_CASE("Signal - += scalar (complex-complex)", "[Signal]") {
	const TimeSignal<std::complex<float>> a = { 1.f + 1.if, 2.f + 1.if, 3.f + 1.if };
	const std::complex<float> b = 3.f + 5.if;
	auto copy = a;
	copy += b;

	for (int i = 0; i < 3; ++i) {
		REQUIRE(copy[i] == a[i] + b);
	}
}

TEST_CASE("Signal - -= scalar (complex-complex)", "[Signal]") {
	const TimeSignal<std::complex<float>> a = { 1.f + 1.if, 2.f + 1.if, 3.f + 1.if };
	const std::complex<float> b = 3.f + 5.if;
	auto copy = a;
	copy -= b;

	for (int i = 0; i < 3; ++i) {
		REQUIRE(copy[i] == a[i] - b);
	}
}

TEST_CASE("Signal - *= scalar (complex-complex)", "[Signal]") {
	const TimeSignal<std::complex<float>> a = { 1.f + 1.if, 2.f + 1.if, 3.f + 1.if };
	const std::complex<float> b = 3.f + 5.if;
	auto copy = a;
	copy *= b;

	for (int i = 0; i < 3; ++i) {
		REQUIRE(copy[i] == a[i] * b);
	}
}

TEST_CASE("Signal - /= scalar (complex-complex)", "[Signal]") {
	const TimeSignal<std::complex<float>> a = { 1.f + 1.if, 2.f + 1.if, 3.f + 1.if };
	const std::complex<float> b = 3.f + 5.if;
	auto copy = a;
	copy /= b;

	for (int i = 0; i < 3; ++i) {
		REQUIRE(Approx(copy[i].real()) == (a[i] / b).real());
		REQUIRE(Approx(copy[i].imag()) == (a[i] / b).imag());
	}
}
*/