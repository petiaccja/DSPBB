#include <dspbb/Primitives/Signal.hpp>

#include <catch2/catch.hpp>
#include <complex>

using namespace dspbb;
using namespace std::complex_literals;

TEST_CASE("Signal - Default construct", "[Signal]") {
	Signal<float> s;
	Signal<std::complex<float>> c;
	REQUIRE(s.Empty());
	REQUIRE(c.Empty());
}


TEST_CASE("Signal - Ilist construct", "[Signal]") {
	Signal<float> s = { 1, 2, 3 };
	Signal<std::complex<float>> c = { 1.f + 4.if, 2.f + 5.if, 3.f + 6.if };
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
	Signal<float> s = { 1, 2, 3 };
	Signal<std::complex<float>> c = { 1.f + 4.if, 2.f + 5.if, 3.f + 6.if };
	for (int i = 0; i < 3; ++i) {
		REQUIRE(s[i] == i + 1.f);
		REQUIRE(c[i] == float(i + 1) + float(i + 4) * 1.if);
	}
}


TEST_CASE("Signal - Conversion construct", "[Signal]") {
	Signal<float> s = { 1, 2, 3 };
	Signal<std::complex<float>> c = { 1.f + 4.if, 2.f + 5.if, 3.f + 6.if };
	Signal<double> d{ s };
	Signal<std::complex<double>> cd{ c };
	Signal<std::complex<double>> cs{ s };
	for (int i = 0; i < 3; ++i) {
		REQUIRE(d[i] == float(i) + 1.f);
		REQUIRE(cs[i] == double(i) + 1.0);
		REQUIRE(cd[i] == double(i + 1) + double(i + 4) * 1.i);
	}
}


TEST_CASE("Signal - Conversion operator=", "[Signal]") {
	Signal<float> s = { 1, 2, 3 };
	Signal<std::complex<float>> c = { 1.f + 4.if, 2.f + 5.if, 3.f + 6.if };
	Signal<double> d;
	d = s;
	Signal<std::complex<double>> cd;
	cd = c;
	Signal<std::complex<double>> cs;
	cs = s;
	for (int i = 0; i < 3; ++i) {
		REQUIRE(d[i] == i + 1);
		REQUIRE(cs[i] == i + 1.);
		REQUIRE(cd[i] == double(i + 1) + double(i + 4) * 1.i);
	}
}


TEST_CASE("Signal - Reserve", "[Signal]") {
	Signal<float> s = { 1, 2, 3 };
	Signal<std::complex<float>> c = { 1.f + 4.if, 2.f + 5.if, 3.f + 6.if };
	s.Reserve(1024);
	c.Reserve(1024);
	REQUIRE(s.Capacity() >= 1024);
	REQUIRE(s.Size() == 3);
	REQUIRE(c.Capacity() >= 1024);
	REQUIRE(c.Size() == 3);
}


TEST_CASE("Signal - Resize", "[Signal]") {
	Signal<float> s = { 1, 2, 3 };
	Signal<std::complex<float>> c = { 1.f + 4.if, 2.f + 5.if, 3.f + 6.if };
	s.Resize(1024);
	c.Resize(1024);
	REQUIRE(s.Capacity() >= 1024);
	REQUIRE(s.Size() == 1024);
	REQUIRE(c.Capacity() >= 1024);
	REQUIRE(c.Size() == 1024);
}


TEST_CASE("Signal - Append", "[Signal]") {
	Signal<float> s1 = { 1, 2, 3 };
	Signal<float> s2 = { 4, 5, 6 };
	s1.Append(s2);
	REQUIRE(s2.Size() == 3);
	REQUIRE(s1.Size() == 6);
	REQUIRE(s1[2] == 3);
	REQUIRE(s1[3] == 4);
}


TEST_CASE("Signal - Prepend", "[Signal]") {
	Signal<float> s1 = { 1, 2, 3 };
	Signal<float> s2 = { 4, 5, 6 };
	s1.Prepend(s2);
	REQUIRE(s2.Size() == 3);
	REQUIRE(s1.Size() == 6);
	REQUIRE(s1[2] == 6);
	REQUIRE(s1[3] == 1);
}


TEST_CASE("Signal - ExtractFront", "[Signal]") {
	Signal<float> s = { 1, 2, 3, 4, 5, 6 };
	Signal<float> part = s.ExtractFront(2);
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
	Signal<float> s = { 1, 2, 3, 4, 5, 6 };
	Signal<float> part = s.ExtractBack(4);
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
	Signal<float> s = { 1, 2, 3, 4, 5, 6 };
	s.Erase(s.begin() + 3);
	REQUIRE(s.Size() == 5);
	REQUIRE(s[2] == 3);
	REQUIRE(s[3] == 5);
}

TEST_CASE("Signal - Erase range", "[Signal]") {
	Signal<float> s = { 1, 2, 3, 4, 5, 6 };
	s.Erase(s.begin() + 1, s.begin() + 5);
	REQUIRE(s.Size() == 2);
	REQUIRE(s[0] == 1);
	REQUIRE(s[1] == 6);
}


TEST_CASE("Signal - Append (complex)", "[Signal]") {
	Signal<std::complex<float>> s1 = { 1, 2, 3 };
	Signal<std::complex<float>> s2 = { 4, 5, 6 };
	s1.Append(s2);
	REQUIRE(s2.Size() == 3);
	REQUIRE(s1.Size() == 6);
	REQUIRE(s1[2] == 3.f);
	REQUIRE(s1[3] == 4.f);
}


TEST_CASE("Signal - Prepend (complex)", "[Signal]") {
	Signal<std::complex<float>> s1 = { 1, 2, 3 };
	Signal<std::complex<float>> s2 = { 4, 5, 6 };
	s1.Prepend(s2);
	REQUIRE(s2.Size() == 3);
	REQUIRE(s1.Size() == 6);
	REQUIRE(s1[2] == 6.f);
	REQUIRE(s1[3] == 1.f);
}


TEST_CASE("Signal - ExtractFront (complex)", "[Signal]") {
	Signal<std::complex<float>> s = { 1.f * (1.f + 1.if), 2.f * (1.f + 1.if), 3.f * (1.f + 1.if), 4.f * (1.f + 1.if), 5.f * (1.f + 1.if), 6.f * (1.f + 1.if) };
	Signal<std::complex<float>> part = s.ExtractFront(2);
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
	Signal<std::complex<float>> s = { 1.f * (1.f + 1.if), 2.f * (1.f + 1.if), 3.f * (1.f + 1.if), 4.f * (1.f + 1.if), 5.f * (1.f + 1.if), 6.f * (1.f + 1.if) };
	Signal<std::complex<float>> part = s.ExtractBack(4);
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
	Signal<std::complex<float>> s = { 1, 2, 3, 4, 5, 6 };
	s.Erase(s.begin() + 3);
	REQUIRE(s.Size() == 5);
	REQUIRE(s[2] == 3.f);
	REQUIRE(s[3] == 5.f);
}


TEST_CASE("Signal - Erase range (complex)", "[Signal]") {
	Signal<std::complex<float>> s = { 1, 2, 3, 4, 5, 6 };
	s.Erase(s.begin() + 1, s.begin() + 5);
	REQUIRE(s.Size() == 2);
	REQUIRE(s[0] == 1.f);
	REQUIRE(s[1] == 6.f);
}


TEST_CASE("Signal - Iteration", "[Signal]") {
	Signal<float> s = { 1, 2, 3, 4, 5, 6 };
	float expected = 1.0f;
	for (const auto& v : s) {
		REQUIRE(v == expected);
		expected += 1.0f;
	}
}


TEST_CASE("Signal - Iteration (complex)", "[Signal]") {
	Signal<std::complex<float>> s = { 1, 2, 3, 4, 5, 6 };
	std::complex<float> expected = 1.0f;
	for (const auto& v : s) {
		REQUIRE(v == expected);
		expected += 1.0f;
	}
}


TEST_CASE("Signal - Const iteration", "[Signal]") {
	const Signal<float> s = { 1, 2, 3, 4, 5, 6 };
	float expected = 1.0f;
	for (const auto& v : s) {
		REQUIRE(v == expected);
		expected += 1.0f;
	}
}


TEST_CASE("Signal - Const iteration (complex)", "[Signal]") {
	const Signal<std::complex<float>> s = { 1, 2, 3, 4, 5, 6 };
	std::complex<float> expected = 1.0f;
	for (const auto& v : s) {
		REQUIRE(v == expected);
		expected += 1.0f;
	}
}