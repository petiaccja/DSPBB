#include <dspbb/Math/Rational.hpp>

#include <catch2/catch.hpp>
#include <random>


using namespace dspbb;

TEST_CASE("Rational create whole", "[Rational]") {
	constexpr Rational r{ 6 };
	REQUIRE(r.Numerator() == 6);
	REQUIRE(r.Denominator() == 1);
}

TEST_CASE("Rational create num,den", "[Rational]") {
	constexpr Rational r{ 1, 3 };
	REQUIRE(r.Numerator() == 1);
	REQUIRE(r.Denominator() == 3);
}


TEST_CASE("Rational increment postfix", "[Rational]") {
	Rational r{ 5, 7 };
	const auto q = r++;
	REQUIRE(q == Rational{ 5, 7 });
	REQUIRE(r == Rational{ 12, 7 });
}

TEST_CASE("Rational decrement postfix", "[Rational]") {
	Rational r{ 5, 7 };
	const auto q = r--;
	REQUIRE(q == Rational{ 5, 7 });
	REQUIRE(r == Rational{ -2, 7 });
}

TEST_CASE("Rational increment prefix", "[Rational]") {
	Rational r{ 5, 7 };
	const auto q = ++r;
	REQUIRE(q == r);
	REQUIRE(r == Rational{ 12, 7 });
}

TEST_CASE("Rational decrement prefix", "[Rational]") {
	Rational r{ 5, 7 };
	const auto q = --r;
	REQUIRE(q == r);
	REQUIRE(r == Rational{ -2, 7 });
}


TEST_CASE("Rational add", "[Rational]") {
	constexpr Rational a{ 5, 7 };
	constexpr Rational b{ 125, 98 };
	constexpr auto r = a + b;
	REQUIRE(float(a) + float(b) == Approx(r));
}

TEST_CASE("Rational sub", "[Rational]") {
	constexpr Rational a{ 5, 7 };
	constexpr Rational b{ 125, 98 };
	constexpr auto r = a - b;
	REQUIRE(float(a) - float(b) == Approx(r));
}

TEST_CASE("Rational mul", "[Rational]") {
	constexpr Rational a{ 5, 7 };
	constexpr Rational b{ 125, 98 };
	constexpr auto r = a * b;
	REQUIRE(float(a) * float(b) == Approx(r));
}

TEST_CASE("Rational div", "[Rational]") {
	constexpr Rational a{ 5, 7 };
	constexpr Rational b{ 125, 98 };
	constexpr auto r = a / b;
	REQUIRE(float(a) / float(b) == Approx(r));
}


TEST_CASE("Rational add/int", "[Rational]") {
	constexpr Rational a{ 5, 7 };
	constexpr auto b = 125;
	constexpr auto r = a + b;
	REQUIRE(float(a) + float(b) == Approx(r));
}

TEST_CASE("Rational sub/int", "[Rational]") {
	constexpr Rational a{ 5, 7 };
	constexpr auto b = 125;
	constexpr auto r = a - b;
	REQUIRE(float(a) - float(b) == Approx(r));
}

TEST_CASE("Rational mul/int", "[Rational]") {
	constexpr Rational a{ 5, 7 };
	constexpr auto b = 125;
	constexpr auto r = a * b;
	REQUIRE(float(a) * float(b) == Approx(r));
}

TEST_CASE("Rational div/int", "[Rational]") {
	constexpr Rational a{ 5, 7 };
	constexpr auto b = 125;
	constexpr auto r = a / b;
	REQUIRE(float(a) / float(b) == Approx(r));
}

TEST_CASE("Rational int/add", "[Rational]") {
	constexpr Rational b{ 5, 7 };
	constexpr auto a = 125;
	constexpr auto r = a + b;
	REQUIRE(float(a) + float(b) == Approx(r));
}

TEST_CASE("Rational int/sub", "[Rational]") {
	constexpr Rational b{ 5, 7 };
	constexpr auto a = 125;
	constexpr auto r = a - b;
	REQUIRE(float(a) - float(b) == Approx(r));
}

TEST_CASE("Rational int/mul", "[Rational]") {
	constexpr Rational b{ 5, 7 };
	constexpr auto a = 125;
	constexpr auto r = a * b;
	REQUIRE(float(a) * float(b) == Approx(r));
}

TEST_CASE("Rational int/div", "[Rational]") {
	constexpr Rational b{ 5, 7 };
	constexpr auto a = 125;
	constexpr auto r = a / b;
	REQUIRE(float(a) / float(b) == Approx(r));
}


TEST_CASE("Rational add assign", "[Rational]") {
	constexpr Rational a{ 5, 7 };
	constexpr Rational b{ 125, 98 };
	auto r = a;
	r += b;
	REQUIRE(float(a) + float(b) == Approx(r));
}

TEST_CASE("Rational sub assign", "[Rational]") {
	constexpr Rational a{ 5, 7 };
	constexpr Rational b{ 125, 98 };
	auto r = a;
	r -= b;
	REQUIRE(float(a) - float(b) == Approx(r));
}

TEST_CASE("Rational mul assign", "[Rational]") {
	constexpr Rational a{ 5, 7 };
	constexpr Rational b{ 125, 98 };
	auto r = a;
	r *= b;
	REQUIRE(float(a) * float(b) == Approx(r));
}

TEST_CASE("Rational div assign", "[Rational]") {
	constexpr Rational a{ 5, 7 };
	constexpr Rational b{ 125, 98 };
	auto r = a;
	r /= b;
	REQUIRE(float(a) / float(b) == Approx(r));
}


TEST_CASE("Rational add/int assign", "[Rational]") {
	constexpr Rational a{ 5, 7 };
	constexpr auto b = 125;
	auto r = a;
	r += b;
	REQUIRE(float(a) + float(b) == Approx(r));
}

TEST_CASE("Rational sub/int assign", "[Rational]") {
	constexpr Rational a{ 5, 7 };
	constexpr auto b = 125;
	auto r = a;
	r -= b;
	REQUIRE(float(a) - float(b) == Approx(r));
}

TEST_CASE("Rational mul/int assign", "[Rational]") {
	constexpr Rational a{ 5, 7 };
	constexpr auto b = 125;
	auto r = a;
	r *= b;
	REQUIRE(float(a) * float(b) == Approx(r));
}

TEST_CASE("Rational div/int assign", "[Rational]") {
	constexpr Rational a{ 5, 7 };
	constexpr auto b = 125;
	auto r = a;
	r /= b;
	REQUIRE(float(a) / float(b) == Approx(r));
}


TEST_CASE("Rational random hammer", "[Rational]") {
	std::mt19937_64 rne(7354634985);
	std::uniform_int_distribution<int> num(-1000, 1000);
	std::uniform_int_distribution<int> den(1, 1000);
	std::uniform_int_distribution<int> op(0, 3);

	for (int i = 0; i < 200; ++i) {
		const Rational a = { num(rne), den(rne) };
		const Rational b = { num(rne), den(rne) };
		const auto opChoice = op(rne);
		const auto [result, expected] = [&] {
			switch (opChoice) {
				case 0: return std::make_pair(a + b, float(a) + float(b));
				case 1: return std::make_pair(a - b, float(a) - float(b));
				case 2: return std::make_pair(a * b, float(a) * float(b));
				case 3: return std::make_pair(a / b, float(a) / float(b));
			}
			std::terminate();
		}();
		REQUIRE(std::gcd(result.Numerator(), result.Denominator()) == 1);
		REQUIRE(float(result) == Approx(expected));
	}
}


TEST_CASE("Rational ==", "[Rational]") {
	REQUIRE(Rational{ 4, 6 } == Rational{ 6, 9 });
	REQUIRE(Rational{ -4, 6 } == Rational{ -6, 9 });
	REQUIRE(Rational{ 23, 9 } == Rational{ 23, 9 });
	REQUIRE_FALSE(Rational{ 23, 5 } == Rational{ 23, 9 });
	REQUIRE_FALSE(Rational{ -23, 5 } == Rational{ 23, 9 });
	REQUIRE_FALSE(Rational{ 21, 63 } == Rational{ 82, 26 });
}

TEST_CASE("Rational !=", "[Rational]") {
	REQUIRE_FALSE(Rational{ 4, 6 } != Rational{ 6, 9 });
	REQUIRE_FALSE(Rational{ -4, 6 } != Rational{ -6, 9 });
	REQUIRE_FALSE(Rational{ 23, 9 } != Rational{ 23, 9 });
	REQUIRE(Rational{ 23, 5 } != Rational{ 23, 9 });
	REQUIRE(Rational{ -23, 5 } != Rational{ 23, 9 });
	REQUIRE(Rational{ 21, 63 } != Rational{ 82, 26 });
}

TEST_CASE("Rational <", "[Rational]") {
	REQUIRE(Rational{ 314, 100 } < Rational{ 355, 113 });
	REQUIRE(Rational{ 27, 61 } < Rational{ 28, 62 });
	REQUIRE_FALSE(Rational{ 355, 113 } < Rational{ 314, 100 });
	REQUIRE_FALSE(Rational{ 28, 62 } < Rational{ 27, 61 });
}

TEST_CASE("Rational >", "[Rational]") {
	REQUIRE_FALSE(Rational{ 314, 100 } > Rational{ 355, 113 });
	REQUIRE_FALSE(Rational{ 27, 61 } > Rational{ 28, 62 });
	REQUIRE(Rational{ 355, 113 } > Rational{ 314, 100 });
	REQUIRE(Rational{ 28, 62 } > Rational{ 27, 61 });
}

TEST_CASE("Rational <=", "[Rational]") {
	REQUIRE(Rational{ 314, 100 } <= Rational{ 355, 113 });
	REQUIRE(Rational{ 27, 61 } <= Rational{ 28, 62 });
	REQUIRE(Rational{ 20, 60 } <= Rational{ 15, 45 });
	REQUIRE_FALSE(Rational{ 355, 113 } <= Rational{ 314, 100 });
	REQUIRE_FALSE(Rational{ 28, 62 } <= Rational{ 27, 61 });
}

TEST_CASE("Rational >=", "[Rational]") {
	REQUIRE_FALSE(Rational{ 314, 100 } >= Rational{ 355, 113 });
	REQUIRE_FALSE(Rational{ 27, 61 } >= Rational{ 28, 62 });
	REQUIRE(Rational{ 15, 45 } >= Rational{ 20, 60 });
	REQUIRE(Rational{ 355, 113 } >= Rational{ 314, 100 });
	REQUIRE(Rational{ 28, 62 } >= Rational{ 27, 61 });
}