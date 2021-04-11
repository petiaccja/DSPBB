#include <catch2/catch.hpp>
#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/Primitives/SignalArithmetic.hpp>
#include <dspbb/Primitives/SignalView.hpp>

using namespace dspbb;

static_assert(is_mutable_signal_v<Signal<float, TIME_DOMAIN>>, "");
static_assert(is_mutable_signal_v<Signal<float, TIME_DOMAIN>&>, "");
static_assert(is_mutable_signal_v<Signal<float, TIME_DOMAIN>&&>, "");
static_assert(!is_mutable_signal_v<const Signal<float, TIME_DOMAIN>>, "");
static_assert(!is_mutable_signal_v<const Signal<float, TIME_DOMAIN>&>, "");
static_assert(!is_mutable_signal_v<const Signal<float, TIME_DOMAIN>&&>, "");

static_assert(is_mutable_signal_v<SignalView<float, TIME_DOMAIN>>, "");
static_assert(is_mutable_signal_v<SignalView<float, TIME_DOMAIN>&>, "");
static_assert(is_mutable_signal_v<SignalView<float, TIME_DOMAIN>&&>, "");
static_assert(is_mutable_signal_v<const SignalView<float, TIME_DOMAIN>>, "");
static_assert(is_mutable_signal_v<const SignalView<float, TIME_DOMAIN>&>, "");
static_assert(is_mutable_signal_v<const SignalView<float, TIME_DOMAIN>&&>, "");
static_assert(!is_mutable_signal_v<SignalView<const float, TIME_DOMAIN>>, "");
static_assert(!is_mutable_signal_v<SignalView<const float, TIME_DOMAIN>&>, "");
static_assert(!is_mutable_signal_v<SignalView<const float, TIME_DOMAIN>&&>, "");
static_assert(!is_mutable_signal_v<const SignalView<const float, TIME_DOMAIN>>, "");
static_assert(!is_mutable_signal_v<const SignalView<const float, TIME_DOMAIN>&>, "");
static_assert(!is_mutable_signal_v<const SignalView<const float, TIME_DOMAIN>&&>, "");

constexpr float a0 = 3.5f, a1 = 2.9f;
constexpr float b0 = 9.3f, b1 = 2.5f;
constexpr float bs = 2.63f;
constexpr float as = 1.75f;

TEST_CASE("Multiply", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	TimeSignal<float> b = { b0, b1 };
	auto r = a * b;
	REQUIRE(r[0] == Approx(a0 * b0));
	REQUIRE(r[1] == Approx(a1 * b1));
}

TEST_CASE("Divide", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	TimeSignal<float> b = { b0, b1 };
	auto r = a / b;
	;
	REQUIRE(r[1] == Approx(a1 / b1));
	REQUIRE(r[0] == Approx(a0 / b0));
}

TEST_CASE("Add", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	TimeSignal<float> b = { b0, b1 };
	auto r = a + b;
	auto rv = AsView(a) * b;
	REQUIRE(r[0] == Approx(a0 + b0));
	REQUIRE(r[1] == Approx(a1 + b1));
}

TEST_CASE("Subtract", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	TimeSignal<float> b = { b0, b1 };
	auto r = a - b;
	REQUIRE(r[0] == Approx(a0 - b0));
	REQUIRE(r[1] == Approx(a1 - b1));
}

TEST_CASE("Multiply mix", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	TimeSignal<float> b = { b0, b1 };
	auto r1 = AsView(a) * b;
	auto r2 = AsView(a) * AsView(b);
	auto r3 = a * AsView(b);
	auto r4 = AsView(a) * AsConstView(b);
	REQUIRE(r1[0] == r2[0]);
	REQUIRE(r1[0] == r3[0]);
	REQUIRE(r1[0] == r4[0]);
	REQUIRE(r1[1] == r2[1]);
	REQUIRE(r1[1] == r3[1]);
	REQUIRE(r1[1] == r4[1]);
}

TEST_CASE("Compound multiply", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	TimeSignal<double> b = { b0, b1 };
	a *= b;
	REQUIRE(a[0] == Approx(a0 * b0));
	REQUIRE(a[1] == Approx(a1 * b1));
}

TEST_CASE("Compound divide", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	TimeSignal<double> b = { b0, b1 };
	a /= b;
	REQUIRE(a[0] == Approx(a0 / b0));
	REQUIRE(a[1] == Approx(a1 / b1));
}

TEST_CASE("Compound add", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	TimeSignal<double> b = { b0, b1 };
	a += b;
	REQUIRE(a[0] == Approx(a0 + b0));
	REQUIRE(a[1] == Approx(a1 + b1));
}

TEST_CASE("Compound subtract", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	TimeSignal<double> b = { b0, b1 };
	a -= b;
	REQUIRE(a[0] == Approx(a0 - b0));
	REQUIRE(a[1] == Approx(a1 - b1));
}


TEST_CASE("Scalar multiply", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	auto r = a * bs;
	REQUIRE(r[0] == Approx(a0 * bs));
	REQUIRE(r[1] == Approx(a1 * bs));
}

TEST_CASE("Scalar divide", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	auto r = a / bs;
	REQUIRE(r[0] == Approx(a0 / bs));
	REQUIRE(r[1] == Approx(a1 / bs));
}

TEST_CASE("Scalar add", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	auto r = a + 2;
	REQUIRE(r[0] == Approx(a0 + 2));
	REQUIRE(r[1] == Approx(a1 + 2));
}

TEST_CASE("Scalar substract", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	auto r = a - bs;
	REQUIRE(r[0] == Approx(a0 - bs));
	REQUIRE(r[1] == Approx(a1 - bs));
}


TEST_CASE("Scalar reverse multiply", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	auto r = bs * a;
	REQUIRE(r[0] == Approx(bs * a0));
	REQUIRE(r[1] == Approx(bs * a1));
}

TEST_CASE("Scalar reverse divide", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	auto r = bs / a;
	REQUIRE(r[0] == Approx(bs / a0));
	REQUIRE(r[1] == Approx(bs / a1));
}

TEST_CASE("Scalar reverse add", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	auto r = 2 + a;
	REQUIRE(r[0] == Approx(2 + a0));
	REQUIRE(r[1] == Approx(2 + a1));
}

TEST_CASE("Scalar reverse substract", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	auto r = bs - a;
	REQUIRE(r[0] == Approx(bs - a0));
	REQUIRE(r[1] == Approx(bs - a1));
}


TEST_CASE("Scalar compound multiply", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	a *= bs;
	REQUIRE(a[0] == Approx(a0 * bs));
	REQUIRE(a[1] == Approx(a1 * bs));
}

TEST_CASE("Scalar compound divide", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	a /= bs;
	REQUIRE(a[0] == Approx(a0 / bs));
	REQUIRE(a[1] == Approx(a1 / bs));
}

TEST_CASE("Scalar compound add", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	a += 2;
	REQUIRE(a[0] == Approx(a0 + 2));
	REQUIRE(a[1] == Approx(a1 + 2));
}

TEST_CASE("Scalar compound substract", "[Signal Arithmetic]") {
	TimeSignal<float> a = { a0, a1 };
	a -= bs;
	REQUIRE(a[0] == Approx(a0 - bs));
	REQUIRE(a[1] == Approx(a1 - bs));
}