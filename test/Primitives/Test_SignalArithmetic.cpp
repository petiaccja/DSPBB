#include "../TestUtils.hpp"

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

constexpr std::array<size_t, 2> sizes = { 1, 137 };


#define TEST_SIGNAL_BINARY_OPERATOR(NAME, OPERATOR)                                                              \
	TEMPLATE_PRODUCT_TEST_CASE("Signal binary " NAME, "[Signal Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) { \
		for (auto size : sizes) {                                                                                \
			SECTION(std::string("Size ") + std::to_string(size)) {                                               \
				using TestType0 = std::tuple_element_t<0, TestType>;                                             \
				using TestType1 = std::tuple_element_t<1, TestType>;                                             \
                                                                                                                 \
				const auto a = RandomPositiveSignal<TestType0>(size);                                            \
				const auto b = RandomPositiveSignal<TestType1>(size);                                            \
				const auto r0 = a OPERATOR b;                                                                    \
				const auto r1 = AsView(a) OPERATOR b;                                                            \
				const auto r2 = a OPERATOR AsView(b);                                                            \
				const auto r3 = AsView(a) OPERATOR AsView(b);                                                    \
                                                                                                                 \
				size_t numEqual0 = 0;                                                                            \
				size_t numEqual1 = 0;                                                                            \
				size_t numEqual2 = 0;                                                                            \
				size_t numEqual3 = 0;                                                                            \
				for (size_t i = 0; i < size; ++i) {                                                              \
					const auto expected = ApproxComplex(a[i] OPERATOR b[i]);                                     \
					numEqual0 += size_t(r0[i] == expected);                                                      \
					numEqual1 += size_t(r1[i] == expected);                                                      \
					numEqual2 += size_t(r2[i] == expected);                                                      \
					numEqual3 += size_t(r3[i] == expected);                                                      \
				}                                                                                                \
				REQUIRE(numEqual0 == size);                                                                      \
				REQUIRE(numEqual1 == size);                                                                      \
				REQUIRE(numEqual2 == size);                                                                      \
				REQUIRE(numEqual3 == size);                                                                      \
			}                                                                                                    \
		}                                                                                                        \
	}

#define TEST_SIGNAL_COMPOUND_OPERATOR(NAME, OPERATOR, VERIFY_OPERATOR)                                          \
	TEMPLATE_PRODUCT_TEST_CASE("Signal compound " NAME, "[Signal Arithmetic]", std::tuple, TYPES_BINARY_REAL) { \
		for (auto size : sizes) {                                                                               \
			SECTION(std::string("Size ") + std::to_string(size)) {                                              \
				using TestType0 = std::tuple_element_t<0, TestType>;                                            \
				using TestType1 = std::tuple_element_t<1, TestType>;                                            \
                                                                                                                \
				const auto a = RandomPositiveSignal<TestType0>(size);                                           \
				auto a0 = a;                                                                                    \
				auto a1 = a;                                                                                    \
				auto a2 = a;                                                                                    \
				auto a3 = a;                                                                                    \
				const auto b = RandomPositiveSignal<TestType1>(size);                                           \
				a0 OPERATOR b;                                                                                  \
				AsView(a1) OPERATOR b;                                                                          \
				a2 OPERATOR AsView(b);                                                                          \
				AsView(a3) OPERATOR AsView(b);                                                                  \
                                                                                                                \
				size_t numEqual0 = 0;                                                                           \
				size_t numEqual1 = 0;                                                                           \
				size_t numEqual2 = 0;                                                                           \
				size_t numEqual3 = 0;                                                                           \
				for (size_t i = 0; i < size; ++i) {                                                             \
					const auto expected = ApproxComplex(a[i] VERIFY_OPERATOR b[i]);                             \
					numEqual0 += size_t(a0[i] == expected);                                                     \
					numEqual1 += size_t(a1[i] == expected);                                                     \
					numEqual2 += size_t(a2[i] == expected);                                                     \
					numEqual3 += size_t(a3[i] == expected);                                                     \
				}                                                                                               \
				REQUIRE(numEqual0 == size);                                                                     \
				REQUIRE(numEqual1 == size);                                                                     \
				REQUIRE(numEqual2 == size);                                                                     \
				REQUIRE(numEqual3 == size);                                                                     \
			}                                                                                                   \
		}                                                                                                       \
	}


#define TEST_SIGNAL_BINARY_SCALAR_OPERATOR(NAME, OPERATOR)                                                              \
	TEMPLATE_PRODUCT_TEST_CASE("Signal binary scalar " NAME, "[Signal Arithmetic]", std::tuple, TYPES_BINARY_COMPLEX) { \
		for (auto size : sizes) {                                                                                       \
			SECTION(std::string("Size ") + std::to_string(size)) {                                                      \
				using TestType0 = std::tuple_element_t<0, TestType>;                                                    \
				using TestType1 = std::tuple_element_t<1, TestType>;                                                    \
                                                                                                                        \
				const auto a = RandomPositiveSignal<TestType0>(size);                                                   \
				const auto b = TestType1(1.55f);                                                                        \
				const auto r0 = a OPERATOR b;                                                                           \
				const auto r1 = AsView(a) OPERATOR b;                                                                   \
				const auto r2 = b OPERATOR a;                                                                           \
				const auto r3 = b OPERATOR AsView(a);                                                                   \
                                                                                                                        \
				size_t numEqual0 = 0;                                                                                   \
				size_t numEqual1 = 0;                                                                                   \
				size_t numEqual2 = 0;                                                                                   \
				size_t numEqual3 = 0;                                                                                   \
				for (size_t i = 0; i < size; ++i) {                                                                     \
					const auto expectedNormal = ApproxComplex(a[i] OPERATOR b);                                         \
					const auto expectedTranspose = ApproxComplex(b OPERATOR a[i]);                                      \
					numEqual0 += size_t(r0[i] == expectedNormal);                                                       \
					numEqual1 += size_t(r1[i] == expectedNormal);                                                       \
					numEqual2 += size_t(r2[i] == expectedTranspose);                                                    \
					numEqual3 += size_t(r3[i] == expectedTranspose);                                                    \
				}                                                                                                       \
				REQUIRE(numEqual0 == size);                                                                             \
				REQUIRE(numEqual1 == size);                                                                             \
				REQUIRE(numEqual2 == size);                                                                             \
				REQUIRE(numEqual3 == size);                                                                             \
			}                                                                                                           \
		}                                                                                                               \
	}


#define TEST_SIGNAL_COMPOUND_SCALAR_OPERATOR(NAME, OPERATOR, VERIFY_OPERATOR)                                          \
	TEMPLATE_PRODUCT_TEST_CASE("Signal compound scalar " NAME, "[Signal Arithmetic]", std::tuple, TYPES_BINARY_REAL) { \
		for (auto size : sizes) {                                                                                      \
			SECTION(std::string("Size ") + std::to_string(size)) {                                                     \
				using TestType0 = std::tuple_element_t<0, TestType>;                                                   \
				using TestType1 = std::tuple_element_t<1, TestType>;                                                   \
                                                                                                                       \
				const auto a = RandomPositiveSignal<TestType0>(size);                                                  \
				auto a0 = a;                                                                                           \
				auto a1 = a;                                                                                           \
				const auto b = TestType1(size);                                                                        \
				a0 OPERATOR b;                                                                                         \
				AsView(a1) OPERATOR b;                                                                                 \
                                                                                                                       \
				size_t numEqual0 = 0;                                                                                  \
				size_t numEqual1 = 0;                                                                                  \
				for (size_t i = 0; i < size; ++i) {                                                                    \
					const auto expected = ApproxComplex(a[i] VERIFY_OPERATOR b);                                       \
					numEqual0 += size_t(a0[i] == expected);                                                            \
					numEqual1 += size_t(a1[i] == expected);                                                            \
				}                                                                                                      \
				REQUIRE(numEqual0 == size);                                                                            \
				REQUIRE(numEqual1 == size);                                                                            \
			}                                                                                                          \
		}                                                                                                              \
	}

TEST_SIGNAL_BINARY_OPERATOR("multiply", *)
TEST_SIGNAL_BINARY_OPERATOR("divide", /)
TEST_SIGNAL_BINARY_OPERATOR("add", +)
TEST_SIGNAL_BINARY_OPERATOR("subtract", -)

TEST_SIGNAL_COMPOUND_OPERATOR("multiply", *=, *)
TEST_SIGNAL_COMPOUND_OPERATOR("divide", /=, /)
TEST_SIGNAL_COMPOUND_OPERATOR("add", +=, +)
TEST_SIGNAL_COMPOUND_OPERATOR("subtract", -=, -)

TEST_SIGNAL_BINARY_SCALAR_OPERATOR("multiply", *)
TEST_SIGNAL_BINARY_SCALAR_OPERATOR("divide", /)
TEST_SIGNAL_BINARY_SCALAR_OPERATOR("add", +)
TEST_SIGNAL_BINARY_SCALAR_OPERATOR("subtract", -)

TEST_SIGNAL_COMPOUND_SCALAR_OPERATOR("multiply", *=, *)
TEST_SIGNAL_COMPOUND_SCALAR_OPERATOR("divide", /=, /)
TEST_SIGNAL_COMPOUND_SCALAR_OPERATOR("add", +=, +)
TEST_SIGNAL_COMPOUND_SCALAR_OPERATOR("subtract", -=, -)