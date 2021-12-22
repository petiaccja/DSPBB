#include <dspbb/Kernels/Numeric.hpp>
#include <dspbb/Utility/Numbers.hpp>

#include <catch2/catch.hpp>

using namespace dspbb;
using namespace std::complex_literals;


//------------------------------------------------------------------------------
// Reduce
//------------------------------------------------------------------------------

TEST_CASE("Reduce float", "[Kernels - Numeric]") {
	std::array<float, 100> a;
	std::iota(a.begin(), a.end(), 1.0f);

	for (auto it = a.begin(); it != a.end(); ++it) {
		const auto reference = std::reduce(a.begin(), it, 5.0f, std::plus<>{});
		const auto value = kernels::Reduce(a.begin(), it, 5.0f, std::plus<>{});
		REQUIRE(reference == value);
	}
}

TEST_CASE("Reduce double", "[Kernels - Numeric]") {
	std::array<double, 100> a;
	std::iota(a.begin(), a.end(), 1.0);

	for (auto it = a.begin(); it != a.end(); ++it) {
		const auto reference = std::reduce(a.begin(), it, 5.0, std::plus<>{});
		const auto value = kernels::Reduce(a.begin(), it, 5.0, std::plus<>{});
		REQUIRE(reference == value);
	}
}

TEST_CASE("Reduce complex", "[Kernels - Numeric]") {
	std::array<std::complex<float>, 100> a;
	std::iota(a.begin(), a.end(), 1.0f);

	for (auto it = a.begin(); it != a.end(); ++it) {
		const auto reference = std::reduce(a.begin(), it, 5.0f + 5.0if, std::plus<>{});
		const auto value = kernels::Reduce(a.begin(), it, 5.0f + 5.0if, std::plus<>{});
		REQUIRE(reference == value);
	}
}

TEST_CASE("Reduce int", "[Kernels - Numeric]") {
	std::array<int, 100> a;
	std::iota(a.begin(), a.end(), 1);

	for (auto it = a.begin(); it != a.end(); ++it) {
		const auto reference = std::reduce(a.begin(), it, 5, std::plus<>{});
		const auto value = kernels::Reduce(a.begin(), it, 5, std::plus<>{});
		REQUIRE(reference == value);
	}
}

TEST_CASE("Reduce compensated", "[Kernels - Numeric]") {
	std::array<float, 100> a;
	std::iota(a.begin(), a.end(), 1.0f);
	const auto reference = std::reduce(a.begin(), a.end(), 5.0f, std::plus<>{});
	const auto value = kernels::Reduce(a.begin(), a.end(), 5.0f, dspbb::plus_compensated<>{});
	REQUIRE(reference == value);
}

TEST_CASE("Reduce compensation effects", "[Kernels - Numeric]") {
	constexpr size_t count = 1 << 18;
	constexpr float item = 1 + 3.814697265625e-6f;
	std::vector<float> a(count, item);
	const float sumRegular = kernels::Reduce(a.begin(), a.end(), 0.0f, std::plus<>{});
	const float sumCompensated = kernels::Reduce(a.begin(), a.end(), 0.0f, dspbb::plus_compensated<>{});
	const float expected = item * float(count);
	REQUIRE(sumCompensated == expected);
	REQUIRE(sumRegular < expected);
}



//------------------------------------------------------------------------------
// Transform reduce
//------------------------------------------------------------------------------

TEST_CASE("Transform reduce float", "[Kernels - Numeric]") {
	std::array<float, 100> a;
	std::iota(a.begin(), a.end(), 1.0f);

	for (auto it = a.begin(); it != a.end(); ++it) {
		const auto reference = std::transform_reduce(a.begin(), it, 5.0f, std::plus<>{}, [](const auto& arg) { return arg * arg; });
		const auto value = kernels::TransformReduce(a.begin(), it, 5.0f, std::plus<>{}, [](const auto& arg) { return arg * arg; });
		REQUIRE(reference == value);
	}
}

TEST_CASE("Transform reduce double", "[Kernels - Numeric]") {
	std::array<double, 100> a;
	std::iota(a.begin(), a.end(), 1.0);

	for (auto it = a.begin(); it != a.end(); ++it) {
		const auto reference = std::transform_reduce(a.begin(), it, 5.0, std::plus<>{}, [](const auto& arg) { return arg * arg; });
		const auto value = kernels::TransformReduce(a.begin(), it, 5.0, std::plus<>{}, [](const auto& arg) { return arg * arg; });
		REQUIRE(reference == value);
	}
}

TEST_CASE("Transform reduce complex", "[Kernels - Numeric]") {
	std::array<std::complex<float>, 100> a;
	std::iota(a.begin(), a.end(), 1.0f);

	for (auto it = a.begin(); it != a.end(); ++it) {
		const auto reference = std::transform_reduce(a.begin(), it, 5.0f + 5.0if, std::plus<>{}, [](const auto& arg) { return arg * arg; });
		const auto value = kernels::TransformReduce(a.begin(), it, 5.0f + 5.0if, std::plus<>{}, [](const auto& arg) { return arg * arg; });
		REQUIRE(reference == value);
	}
}

TEST_CASE("Transform reduce int", "[Kernels - Numeric]") {
	std::array<int, 100> a;
	std::iota(a.begin(), a.end(), 1);

	for (auto it = a.begin(); it != a.end(); ++it) {
		const auto reference = std::transform_reduce(a.begin(), it, 5, std::plus<>{}, [](const auto& arg) { return arg * arg; });
		const auto value = kernels::TransformReduce(a.begin(), it, 5, std::plus<>{}, [](const auto& arg) { return arg * arg; });
		REQUIRE(reference == value);
	}
}


//------------------------------------------------------------------------------
// Inner product
//------------------------------------------------------------------------------

TEST_CASE("InnerProduct float", "[Kernels - Numeric]") {
	std::array<float, 100> a;
	std::array<float, 100> b;
	std::iota(a.begin(), a.end(), 1.0f);
	std::iota(b.begin(), b.end(), 3.0f);

	for (auto it = a.begin(); it != a.end(); ++it) {
		INFO(it - a.begin());
		const auto reference = std::inner_product(a.begin(), it, b.begin(), 5.0f, std::plus<>{}, std::multiplies<>{});
		const auto value = kernels::InnerProduct(a.begin(), it, b.begin(), 5.0f, std::plus<>{}, std::multiplies<>{});
		REQUIRE(reference == value);
	}
}

TEST_CASE("InnerProduct double", "[Kernels - Numeric]") {
	std::array<double, 100> a;
	std::array<double, 100> b;
	std::iota(a.begin(), a.end(), 1.0);
	std::iota(b.begin(), b.end(), 3.0);

	for (auto it = a.begin(); it != a.end(); ++it) {
		INFO(it - a.begin());
		const auto reference = std::inner_product(a.begin(), it, b.begin(), 5.0, std::plus<>{}, std::multiplies<>{});
		const auto value = kernels::InnerProduct(a.begin(), it, b.begin(), 5.0, std::plus<>{}, std::multiplies<>{});
		REQUIRE(reference == value);
	}
}

TEST_CASE("InnerProduct complex", "[Kernels - Numeric]") {
	std::array<std::complex<float>, 100> a;
	std::array<std::complex<float>, 100> b;
	std::iota(a.begin(), a.end(), 1.0f);
	std::iota(b.begin(), b.end(), 3.0f);

	for (auto it = a.begin(); it != a.end(); ++it) {
		const auto reference = std::inner_product(a.begin(), it, b.begin(), std::complex<float>(5.0f), std::plus<>{}, std::multiplies<>{});
		const auto value = kernels::InnerProduct(a.begin(), it, b.begin(), std::complex<float>(5.0f), std::plus<>{}, std::multiplies<>{});
		REQUIRE(reference == value);
	}
}

TEST_CASE("InnerProduct int", "[Kernels - Numeric]") {
	std::array<int, 100> a;
	std::array<int, 100> b;
	std::iota(a.begin(), a.end(), 1);
	std::iota(b.begin(), b.end(), 3);

	for (auto it = a.begin(); it != a.end(); ++it) {
		const auto reference = std::inner_product(a.begin(), it, b.begin(), 5, std::plus<>{}, std::multiplies<>{});
		const auto value = kernels::InnerProduct(a.begin(), it, b.begin(), 5, std::plus<>{}, std::multiplies<>{});
		REQUIRE(reference == value);
	}
}

//------------------------------------------------------------------------------
// Transform
//------------------------------------------------------------------------------

TEST_CASE("Transform unary float", "[Kernels - Numeric]") {
	std::array<float, 100> a;
	std::iota(a.begin(), a.end(), 1.0f);

	std::array<float, 100> reference;
	std::array<float, 100> value;

	std::transform(a.begin(), a.end(), reference.begin(), [](const auto& v) { return -v; });
	const auto endIt = kernels::Transform(a.begin(), a.end(), value.begin(), [](const auto& v) { return -v; });
	REQUIRE(endIt == value.end());
	REQUIRE(reference == value);
}

TEST_CASE("Transform unary int", "[Kernels - Numeric]") {
	std::array<int, 100> a;
	std::iota(a.begin(), a.end(), 1);

	std::array<int, 100> reference;
	std::array<int, 100> value;

	std::transform(a.begin(), a.end(), reference.begin(), [](const auto& v) { return -v; });
	const auto endIt = kernels::Transform(a.begin(), a.end(), value.begin(), [](const auto& v) { return -v; });
	REQUIRE(endIt == value.end());
	REQUIRE(reference == value);
}


TEST_CASE("Transform binary float", "[Kernels - Numeric]") {
	std::array<float, 100> a;
	std::array<float, 100> b;
	std::iota(a.begin(), a.end(), 1.0f);
	std::iota(b.begin(), b.end(), 3.0f);

	std::array<float, 100> reference;
	std::array<float, 100> value;

	std::transform(a.begin(), a.end(), b.begin(), reference.begin(), std::multiplies<>{});
	const auto endIt = kernels::Transform(a.begin(), a.end(), b.begin(), value.begin(), std::multiplies<>{});
	REQUIRE(endIt == value.end());
	REQUIRE(reference == value);
}

TEST_CASE("Transform binary int", "[Kernels - Numeric]") {
	std::array<int, 100> a;
	std::array<int, 100> b;
	std::iota(a.begin(), a.end(), 1);
	std::iota(b.begin(), b.end(), 3);

	std::array<int, 100> reference;
	std::array<int, 100> value;

	std::transform(a.begin(), a.end(), b.begin(), reference.begin(), std::multiplies<>{});
	const auto endIt = kernels::Transform(a.begin(), a.end(), b.begin(), value.begin(), std::multiplies<>{});
	REQUIRE(endIt == value.end());
	REQUIRE(reference == value);
}

TEST_CASE("Transform unary self-assign", "[Kernels - Numeric]") {
	std::array<float, 100> a;
	std::iota(a.begin(), a.end(), 1.0f);

	std::array<float, 100> reference;

	std::transform(a.begin(), a.end(), reference.begin(), [](const auto& v) { return -v; });
	const auto endIt = kernels::Transform(a.begin(), a.end(), a.begin(), [](const auto& v) { return -v; });
	REQUIRE(endIt == a.end());
	REQUIRE(reference == a);
}

TEST_CASE("Transform binary self-assign", "[Kernels - Numeric]") {
	std::array<float, 100> a;
	std::array<float, 100> b;
	std::iota(a.begin(), a.end(), 1.0f);
	std::iota(b.begin(), b.end(), 3.0f);

	std::array<float, 100> reference;

	std::transform(a.begin(), a.end(), b.begin(), reference.begin(), std::multiplies<>{});
	const auto endIt = kernels::Transform(a.begin(), a.end(), b.begin(), a.begin(), std::multiplies<>{});
	REQUIRE(endIt == a.end());
	REQUIRE(reference == a);
}
