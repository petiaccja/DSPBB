#include <dspbb/ComputeKernels/VectorizedAlgorithms.hpp>
#include <dspbb/Utility/Numbers.hpp>

#include <catch2/catch.hpp>

using namespace dspbb;
using namespace std::complex_literals;

TEST_CASE("Reduce float", "[Kernels]") {
	std::array<float, 100> a;
	std::iota(a.begin(), a.end(), 1.0f);

	for (auto it = a.begin(); it != a.end(); ++it) {
		const auto reference = std::reduce(a.begin(), a.end(), 5.0f, std::plus<>{});
		const auto value = kernels::Reduce(a.begin(), a.end(), 5.0f, std::plus<>{});
		REQUIRE(reference == value);
	}
}

TEST_CASE("Reduce double", "[Kernels]") {
	std::array<double, 100> a;
	std::iota(a.begin(), a.end(), 1.0);

	for (auto it = a.begin(); it != a.end(); ++it) {
		const auto reference = std::reduce(a.begin(), a.end(), 5.0, std::plus<>{});
		const auto value = kernels::Reduce(a.begin(), a.end(), 5.0, std::plus<>{});
		REQUIRE(reference == value);
	}
}

TEST_CASE("Reduce complex", "[Kernels]") {
	std::array<std::complex<float>, 100> a;
	std::iota(a.begin(), a.end(), 1.0f);

	for (auto it = a.begin(); it != a.end(); ++it) {
		const auto reference = std::reduce(a.begin(), a.end(), 5.0f + 5.0if, std::plus<>{});
		const auto value = kernels::Reduce(a.begin(), a.end(), 5.0f + 5.0if, std::plus<>{});
		REQUIRE(reference == value);
	}
}

TEST_CASE("Reduce int", "[Kernels]") {
	std::array<int, 100> a;
	std::iota(a.begin(), a.end(), 1);

	for (auto it = a.begin(); it != a.end(); ++it) {
		const auto reference = std::reduce(a.begin(), a.end(), 5, std::plus<>{});
		const auto value = kernels::Reduce(a.begin(), a.end(), 5, std::plus<>{});
		REQUIRE(reference == value);
	}
}

TEST_CASE("Reduce compensated", "[Kernels]") {
	std::array<float, 100> a;
	std::iota(a.begin(), a.end(), 1.0f);
	const auto reference = std::reduce(a.begin(), a.end(), 5.0f, std::plus<>{});
	const auto value = kernels::Reduce(a.begin(), a.end(), 5.0f, kernels::plus_compensated<>{});
	REQUIRE(reference == value);
}

TEST_CASE("Reduce compensation effects", "[Kernels]") {
	constexpr size_t count = 1 << 18;
	constexpr float item = 1 + 3.814697265625e-6f;
	std::vector<float> a(count, item);
	const float sumRegular = kernels::Reduce(a.begin(), a.end(), 0.0f, std::plus<>{});
	const float sumCompensated = kernels::Reduce(a.begin(), a.end(), 0.0f, dspbb::kernels::plus_compensated<>{});
	const float expected = item * float(count);
	REQUIRE(sumCompensated == expected);
	REQUIRE(sumRegular < expected);
}


TEST_CASE("MapReduce", "[Kernels]") {
	const auto reduceOp = [](const auto& a, const auto& b) { return a + b; };
	const auto mapOp = [](const auto& a) { return 1.0 / (a * a); };

	std::vector<double> a(50000);
	std::iota(a.begin(), a.end(), 1);
	const auto sum = kernels::MapReduce(a.data(), a.size(), 10.0, reduceOp, mapOp);
	REQUIRE(std::sqrt((sum - 10.0) * 6) == Approx(pi_v<double>).margin(0.001));
}

TEST_CASE("MapReduceVectorized", "[Kernels]") {
	const auto reduceOp = [](const auto& a, const auto& b) { return a + b; };
	const auto mapOp = [](const auto& a) { return 1.0 / (a * a); };

	std::vector<double> a(50000);
	std::iota(a.begin(), a.end(), 1);
	const auto sum = kernels::MapReduceVectorized(a.data(), a.size(), 10.0, reduceOp, mapOp);
	REQUIRE(std::sqrt((sum - 10.0) * 6) == Approx(pi_v<double>).margin(0.001));
}

TEST_CASE("InnerProduct", "[Kernels]") {
	const auto reduceOp = [](const auto& a, const auto& b) { return a + b; };
	const auto productOp = [](const auto& a, const auto& b) { return 1.0 / (a * b); };

	std::array<double, 50000> a;
	std::iota(a.begin(), a.end(), 1);
	const auto sum = kernels::InnerProduct(a.data(), a.data(), a.size(), 10.0, productOp, reduceOp);
	REQUIRE(std::sqrt((sum - 10.0) * 6) == Approx(pi_v<double>).margin(0.001));
}

TEST_CASE("InnerProductVectorized", "[Kernels]") {
	const auto reduceOp = [](const auto& a, const auto& b) { return a + b; };
	const auto productOp = [](const auto& a, const auto& b) { return 1.0 / (a * b); };

	std::array<double, 50000> a;
	std::iota(a.begin(), a.end(), 1);
	const auto sum = kernels::InnerProductVectorized(a.data(), a.data(), a.size(), 10.0, productOp, reduceOp);
	REQUIRE(std::sqrt((sum - 10.0) * 6) == Approx(pi_v<double>).margin(0.001));
}