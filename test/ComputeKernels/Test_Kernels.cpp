#include <dspbb/ComputeKernels/VectorizedAlgorithms.hpp>
#include <dspbb/Utility/Numbers.hpp>

#include <catch2/catch.hpp>

using namespace dspbb;

TEST_CASE("Reduce large", "[Kernels]") {
	std::array<double, 100> a;
	std::iota(a.begin(), a.end(), 1);
	const auto sum = kernels::Reduce(a.data(), a.size(), 1000.0, [](const auto& a, const auto& b) { return a + b; });
	REQUIRE(sum == Approx(6050.0));
}

TEST_CASE("Reduce small", "[Kernels]") {
	std::array<double, 7> a;
	std::iota(a.begin(), a.end(), 1);
	const auto prod = kernels::Reduce(a.data(), a.size(), 8.0, [](const auto& a, const auto& b) { return a * b; });
	REQUIRE(prod == Approx(40320.0));
}

TEST_CASE("Reduce compensated", "[Kernels]") {
	std::array<double, 100> a;
	std::iota(a.begin(), a.end(), 1);
	const auto sum = kernels::Reduce(a.data(), a.size(), 1000.0, dspbb::kernels::plus_compensated<>{});
	REQUIRE(sum == Approx(6050.0));
}

TEST_CASE("Reduce compensation effects", "[Kernels]") {
	constexpr size_t count = 1 << 18;
	constexpr float item = 1 + 3.814697265625e-6f;
	std::vector<float> a(count, item);
	const float sumRegular = kernels::Reduce(a.data(), a.size(), 0.0f, std::plus<>{});
	const float sumCompensated = kernels::Reduce(a.data(), a.size(), 0.0f, dspbb::kernels::plus_compensated<>{});
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