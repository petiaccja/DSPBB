#include <catch2/catch.hpp>
#include <dspbb/Utility/Numbers.hpp>
#include <dspbb/Vectorization/Kernels.hpp>

using namespace dspbb;

TEST_CASE("Reduce", "[Kernels]") {
	std::array<float, 100> a;
	std::iota(a.begin(), a.end(), 1);
	const auto sum = Reduce(a.data(), a.size(), 1000.f, [](const auto& a, const auto& b) { return a + b; });
	REQUIRE(sum == Approx(6050.f));
}

TEST_CASE("ReduceVectorized", "[Kernels]") {
	std::array<double, 7> a;
	std::iota(a.begin(), a.end(), 1);
	const auto prod = ReduceVectorized(a.data(), a.size(), 8.0, [](const auto& a, const auto& b) { return a * b; });
	REQUIRE(prod == Approx(40320.0));
}

TEST_CASE("MapReduce", "[Kernels]") {
	const auto reduceOp = [](const auto& a, const auto& b) { return a + b; };
	const auto mapOp = [](const auto& a) { return 1.0 / (a * a); };

	std::array<double, 50000> a;
	std::iota(a.begin(), a.end(), 1);
	const auto sum = MapReduce(a.data(), a.size(), 10.0, reduceOp, mapOp);
	REQUIRE(std::sqrt((sum - 10.0) * 6) == Approx(pi_v<double>).margin(0.001));
}

TEST_CASE("MapReduceVectorized", "[Kernels]") {
	const auto reduceOp = [](const auto& a, const auto& b) { return a + b; };
	const auto mapOp = [](const auto& a) { return 1.0 / (a * a); };

	std::array<double, 50000> a;
	std::iota(a.begin(), a.end(), 1);
	const auto sum = MapReduceVectorized(a.data(), a.size(), 10.0, reduceOp, mapOp);
	REQUIRE(std::sqrt((sum - 10.0) * 6) == Approx(pi_v<double>).margin(0.001));
}
