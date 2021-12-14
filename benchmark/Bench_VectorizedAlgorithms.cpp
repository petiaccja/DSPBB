#include "dspbb/Utility/TypeTraits.hpp"
#include <dspbb/ComputeKernels/VectorizedAlgorithms.hpp>

#include <array>
#include <celero/Celero.h>
#include <iostream>
#include <iterator>
#include <memory_resource>
#include <random>
#include <stack>

using namespace dspbb;

float LinearSum(const float* u, size_t count) {
	return std::accumulate(u, u + count, 0.0f);
}

std::minstd_rand rne;
std::uniform_real_distribution<float> randomFloat(-1, 1);

class RandomArrayFicture : public celero::TestFixture {
public:
	static constexpr std::array sizes = { 32, 256, 1024, 4096, 16384, 131072, 1048576, 8388608, 67108864 };

	std::vector<ExperimentValue> getExperimentValues() const override {
		std::vector<ExperimentValue> experimentValues;
		std::transform(sizes.begin(), sizes.end(), std::back_inserter(experimentValues), [](size_t size) {
			return ExperimentValue{ int64_t(size), std::max(int64_t(1), int64_t(8388608 / size)) };
		});
		return experimentValues;
	}

	virtual void setUp(const ExperimentValue& experimentValue) {
		size_t arraySize = experimentValue.Value;
		array = std::vector<float>(arraySize);
		std::array<float, 16> pattern;
		for (auto& v : pattern) {
			v = randomFloat(rne);
		}
		size_t index = 0;
		for (auto& v : array) {
			v = pattern[index];
			index = (index + 1) % pattern.size();
		}
	}

	std::vector<float> array;
};

BASELINE_F(Reduce, Linear, RandomArrayFicture, 25, 500) {
	const float result = LinearSum(array.data(), array.size());
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(Reduce, Vectorized, RandomArrayFicture, 25, 500) {
	const float result = kernels::ReduceVectorized(array.data(), array.size(), 0.0f, [](const auto& acc, const auto& x) { return acc + x; });
	celero::DoNotOptimizeAway(result);
}