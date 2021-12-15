#include "dspbb/Utility/TypeTraits.hpp"
#include <dspbb/ComputeKernels/VectorizedAlgorithms.hpp>

#include <array>
#include <vector>
#include <celero/Celero.h>
#include <iterator>
#include <numeric>
#include <random>

using namespace dspbb;


static constexpr std::array reductionSizes = {
	//2,
	//4,
	//6,
	//8,
	//10,
	//16,
	//32,
	//64,
	//96,
	//384,
	//1536,
	6144,
	//24576,
	//98304,
	//393216,
	//1572864,
	//6291456,
	//25165824,
	//100663296,
};

std::minstd_rand rne;
std::uniform_real_distribution<float> randomFloat(-1, 1);

template <class T>
class RandomArrayFixture : public celero::TestFixture {
public:
	std::vector<ExperimentValue> getExperimentValues() const override {
		std::vector<ExperimentValue> experimentValues;
		std::transform(reductionSizes.begin(), reductionSizes.end(), std::back_inserter(experimentValues), [](size_t size) {
			return ExperimentValue{ int64_t(size), std::max(int64_t(1), int64_t(10485760 / size)) };
		});
		return experimentValues;
	}

	virtual void setUp(const ExperimentValue& experimentValue) {
		size_t arraySize = experimentValue.Value;
		array = std::vector<T>(arraySize);
		std::array<T, 16> pattern;
		for (auto& v : pattern) {
			v = static_cast<T>(randomFloat(rne));
		}
		size_t index = 0;
		for (auto& v : array) {
			v = pattern[index];
			index = (index + 1) % pattern.size();
		}
	}

	std::vector<T> array;
};

template <class T>
T AccumulateSum(const T* u, size_t count) {
	return std::accumulate(u, u + count, static_cast<T>(0.0f));
}
template <class T>
T ReduceSum(const T* u, size_t count) {
	return std::reduce(u, u + count, static_cast<T>(0.0f));
}

BASELINE_F(Reduce_Float, std_accumulate, RandomArrayFixture<float>, 25, 500) {
	const auto result = AccumulateSum(array.data(), array.size());
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(Reduce_Float, std_reduce, RandomArrayFixture<float>, 25, 500) {
	const auto result = ReduceSum(array.data(), array.size());
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(Reduce_Float, dspbb_reduce, RandomArrayFixture<float>, 25, 500) {
	const auto result = kernels::Reduce(array.data(), array.size(), 0.0f, std::plus<>{});
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(Reduce_Float, dspbb_reduce_comp, RandomArrayFixture<float>, 25, 500) {
	const auto result = kernels::Reduce(array.data(), array.size(), 0.0f, kernels::plus_compensated<>{});
	celero::DoNotOptimizeAway(result);
}

BASELINE_F(Reduce_ComplexFloat, std_accumulate, RandomArrayFixture<std::complex<float>>, 25, 500) {
	const auto result = AccumulateSum(array.data(), array.size());
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(Reduce_ComplexFloat, std_reduce, RandomArrayFixture<std::complex<float>>, 25, 500) {
	const auto result = ReduceSum(array.data(), array.size());
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(Reduce_ComplexFloat, dspbb_reduce, RandomArrayFixture<std::complex<float>>, 25, 500) {
	const auto result = kernels::Reduce(array.data(), array.size(), std::complex<float>(0.0f), std::plus<>{});
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(Reduce_ComplexFloat, dspbb_reduce_comp, RandomArrayFixture<std::complex<float>>, 25, 500) {
	const auto result = kernels::Reduce(array.data(), array.size(), std::complex<float>(0.0f), kernels::plus_compensated<>{});
	celero::DoNotOptimizeAway(result);
}