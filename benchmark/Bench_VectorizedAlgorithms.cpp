#include "dspbb/Utility/TypeTraits.hpp"
#include <dspbb/Kernels/Math.hpp>
#include <dspbb/Kernels/Numeric.hpp>

#include <array>
#include <celero/Celero.h>
#include <iterator>
#include <numeric>
#include <random>
#include <vector>

using namespace dspbb;


//------------------------------------------------------------------------------
// Input sizes for which to benchmark
//------------------------------------------------------------------------------

static constexpr std::array reductionSizes = {
	2,
	4,
	6,
	8,
	10,
	12,
	14,
	16,
	25,
	32,
	36,
	42,
	48,
	56,
	64,
	96,
	384,
	1536,
	6144,
	24576,
	98304,
	393216,
	1572864,
	6291456,
	25165824,
	100663296,
};

//------------------------------------------------------------------------------
// Fixtures to generate random input
//------------------------------------------------------------------------------

static std::minstd_rand rne;
static std::uniform_real_distribution<float> randomFloat(-1, 1);

template <class T>
class RandomArrayFixture : public celero::TestFixture {
public:
	std::vector<ExperimentValue> getExperimentValues() const override {
		std::vector<ExperimentValue> experimentValues;
		std::transform(reductionSizes.begin(), reductionSizes.end(), std::back_inserter(experimentValues), [](size_t size) {
			return ExperimentValue{ int64_t(size), std::max(int64_t(1), int64_t(5 * 1048576 / size)) };
		});
		return experimentValues;
	}

	void setUp(const ExperimentValue& experimentValue) override {
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
class RandomArray2Fixture : public celero::TestFixture {
public:
	std::vector<ExperimentValue> getExperimentValues() const override {
		std::vector<ExperimentValue> experimentValues;
		std::transform(reductionSizes.begin(), reductionSizes.end(), std::back_inserter(experimentValues), [](size_t size) {
			return ExperimentValue{ int64_t(size), std::max(int64_t(1), int64_t(5 * 1048576 / size)) };
		});
		return experimentValues;
	}

	void setUp(const ExperimentValue& experimentValue) override {
		size_t arraySize = experimentValue.Value;
		array1 = std::vector<T>(arraySize);
		array2 = std::vector<T>(arraySize);
		std::array<T, 16> pattern;
		for (auto& v : pattern) {
			v = static_cast<T>(randomFloat(rne));
		}
		size_t index = 0;
		for (auto& v : array1) {
			v = pattern[index];
			index = (index + 1) % pattern.size();
		}
		for (auto& v : array2) {
			v = pattern[index];
			index = (index + 1) % pattern.size();
		}
	}

	std::vector<T> array1;
	std::vector<T> array2;
};

template <class T>
class RandomArray3Fixture : public celero::TestFixture {
public:
	std::vector<ExperimentValue> getExperimentValues() const override {
		std::vector<ExperimentValue> experimentValues;
		std::transform(reductionSizes.begin(), reductionSizes.end(), std::back_inserter(experimentValues), [](size_t size) {
			return ExperimentValue{ int64_t(size), std::max(int64_t(1), int64_t(5 * 1048576 / size)) };
		});
		return experimentValues;
	}

	void setUp(const ExperimentValue& experimentValue) override {
		size_t arraySize = experimentValue.Value;
		array1 = std::vector<T>(arraySize);
		array2 = std::vector<T>(arraySize);
		array3 = std::vector<T>(arraySize, T(0));
		std::array<T, 16> pattern;
		for (auto& v : pattern) {
			v = static_cast<T>(randomFloat(rne));
		}
		size_t index = 0;
		for (auto& v : array1) {
			v = pattern[index];
			index = (index + 1) % pattern.size();
		}
		for (auto& v : array2) {
			v = pattern[index];
			index = (index + 1) % pattern.size();
		}
	}

	std::vector<T> array1;
	std::vector<T> array2;
	std::vector<T> array3;
};


//------------------------------------------------------------------------------
// Reduce
//------------------------------------------------------------------------------

BASELINE_F(Reduce_Float, std_accumulate, RandomArrayFixture<float>, 25, 500) {
	const auto result = std::accumulate(array.begin(), array.end(), 0.0f, std::plus<>{});
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(Reduce_Float, std_reduce, RandomArrayFixture<float>, 25, 500) {
	const auto result = std::reduce(array.begin(), array.end(), 0.0f, std::plus<>{});
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(Reduce_Float, dspbb_reduce, RandomArrayFixture<float>, 25, 500) {
	const auto result = kernels::Reduce(array.begin(), array.end(), 0.0f, std::plus<>{});
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(Reduce_Float, dspbb_reduce_comp, RandomArrayFixture<float>, 25, 500) {
	const auto result = kernels::Reduce(array.begin(), array.end(), 0.0f, dspbb::plus_compensated<>{});
	celero::DoNotOptimizeAway(result);
}


//------------------------------------------------------------------------------
// Reduce std::complex
//------------------------------------------------------------------------------

BASELINE_F(Reduce_ComplexFloat, std_accumulate, RandomArrayFixture<std::complex<float>>, 25, 500) {
	const auto result = std::accumulate(array.begin(), array.end(), std::complex<float>(0.0f), std::plus<>{});
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(Reduce_ComplexFloat, std_reduce, RandomArrayFixture<std::complex<float>>, 25, 500) {
	const auto result = std::reduce(array.begin(), array.end(), std::complex<float>(0.0f), std::plus<>{});
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(Reduce_ComplexFloat, dspbb_reduce, RandomArrayFixture<std::complex<float>>, 25, 500) {
	const auto result = kernels::Reduce(array.begin(), array.end(), std::complex<float>(0.0f), std::plus<>{});
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(Reduce_ComplexFloat, dspbb_reduce_comp, RandomArrayFixture<std::complex<float>>, 25, 500) {
	const auto result = kernels::Reduce(array.begin(), array.end(), std::complex<float>(0.0f), dspbb::plus_compensated<>{});
	celero::DoNotOptimizeAway(result);
}


//------------------------------------------------------------------------------
// Transform reduce
//------------------------------------------------------------------------------

BASELINE_F(TransformReduce, std, RandomArrayFixture<float>, 25, 500) {
	const auto result = std::transform_reduce(array.begin(), array.end(), 0.0f, std::plus<>{}, [](const auto& v) { return v * v; });
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(TransformReduce, dspbb, RandomArrayFixture<float>, 25, 500) {
	const auto result = kernels::TransformReduce(array.begin(), array.end(), 0.0f, std::plus<>{}, [](const auto& v) { return v * v; });
	celero::DoNotOptimizeAway(result);
}


//------------------------------------------------------------------------------
// Inner product
//------------------------------------------------------------------------------

BASELINE_F(InnerProduct, std, RandomArray2Fixture<float>, 25, 500) {
	const auto result = std::inner_product(array1.begin(), array1.end(), array2.begin(), 0.0f, std::plus<>{}, std::multiplies<>{});
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(InnerProduct, dspbb, RandomArray2Fixture<float>, 25, 500) {
	const auto result = kernels::InnerProduct(array1.begin(), array1.end(), array2.begin(), 0.0f, std::plus<>{}, std::multiplies<>{});
	celero::DoNotOptimizeAway(result);
}


//------------------------------------------------------------------------------
// Transform
//------------------------------------------------------------------------------

BASELINE_F(TransformUnaryLight, std, RandomArray2Fixture<float>, 25, 500) {
	const auto result = std::transform(array1.begin(), array1.end(), array2.begin(), [](const auto& v) { return v * v; });
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(TransformUnaryLight, dspbb, RandomArray2Fixture<float>, 25, 500) {
	const auto result = kernels::Transform(array1.begin(), array1.end(), array2.begin(), [](const auto& v) { return v * v; });
	celero::DoNotOptimizeAway(result);
}

BASELINE_F(TransformBinaryLight, std, RandomArray3Fixture<float>, 25, 500) {
	const auto result = std::transform(array1.begin(), array1.end(), array2.begin(), array3.begin(), std::multiplies<>{});
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(TransformBinaryLight, dspbb, RandomArray3Fixture<float>, 25, 500) {
	const auto result = kernels::Transform(array1.begin(), array1.end(), array2.begin(), array3.begin(), std::multiplies<>{});
	celero::DoNotOptimizeAway(result);
}

static const auto heavyFunctor = [](const auto& v) { return kernels::math_functions::sin(v); };
static const auto heavyBinaryFunctor = [](const auto& u, const auto& v) { return kernels::math_functions::sin(u * u + v * v); };

BASELINE_F(TransformUnaryHeavy, std, RandomArray2Fixture<float>, 25, 500) {
	const auto result = std::transform(array1.begin(), array1.end(), array2.begin(), heavyFunctor);
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(TransformUnaryHeavy, dspbb, RandomArray2Fixture<float>, 25, 500) {
	const auto result = kernels::Transform(array1.begin(), array1.end(), array2.begin(), heavyFunctor);
	celero::DoNotOptimizeAway(result);
}

BASELINE_F(TransformBinaryHeavy, std, RandomArray3Fixture<float>, 25, 500) {
	const auto result = std::transform(array1.begin(), array1.end(), array2.begin(), array3.begin(), heavyBinaryFunctor);
	celero::DoNotOptimizeAway(result);
}

BENCHMARK_F(TransformBinaryHeavy, dspbb, RandomArray3Fixture<float>, 25, 500) {
	const auto result = kernels::Transform(array1.begin(), array1.end(), array2.begin(), array3.begin(), heavyBinaryFunctor);
	celero::DoNotOptimizeAway(result);
}