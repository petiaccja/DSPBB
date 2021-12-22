#include <dspbb/Kernels/Convolution.hpp>
#include <dspbb/Math/Convolution.hpp>

#include <array>
#include <celero/Celero.h>
#include <random>
#include <vector>


using namespace dspbb;


//------------------------------------------------------------------------------
// Input sizes for which to benchmark
//------------------------------------------------------------------------------

static constexpr std::array signalSizes = {
	2048,
	262144,
};

static constexpr std::array filterSizes = {
	1,
	2,
	4,
	6,
	8,
	12,
	16,
	32,
	64,
	128,
	256,
	512,
	1024,
	2048,
};

constexpr int complexityLimit = 32 * 1024 * 1024;


//------------------------------------------------------------------------------
// Fixtures to generate random input
//------------------------------------------------------------------------------

static std::minstd_rand rne;
static std::uniform_real_distribution<float> randomFloat(-1, 1);

template <class T, size_t SignalSize>
class ConvolutionFixture : public celero::TestFixture {
public:
	std::vector<ExperimentValue> getExperimentValues() const override {
		std::vector<ExperimentValue> experimentValues;
		for (auto& filterSize : filterSizes) {
			const auto complexity = filterSize * SignalSize;
			const auto iterations = complexityLimit / complexity;
			if (iterations > 0) {
				experimentValues.emplace_back(int64_t(filterSize), std::max(int64_t(1), int64_t(iterations)));
			}
		};
		return experimentValues;
	}

	void setUp(const ExperimentValue& experimentValue) override {
		size_t filterSize = experimentValue.Value;
		out = std::vector<T>(ConvolutionLength(SignalSize, filterSize, CONV_FULL));
		signal = std::vector<T>(SignalSize);
		filter = std::vector<T>(filterSize);
		std::array<T, 16> pattern;
		for (auto& v : pattern) {
			v = static_cast<T>(randomFloat(rne));
		}
		size_t index = 0;
		for (auto& v : signal) {
			v = pattern[index];
			index = (index + 1) % pattern.size();
		}
		for (auto& v : filter) {
			v = pattern[index];
			index = (index + 1) % pattern.size();
		}
	}

	std::vector<T> out;
	std::vector<T> signal;
	std::vector<T> filter;
};


//------------------------------------------------------------------------------
// Benchmarks
//------------------------------------------------------------------------------

using FixtureCache = ConvolutionFixture<float, signalSizes[0]>;

BASELINE_F(ConvolutionCache, naive, FixtureCache, 25, 500) {
	kernels::ConvolutionNaive(signal.begin(), signal.end(), filter.begin(), filter.end(), out.begin(), out.end(), 0);
	celero::DoNotOptimizeAway(out.front());
}

BENCHMARK_F(ConvolutionCache, slide, FixtureCache, 25, 500) {
	kernels::ConvolutionSlide(signal.begin(), signal.end(), filter.begin(), filter.end(), out.begin(), out.end(), 0);
	celero::DoNotOptimizeAway(out.front());
}

using FixtureLarge = ConvolutionFixture<float, signalSizes[1]>;

BASELINE_F(ConvolutionLarge, naive, FixtureLarge, 25, 500) {
	kernels::ConvolutionNaive(signal.begin(), signal.end(), filter.begin(), filter.end(), out.begin(), out.end(), 0);
	celero::DoNotOptimizeAway(out.front());
}

BENCHMARK_F(ConvolutionLarge, slide, FixtureLarge, 25, 500) {
	kernels::ConvolutionSlide(signal.begin(), signal.end(), filter.begin(), filter.end(), out.begin(), out.end(), 0);
	celero::DoNotOptimizeAway(out.front());
}
