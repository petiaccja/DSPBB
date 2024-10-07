#include "dspbb/Filtering/FIR.hpp"
#include "dspbb/Filtering/IIR.hpp"

#include <array>
#include <celero/Celero.h>
#include <random>
#include <vector>

using namespace dspbb;



//------------------------------------------------------------------------------
// Input sizes for which to benchmark
//------------------------------------------------------------------------------

constexpr size_t signalSize = 262144;
constexpr size_t maxFirOrder = 4096;
constexpr size_t maxIirDirectOrder = 8;
constexpr size_t maxIirCascadeOrder = 16;

constexpr size_t complexityLimit = signalSize * maxFirOrder;

//------------------------------------------------------------------------------
// Fixtures to generate random input
//------------------------------------------------------------------------------

static std::minstd_rand rne;
static std::uniform_real_distribution<float> randomFloat(-1, 1);


template <class T, int64_t MinOrder = 1, int64_t MaxOrder = maxFirOrder, int64_t MaxIter = 16384>
class FirFilterFixture : public celero::TestFixture {
public:
	std::vector<std::shared_ptr<ExperimentValue>> getExperimentValues() const override {
		std::vector<std::shared_ptr<ExperimentValue>> experimentValues;
		for (size_t filterOrder = MinOrder; filterOrder <= MaxOrder; filterOrder *= 2) {
			size_t iterations = std::clamp(int64_t(complexityLimit / (filterOrder * signalSize)), int64_t(1), MaxIter);
			experimentValues.emplace_back(std::make_shared<ExperimentValue>(int64_t(filterOrder), iterations));
		};
		return experimentValues;
	}

	void setUp(const ExperimentValue* experimentValue) override {
		size_t filterOrder = experimentValue->Value;
		out = Signal<T>(ConvolutionLength(signalSize, filterOrder + 1, CONV_FULL));
		signal = Signal<T>(signalSize);
		filter = Signal<T>(filterOrder + 1);
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

	Signal<T> out;
	Signal<T> signal;
	Signal<T> filter;
};


template <class T, int64_t MaxOrder>
class DesignFilterFixture : public celero::TestFixture {
public:
	std::vector<std::shared_ptr<ExperimentValue>> getExperimentValues() const override {
		std::vector<std::shared_ptr<ExperimentValue>> experimentValues;
		for (size_t filterOrder = 1; filterOrder <= MaxOrder; filterOrder += 1) {
			size_t iterations = 50 / filterOrder;
			experimentValues.emplace_back(std::make_shared<ExperimentValue>(int64_t(filterOrder), iterations));
		};
		return experimentValues;
	}

	void setUp(const ExperimentValue* experimentValue) override {
		size_t filterOrder = experimentValue->Value;
		out = Signal<T>(signalSize);
		signal = Signal<T>(signalSize);
		std::array<T, 16> pattern;
		for (auto& v : pattern) {
			v = static_cast<T>(randomFloat(rne));
		}
		size_t index = 0;
		for (auto& v : signal) {
			v = pattern[index];
			index = (index + 1) % pattern.size();
		}

		filter.zeros.resize(filterOrder, 0);
		filter.poles.resize(filterOrder, 0);
		filter.gain = T(1.0) + T(0.001) * pattern[0];
		for (auto& v : filter.zeros.real_roots()) {
			v = T(-0.95) + T(0.001) * pattern[index];
			index = (index + 1) % pattern.size();
		}
		for (auto& v : filter.poles.real_roots()) {
			v = T(-0.90) + T(0.001) * pattern[index];
			index = (index + 1) % pattern.size();
		}
	}

	Signal<T> out;
	Signal<T> signal;
	DiscreteZeroPoleGain<T> filter;
};


//------------------------------------------------------------------------------
// Benchmarks
//------------------------------------------------------------------------------
using BaselineFixture = FirFilterFixture<float, 1, 1, 8192>;
using ConvFixture = FirFilterFixture<float, 1, maxFirOrder, 16>;
using OlaFixture = FirFilterFixture<float, 32, maxFirOrder, 16>;
using TfFixture = DesignFilterFixture<float, maxIirDirectOrder>;
using CascadeFixture = DesignFilterFixture<float, maxIirCascadeOrder>;

BASELINE_F(ApplyFilter, gain, BaselineFixture, 25, 1) {
	Multiply(AsView(out).subsignal(0, signal.size()), signal, filter[0]);
	celero::DoNotOptimizeAway(out[0]);
}

BENCHMARK_F(ApplyFilter, fir_conv, FirFilterFixture<float>, 25, 1) {
	Filter(out, signal, filter, CONV_FULL, FILTER_CONV);
	celero::DoNotOptimizeAway(out[0]);
}

BENCHMARK_F(ApplyFilter, fir_ola, OlaFixture, 25, 1) {
	Filter(out, signal, filter, CONV_FULL, FILTER_OLA);
	celero::DoNotOptimizeAway(out[0]);
}

BENCHMARK_F(ApplyFilter, iir_df_i, TfFixture, 25, 1) {
	const auto realization = TransferFunction{ filter };
	DirectFormI<float> state{ realization.order() };
	Filter(out, signal, realization, state);
	celero::DoNotOptimizeAway(out[0]);
}

BENCHMARK_F(ApplyFilter, iir_df_ii, TfFixture, 25, 1) {
	const auto realization = TransferFunction{ filter };
	DirectFormII<float> state{ realization.order() };
	Filter(out, signal, realization, state);
	celero::DoNotOptimizeAway(out[0]);
}

BENCHMARK_F(ApplyFilter, iir_cascade, CascadeFixture, 25, 1) {
	const auto realization = CascadedBiquad{ filter };
	CascadedForm<float> state{ realization.order() };
	Filter(out, signal, realization, state);
	celero::DoNotOptimizeAway(out[0]);
}