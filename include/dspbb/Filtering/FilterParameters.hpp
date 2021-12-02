#pragma once

#include "../Generators/Spaces.hpp"
#include "../LTISystems/Systems.hpp"
#include "../Math/FFT.hpp"
#include "../Math/Functions.hpp"
#include "../Math/Statistics.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"

#include <optional>
#include <variant>


namespace dspbb {

template <class T>
struct LowpassParameters {
	T passbandEdge;
	T stopbandEdge;
	T passbandRipple;
	T stopbandAtten;
};

template <class T>
struct HighpassParameters {
	T stopbandEdge;
	T passbandEdge;
	T stopbandAtten;
	T passbandRipple;
};

template <class T>
struct BandpassParameters {
	T lowerStopbandEdge;
	T passbandLowerEdge;
	T passbandUpperEdge;
	T upperStopbandEdge;
	T lowerStopbandAtten;
	T passbandRipple;
	T upperStopbandAtten;
};

template <class T>
struct BandstopParameters {
	T lowerPassbandEdge;
	T stopbandLowerEdge;
	T stopbandUpperEdge;
	T upperPassbandEdge;
	T lowerPassbandRipple;
	T stopbandAtten;
	T upperPassbandRipple;
};

template <class T>
using FilterParameters = std::variant<
	std::monostate,
	LowpassParameters<T>,
	HighpassParameters<T>,
	BandpassParameters<T>,
	BandstopParameters<T>>;


namespace impl {
	template <class T>
	constexpr T threshold = T(0.5);
	constexpr size_t kernelSize = 5;

	struct Band {
		size_t first;
		bool pass;
	};

	template <class T>
	struct BandParameters {
		T lowerEdge;
		T upperEdge;
		T ripple;
	};

	template <class T>
	std::vector<Band> ExtractFilterBands(const SpectrumView<const T>& spectrum, T threshold, T noise = T(0.0005)) {
		const T upperThreshold = threshold + noise;
		const T lowerThreshold = threshold - noise;
		std::vector<Band> bands;
		auto it = spectrum.begin();
		size_t index = 0;
		do {
			const bool pass = *it > threshold;
			bands.push_back({ index, pass });

			const T triggerThreshold = pass ? upperThreshold : lowerThreshold;
			const T concludeThreshold = pass ? lowerThreshold : upperThreshold;
			const auto triggerIt = std::find_if(it, spectrum.end(), [pass, triggerThreshold](const T& bin) {
				return bin > triggerThreshold != pass;
			});
			const auto concludeIt = std::find_if(it, spectrum.end(), [pass, concludeThreshold](const T& bin) {
				return bin > concludeThreshold != pass;
			});

			const size_t triggerIndex = triggerIt - spectrum.begin();
			const size_t concludeIndex = concludeIt - spectrum.begin();
			index = (concludeIndex + triggerIndex) / 2;
			it = concludeIt;
		} while (it != spectrum.end());

		return bands;
	}

	template <class T, class Func>
	T Bisect(Func f, T a, T b) {
		auto lower = std::min(a, b);
		auto upper = std::max(a, b);
		int i = 100;
		do {
			const T c = (lower + upper) / T(2);
			const auto flower = f(lower);
			const auto fc = f(c);
			if (fc * flower > decltype(fc)(0)) {
				lower = c;
			}
			else {
				upper = c;
			}
		} while (i-- > 0);
		return lower;
	}

	template <class T>
	T FitErrorDerivative(T x, T y, T p) {
		const auto q = decltype(p)(1) / p;
		return (y - std::erf(q * x)) * x * std::exp(-q * q * x * x);
	}

	template <class T>
	T FitLoss(const SpectrumView<const T>& y, T p, bool invertX, bool invertY, T baseY) {
		const size_t count = y.Size();
		T sum = T(0);
		for (size_t i = 0; i < y.Size(); ++i) {
			const T xi = !invertX ? T(i) : T(count - i - 1);
			const T yi = !invertY ? (y[i] - baseY) / (T(1) - baseY) : (baseY - y[i]) / baseY;
			sum += FitErrorDerivative(xi, yi, p);
		}
		return sum;
	}

	template <class T>
	auto FitBand(const SpectrumView<const T>& band, bool pass, bool left, T threshold) {
		const size_t size = band.Size();

		const bool invertX = !left;
		const bool invertY = !pass;
		const auto loss = [invertX, invertY, &band, threshold](T p) {
			return FitLoss(band, p, invertX, invertY, threshold);
		};

		const auto maxIndex = (pass ? std::max_element(band.begin(), band.end()) : std::min_element(band.begin(), band.end())) - band.begin();
		const T pLower = 0.5f;
		const T pUpper = T(0.66) * (!invertX ? T(maxIndex) : T(size - maxIndex));
		const auto pOptimal = Bisect(loss, pLower, pUpper);

		return pOptimal;
	}

	template <class T>
	std::pair<std::optional<T>, std::optional<T>> FindBandEdgesFit(const SpectrumView<const T>& band, bool findLeft, bool findRight, bool passBand, T threshold) {
		constexpr T paramToIndex = T(1.5);

		const std::optional<T> paramLeft = findLeft ? std::optional<T>{ impl::FitBand(band, passBand, true, threshold) } : std::optional<T>{};
		const std::optional<T> paramRight = findRight ? std::optional<T>{ impl::FitBand(band, passBand, false, threshold) } : std::optional<T>{};

		const std::optional<T> indexLeft = paramLeft ? std::optional<T>{ std::min(T(band.Size()), paramToIndex * paramLeft.value()) } : std::optional<T>{};
		const std::optional<T> indexRight = paramRight ? std::optional<T>{ std::max(T(0), T(band.Size()) - paramToIndex * paramRight.value()) } : std::optional<T>{};

		return { indexLeft, indexRight };
	}

	template <class Iter, class Compare>
	Iter FindFirstExtremum(Iter first, Iter last, size_t kernelSize, Compare compare) {
		using T = std::decay_t<decltype(*first)>;
		for (auto baseIt = first; baseIt != last; ++baseIt) {
			Iter strideIt = baseIt;
			Iter centerIt = last;
			size_t strideIdx;

			T extremeValue = *baseIt;
			std::optional<T> centerValue;

			for (strideIdx = 0; strideIdx < kernelSize && strideIt != last; ++strideIdx, ++strideIt) {
				extremeValue = std::max(extremeValue, *strideIt, compare);
				if (strideIdx == kernelSize / 2) {
					centerValue = *strideIt;
					centerIt = strideIt;
				}
			}

			if (centerValue && extremeValue == centerValue.value() && strideIdx == kernelSize) {
				return centerIt;
			}
		}
		return last;
	}

	template <class T>
	std::pair<std::optional<T>, std::optional<T>> FindBandEdgesRipple(const SpectrumView<const T>& band, bool findLeft, bool findRight, bool passBand) {
		std::optional<T> edgeLeft;
		std::optional<T> edgeRight;
		if (findLeft) {
			const auto extremumIt = passBand ? FindFirstExtremum(band.begin(), band.end(), kernelSize, std::less<T>{}) : FindFirstExtremum(band.begin(), band.end(), kernelSize, std::greater<T>{});
			if (extremumIt != band.end()) {
				edgeLeft = T(extremumIt - band.begin());
			}
		}
		if (findRight) {
			const auto extremumIt = passBand ? FindFirstExtremum(band.rbegin(), band.rend(), kernelSize, std::less<T>{}) : FindFirstExtremum(band.rbegin(), band.rend(), kernelSize, std::greater<T>{});
			if (extremumIt != band.rend()) {
				edgeRight = T(extremumIt.base() - band.begin());
			}
		}

		return { edgeLeft, edgeRight };
	}

	template <class T>
	std::optional<T> MeasureBandRipple(const SpectrumView<const T>& band, bool passBand) {
		const T target = T(passBand);
		T maxDifference = T(-1);

		auto it = band.begin();
		while (it != band.end()) {
			it = FindFirstExtremum(it, band.end(), kernelSize, std::less<T>{});
			if (it != band.end()) {
				maxDifference = std::max(maxDifference, std::abs(target - *it));
			}
		}

		it = band.begin();
		while (it != band.end()) {
			it = FindFirstExtremum(it, band.end(), kernelSize, std::greater<T>{});
			if (it != band.end()) {
				maxDifference = std::max(maxDifference, std::abs(target - *it));
			}
		}

		return maxDifference > 0 ? std::optional<T>{ maxDifference } : std::optional<T>{};
	}

	template <class T, class Compare>
	const std::optional<T>& MaxOptional(const std::optional<T>& lhs, const std::optional<T>& rhs, Compare compare) {
		if (lhs && rhs) {
			return compare(lhs.value(), rhs.value()) ? lhs : rhs;
		}
		return lhs ? lhs : rhs;
	}
	template <class T>
	const std::optional<T>& MinOptional(const std::optional<T>& lhs, const std::optional<T>& rhs) {
		return MaxOptional(lhs, rhs, std::less<T>{});
	}
	template <class T>
	const std::optional<T>& MaxOptional(const std::optional<T>& lhs, const std::optional<T>& rhs) {
		return MaxOptional(lhs, rhs, std::greater_equal<T>{});
	}

	template <class T>
	std::vector<BandParameters<T>> ParametrizeFilterBands(const SpectrumView<const T>& spectrum, const std::vector<Band>& bands, T threshold) {
		std::vector<BandParameters<T>> parameters;
		for (auto it = bands.begin(); it != bands.end(); ++it) {
			auto next = it;
			++next;
			const size_t bandFirstIndex = it->first;
			const size_t bandLastIndex = next != bands.end() ? next->first : spectrum.Size();
			const auto bandFirstIt = spectrum.begin() + bandFirstIndex;
			const auto bandLastIt = spectrum.begin() + bandLastIndex;

			const SpectrumView<const T> bandSamples{ bandFirstIt, bandLastIt };
			const bool isPassBand = it->pass;
			const bool findLeft = it != bands.begin();
			const bool findRight = next != bands.end();

			const auto [edgeFitLeft, edgeFitRight] = impl::FindBandEdgesFit(bandSamples, findLeft, findRight, isPassBand, threshold);
			const auto [edgeRippleLeft, edgeRippleRight] = impl::FindBandEdgesRipple(bandSamples, findLeft, findRight, isPassBand);
			const auto ripple = impl::MeasureBandRipple(bandSamples, isPassBand);
			const auto edgeLeft = MinOptional(edgeFitLeft, edgeRippleLeft);
			const auto edgeRight = MaxOptional(edgeFitRight, edgeRippleRight);

			const T normalizedEdgeLeft = (T(bandFirstIndex) + edgeLeft.value_or(T(0))) / T(spectrum.Size());
			const T normalizedEdgeRight = (T(bandFirstIndex) + edgeRight.value_or(T(bandSamples.Size()))) / T(spectrum.Size());
			parameters.push_back({ normalizedEdgeLeft, normalizedEdgeRight, ripple.value_or(T(0)) });
		}
		return parameters;
	}

	template <class... BandResponses>
	bool HasBands(const std::vector<Band>& bands, BandResponses... desiredPasses) {
		std::array<bool, sizeof...(BandResponses)> desiredPassesArray{ desiredPasses... };
		if (desiredPassesArray.size() != bands.size()) {
			return false;
		}
		for (size_t i = 0; i < desiredPassesArray.size(); ++i) {
			if (desiredPassesArray[i] != bands[i].pass) {
				return false;
			}
		}
		return true;
	}

	template <class T>
	std::optional<LowpassParameters<T>> ExtractLowpassParameters(const std::vector<Band>& bands, const std::vector<BandParameters<T>>& parameters) {
		if (HasBands(bands, true, false)) {
			return { LowpassParameters<T>{ parameters[0].upperEdge,
										   parameters[1].lowerEdge,
										   parameters[0].ripple,
										   parameters[1].ripple } };
		}
		return {};
	}

	template <class T>
	std::optional<HighpassParameters<T>> ExtractHighpassParameters(const std::vector<Band>& bands, const std::vector<BandParameters<T>>& parameters) {
		if (HasBands(bands, false, true)) {
			return { HighpassParameters<T>{ parameters[0].upperEdge,
											parameters[1].lowerEdge,
											parameters[0].ripple,
											parameters[1].ripple } };
		}
		return {};
	}

	template <class T>
	std::optional<BandpassParameters<T>> ExtractBandpassParameters(const std::vector<Band>& bands, const std::vector<BandParameters<T>>& parameters) {
		if (HasBands(bands, false, true, false)) {
			return { BandpassParameters<T>{ parameters[0].upperEdge,
											parameters[1].lowerEdge,
											parameters[1].upperEdge,
											parameters[2].lowerEdge,
											parameters[0].ripple,
											parameters[1].ripple,
											parameters[2].ripple } };
		}
		return {};
	}

	template <class T>
	std::optional<BandstopParameters<T>> ExtractBandstopParameters(const std::vector<Band>& bands, const std::vector<BandParameters<T>>& parameters) {
		if (HasBands(bands, true, false, true)) {
			return { BandstopParameters<T>{ parameters[0].upperEdge,
											parameters[1].lowerEdge,
											parameters[1].upperEdge,
											parameters[2].lowerEdge,
											parameters[0].ripple,
											parameters[1].ripple,
											parameters[2].ripple } };
		}
		return {};
	}

	template <class T>
	FilterParameters<T> ExtractFilterParameters(const std::vector<Band>& bands, const std::vector<BandParameters<T>>& parameters) {
		if (auto lowpassParameters = ExtractLowpassParameters(bands, parameters)) {
			return lowpassParameters.value();
		}
		if (auto highpassParameters = ExtractHighpassParameters(bands, parameters)) {
			return highpassParameters.value();
		}
		if (auto bandpassParameters = ExtractBandpassParameters(bands, parameters)) {
			return bandpassParameters.value();
		}
		if (auto bandstopParameters = ExtractBandstopParameters(bands, parameters)) {
			return bandstopParameters.value();
		}
		return {};
	}
} // namespace impl


template <class T>
LowpassParameters<T> ParametrizeLowpassFilter(const SpectrumView<const T>& response) {
	const std::vector<impl::Band> bands = impl::ExtractFilterBands(response, impl::threshold<T>);
	const std::vector<impl::BandParameters<T>> bandParams = impl::ParametrizeFilterBands(response, bands, impl::threshold<T>);
	auto filterParams = impl::ExtractLowpassParameters(bands, bandParams);
	if (!filterParams) {
		throw std::invalid_argument("Not a low-pass filter.");
	}
	return filterParams.value();
}

template <class T>
HighpassParameters<T> ParametrizeHighpassFilter(const SpectrumView<const T>& response) {
	const std::vector<impl::Band> bands = impl::ExtractFilterBands(response, impl::threshold<T>);
	const std::vector<impl::BandParameters<T>> bandParams = impl::ParametrizeFilterBands(response, bands, impl::threshold<T>);
	auto filterParams = impl::ExtractHighpassParameters(bands, bandParams);
	if (!filterParams) {
		throw std::invalid_argument("Not a high-pass filter.");
	}
	return filterParams.value();
}

template <class T>
BandpassParameters<T> ParametrizeBandpassFilter(const SpectrumView<const T>& response) {
	const std::vector<impl::Band> bands = impl::ExtractFilterBands(response, impl::threshold<T>);
	const std::vector<impl::BandParameters<T>> bandParams = impl::ParametrizeFilterBands(response, bands, impl::threshold<T>);
	auto filterParams = impl::ExtractBandpassParameters(bands, bandParams);
	if (!filterParams) {
		throw std::invalid_argument("Not a band-pass filter.");
	}
	return filterParams.value();
}

template <class T>
BandstopParameters<T> ParametrizeBandstopFilter(const SpectrumView<const T>& response) {
	const std::vector<impl::Band> bands = impl::ExtractFilterBands(response, impl::threshold<T>);
	const std::vector<impl::BandParameters<T>> bandParams = impl::ParametrizeFilterBands(response, bands, impl::threshold<T>);
	auto filterParams = impl::ExtractBandstopParameters(bands, bandParams);
	if (!filterParams) {
		throw std::invalid_argument("Not a band-stop filter.");
	}
	return filterParams.value();
}

template <class T>
FilterParameters<T> ParametrizeFilter(const SpectrumView<const T>& response) {
	const std::vector<impl::Band> bands = impl::ExtractFilterBands(response, impl::threshold<T>);
	const std::vector<impl::BandParameters<T>> bandParams = impl::ParametrizeFilterBands(response, bands, impl::threshold<T>);
	auto filterParams = impl::ExtractFilterParameters(bands, bandParams);
	return filterParams;
}


template <class T>
LowpassParameters<T> ParametrizeLowpassFilter(const Spectrum<T>& response) {
	return ParametrizeLowpassFilter(AsView(response));
}

template <class T>
HighpassParameters<T> ParametrizeHighpassFilter(const Spectrum<T>& response) {
	return ParametrizeHighpassFilter(AsView(response));
}

template <class T>
BandpassParameters<T> ParametrizeBandpassFilter(const Spectrum<T>& response) {
	return ParametrizeBandpassFilter(AsView(response));
}

template <class T>
BandstopParameters<T> ParametrizeBandstopFilter(const Spectrum<T>& response) {
	return ParametrizeBandstopFilter(AsView(response));
}

template <class T>
FilterParameters<T> ParametrizeFilter(const Spectrum<T>& response) {
	return ParametrizeFilter(AsView(response));
}


namespace impl {
	inline size_t FrequencyResponseFftSize(size_t impulseSize, size_t desiredGridSize) {
		return std::max(impulseSize, 2 * desiredGridSize - 1);
	}
} // namespace impl


template <class T>
auto FrequencyResponse(const BasicSignalView<const T, TIME_DOMAIN>& impulse, size_t gridSizeHint = 0) {
	const size_t gridSize = gridSizeHint > 0 ? gridSizeHint : impulse.Size() * 10;
	const size_t paddedSize = impl::FrequencyResponseFftSize(impulse.Size(), gridSize);

	BasicSignal<T, TIME_DOMAIN> padded(paddedSize, T(0));
	std::copy(impulse.begin(), impulse.end(), padded.begin());

	auto spectrum = Fft(padded, FFT_HALF);
	auto amplitude = Abs(spectrum);
	std::replace(spectrum.begin(), spectrum.end(), std::complex<T>{ T(0) }, std::complex<T>{ T(1) });
	auto phase = Arg(spectrum);
	return std::make_pair(std::move(amplitude), std::move(phase));
}

template <class T>
auto FrequencyResponse(const BasicSignal<T, TIME_DOMAIN>& impulse, size_t gridSizeHint = 0) {
	return FrequencyResponse(AsView(impulse), gridSizeHint);
}


namespace impl {

	template <class T, class System>
	auto FrequencyResponse(const System& sys, size_t gridSizeHint = 0) -> std::pair<Spectrum<T>, Spectrum<T>> {
		const size_t order = sys.Order();
		const size_t gridSize = gridSizeHint > 0 ? gridSizeHint : (1 + order) * 20;

		Spectrum<T> amplitude(gridSize);
		Spectrum<T> phase(gridSize);

		LinSpace(amplitude, 0.0f, pi_v<T>, true);

		for (size_t i = 0; i < gridSize; ++i) {
			const T freq = amplitude[i];
			const std::complex<T> point = std::polar(T(1), freq);
			const std::complex<T> response = sys(point);
			amplitude[i] = std::abs(response);
			phase[i] = std::arg(response);
		}

		return { std::move(amplitude), std::move(phase) };
	}

} // namespace impl


template <class T>
auto FrequencyResponse(const DiscreteZeroPoleGain<T>& zpk, size_t gridSizeHint = 0) {
	return impl::FrequencyResponse<T>(zpk, gridSizeHint);
}

template <class T>
auto FrequencyResponse(const CascadedBiquad<T>& biquad, size_t gridSizeHint = 0) {
	return impl::FrequencyResponse<T>(biquad, gridSizeHint);
}

template <class T>
auto FrequencyResponse(const DiscreteTransferFunction<T>& tf, size_t gridSizeHint = 0) {
	return impl::FrequencyResponse<T>(tf, gridSizeHint);
}


} // namespace dspbb