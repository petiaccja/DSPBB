#pragma once

#include "../Primitives/Signal.hpp"
#include "../Utility/Numbers.hpp"

#include <Eigen/Dense>
#include <Eigen/QR>
#include <iostream>

namespace dspbb {


namespace impl {

	template <class T>
	auto CoefficientMatrix(size_t filterLength, size_t gridSize) {
		Eigen::Matrix<T, Eigen::Dynamic, 1> col;
		Eigen::Matrix<T, 1, Eigen::Dynamic> row;

		col.resize(gridSize);
		row.resize(filterLength);

		for (size_t i = 0; i < gridSize; ++i) {
			col(i) = T(i);
		}
		col *= T(1) / T(gridSize - 1) * pi_v<T>;

		for (size_t i = 0; i < filterLength; ++i) {
			row(i) = T(i);
		}

		Eigen::MatrixX<T> coefficientMatrix = col * row;

		for (size_t row = 0; row < gridSize; ++row) {
			coefficientMatrix(row, 0) = T(1);
			for (size_t col = 1; col < filterLength; ++col) {
				coefficientMatrix(row, col) = 2 * std::cos(coefficientMatrix(row, col));
			}
		}

		return coefficientMatrix;
	}

	template <class T, class Func>
	auto WeightMatrix(size_t gridSize, const Func& weight) {
		Eigen::DiagonalMatrix<T, Eigen::Dynamic> weightMatrix;
		weightMatrix.resize(gridSize);
		for (size_t i = 0; i < gridSize; ++i) {
			weightMatrix.diagonal()(i) = weight(T(i) / T(gridSize - 1));
		}
		return weightMatrix;
	}

	template <class T, class Func>
	auto ResponseVector(size_t gridSize, const Func& response) {
		Eigen::Matrix<T, Eigen::Dynamic, 1> responseVector;
		responseVector.resize(gridSize);
		for (size_t i = 0; i < gridSize; ++i) {
			responseVector(i) = T(i);
		}
		responseVector *= T(1) / T(gridSize - 1);
		for (size_t i = 0; i < gridSize; ++i) {
			responseVector(i) = response(responseVector(i));
		}
		return responseVector;
	}

} // namespace impl


template <class SignalR, class ResponseFunc, class WeightFunc, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirLeastSquares(SignalR&& coefficients, ResponseFunc responseFunc, WeightFunc weightFunc, size_t gridSize = 0) {
	using R = typename std::decay_t<SignalR>::value_type;
	using T = remove_complex_t<R>;

	const size_t filterLength = (coefficients.Size() + 1) / 2;
	gridSize = gridSize == 0 ? 4 * filterLength : std::min(filterLength, gridSize);

	const auto coefficientMatrix = impl::CoefficientMatrix<T>(filterLength, gridSize);
	const auto weightMatrix = impl::WeightMatrix<T>(gridSize, weightFunc);
	const auto responseVector = impl::ResponseVector<T>(gridSize, responseFunc);

	Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixX<T>> decomp{ weightMatrix * coefficientMatrix };
	const Eigen::VectorX<T> halfFilter = decomp.solve(weightMatrix * responseVector);
	for (size_t i = 0; i < filterLength; ++i) {
		coefficients[i] = halfFilter(filterLength - i - 1);
		coefficients[i + filterLength - 1] = halfFilter[i];
	}
}

} // namespace dspbb