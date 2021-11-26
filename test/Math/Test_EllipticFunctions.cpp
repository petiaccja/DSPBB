#include "../TestUtils.hpp"

#include <dspbb/Math/EllipticFunctions.hpp>
#include <dspbb/Utility/Numbers.hpp>

#include <catch2/catch.hpp>


using namespace dspbb;
using namespace std::complex_literals;

//------------------------------------------------------------------------------
// Carlson symmetric forms
//------------------------------------------------------------------------------

TEST_CASE("Carlson RF real", "[Elliptic functions]") {
	const auto r = CarlsonRF(0.0, 1.0, 2.0);
	REQUIRE(r == Approx(1.31102877714605990523).epsilon(1e-15));
}

TEST_CASE("Carlson RF complex", "[Elliptic functions]") {
	const auto r = CarlsonRF(1.0 + 2.0i, 0.6 + 3.0i, 1.55);
	REQUIRE(r == ApproxComplex(0.6252032195312656 - 0.2990455989714608i));
}



//------------------------------------------------------------------------------
// Jacobi amplitude function
//------------------------------------------------------------------------------

TEST_CASE("Elliptic AM x=0", "[Elliptic functions]") {
	const auto r = EllipticAM(0, 0.4);
	REQUIRE(r == Approx(0).epsilon(1e-15f));
}

TEST_CASE("Elliptic AM x=anything", "[Elliptic functions]") {
	const auto r = EllipticAM(3.6, 0.4);
	REQUIRE(r == Approx(3.460737469702738).epsilon(1e-15f));
}

TEST_CASE("Elliptic AM x=K(k)", "[Elliptic functions]") {
	const auto k = 0.4;
	const auto r = EllipticAM(EllipticK(k), k);
	REQUIRE(r == Approx(pi_v<double> / 2).epsilon(1e-15f));
}

TEST_CASE("Elliptic AM k=0", "[Elliptic functions]") {
	const auto k = 0.0;
	const auto x = 1.3453452364625342563542;
	const auto r = EllipticAM(x, k);
	REQUIRE(r == Approx(x).epsilon(1e-15f));
}

TEST_CASE("Elliptic AM k=1", "[Elliptic functions]") {
	const auto k = 1.0;
	const auto x = 1.3453452364625342563542;
	const auto r = EllipticAM(x, k);
	REQUIRE(r == Approx(1.0612177100827493027662).epsilon(1e-15f));
}

TEST_CASE("Elliptic AM k=tiny", "[Elliptic functions]") {
	const auto k = std::numeric_limits<double>::epsilon();
	const auto x = 1.3453452364625342563542;
	const auto r = EllipticAM(x, k);
	REQUIRE(r == Approx(1.34534523646253419376).epsilon(1e-15f));
}

TEST_CASE("Elliptic AM k=denormal", "[Elliptic functions]") {
	const auto k = std::nextafter(0.0, 1.0);
	const auto x = 1.3453452364625342563542;
	const auto r = EllipticAM(x, k);
	REQUIRE(r == Approx(1.34534523646253419376).epsilon(1e-15f));
}

TEST_CASE("Elliptic AM k=1-tiny", "[Elliptic functions]") {
	const auto k = std::nextafter(1.0, 0.0);
	const auto x = 1.3453452364625342563542;
	const auto r = EllipticAM(x, k);
	REQUIRE(r == Approx(1.06121771008274933422041).epsilon(1e-15f));
}


//------------------------------------------------------------------------------
// Jacobi elliptic functions
//------------------------------------------------------------------------------


TEST_CASE("Elliptic SNCNDN inverse definition real", "[Elliptic functions]") {
	const auto k = 0.1;
	const auto x = 0.8;
	const auto [s, c, d] = EllipticSNCNDN(x, k);
	const auto xs = EllipticArcSN(s, k);
	const auto xc = EllipticArcCN(c, k);
	const auto xd = EllipticArcDN(d, k);
	REQUIRE(xs == Approx(x).epsilon(1e-15f));
	REQUIRE(xc == Approx(x).epsilon(1e-15f));
	REQUIRE(xd == Approx(x).epsilon(3e-13f));
}

TEST_CASE("Elliptic SNCNDN inverse definition complex", "[Elliptic functions]") {
	const auto k = 0.4;
	const auto x = 0.8 + 0.1i;
	const auto [s, c, d] = EllipticSNCNDN(x, k);
	const auto xs = EllipticArcSN(s, k);
	const auto xc = EllipticArcCN(c, k);
	const auto xd = EllipticArcDN(d, k);
	REQUIRE(xs == ApproxComplex(x).epsilon(1e-15f));
	REQUIRE(xc == ApproxComplex(x).epsilon(1e-15f));
	REQUIRE(xd == ApproxComplex(x).epsilon(3e-13f));
}

TEST_CASE("Elliptic SNCNDN special values", "[Elliptic functions]") {
	struct Record {
		std::complex<double> z;
		double k;
		std::complex<double> ysn;
		std::complex<double> ycn;
		std::complex<double> ydn;
	};

	const double k = 0.17;
	const double kp = std::sqrt(1 - k * k);
	const double K = EllipticK(k);
	const double Kp = EllipticK(kp);
	constexpr double inf = std::numeric_limits<double>::infinity();

	const std::array values = {
		Record{ 0.0, k, 0.0, 1.0, 1.0 },
		Record{ K, k, 1.0, 0.0, kp },
		Record{ K + 1.0i * Kp, k, 1 / k, -1.0i * kp / k, 0 },
		Record{ 1.0i * Kp, k, inf, inf, inf },
		Record{ 2.0 * K, k, 0, -1, 1 },
		Record{ 2.0 * K + 2.0i * Kp, k, 0, 1, -1 },
		Record{ 2.0i * Kp, k, 0, -1, -1 },
	};

	for (auto& value : values) {
		INFO("z=" << value.z << "|"
				  << "k=" << value.k)
		auto [ysn, ycn, ydn] = EllipticSNCNDN(value.z, value.k);
		ysn = abs(ysn) > 1e+8 ? inf : ysn;
		ycn = abs(ysn) > 1e+8 ? inf : ycn;
		ydn = abs(ysn) > 1e+8 ? inf : ydn;
		REQUIRE(ysn == ApproxComplex(value.ysn).margin(1e-15f));
		REQUIRE(ycn == ApproxComplex(value.ycn).margin(1e-15f));
		REQUIRE(ydn == ApproxComplex(value.ydn).margin(1e-15f));
	}
}


const std::array<std::complex<double>, 6> zs = {
	0.0,
	0.0 + 1.0i,
	1.0,
	pi_v<double>,
	pi_v<double> + 2.0i,
	8.7 - 2.15i
};

TEST_CASE("Elliptic SNCNDN degeneration at k=0", "[Elliptic functions]") {
	for (auto& z : zs) {
		const auto [ysn, ycn, ydn] = EllipticSNCNDN(z, 0.0);
		REQUIRE(ysn == ApproxComplex(std::sin(z)).margin(1e-15f));
		REQUIRE(ycn == ApproxComplex(std::cos(z)).margin(1e-15f));
		REQUIRE(ydn == ApproxComplex(1.0).margin(1e-15f));
	}
}

TEST_CASE("Elliptic SNCNDN degeneration at k=1", "[Elliptic functions]") {
	for (auto& z : zs) {
		const auto [ysn, ycn, ydn] = EllipticSNCNDN(z, 1.0);
		REQUIRE(ysn == ApproxComplex(std::tanh(z)).margin(1e-15f));
		REQUIRE(ycn == ApproxComplex(1.0 / std::cosh(z)).margin(1e-15f));
		REQUIRE(ydn == ApproxComplex(1.0 / std::cosh(z)).margin(1e-15f));
	}
}

TEST_CASE("Elliptic SNCNDN degeneration at k=0 + tiny", "[Elliptic functions]") {
	const double k = std::nextafter(0.0, 1.0);
	for (auto& z : zs) {
		const auto [ysn, ycn, ydn] = EllipticSNCNDN(z, k);
		REQUIRE(ysn == ApproxComplex(std::sin(z)).margin(1e-15f));
		REQUIRE(ycn == ApproxComplex(std::cos(z)).margin(1e-15f));
		REQUIRE(ydn == ApproxComplex(1.0).margin(1e-15f));
	}
}

TEST_CASE("Elliptic SNCNDN degeneration at k=1 - tiny", "[Elliptic functions]") {
	const double k = std::nextafter(1.0, 0.0);
	for (auto& z : zs) {
		const auto [ysn, ycn, ydn] = EllipticSNCNDN(z, k);
		REQUIRE(ysn == ApproxComplex(std::tanh(z)).margin(1e-15f));
		REQUIRE(ycn == ApproxComplex(1.0 / std::cosh(z)).margin(1e-15f));
		REQUIRE(ydn == ApproxComplex(1.0 / std::cosh(z)).margin(1e-15f));
	}
}

TEST_CASE("Elliptic SNCNDN half argument identity", "[Elliptic functions]") {
	const double k = 0.17;
	for (auto& z : zs) {
		INFO("z=" << z)
		const auto [halfsn, halfcn, halfdn] = EllipticSNCNDN(z / 2.0, k);
		const auto [fullsn, fullcn, fulldn] = EllipticSNCNDN(z, k);
		REQUIRE(halfsn * halfsn == ApproxComplex((1.0 - fullcn) / (1.0 + fulldn)).margin(1e-15));
		REQUIRE(halfcn * halfcn == ApproxComplex((fulldn + k * k * fullcn - (1.0 - k * k)) / (k * k * (1.0 + fullcn))).margin(1e-15));
		REQUIRE(halfdn * halfdn == ApproxComplex((k * k * fullcn + fulldn + (1.0 - k * k)) / (1.0 + fulldn)).margin(1e-15));
	}
}

TEST_CASE("Elliptic SNCNDN imaginary transformations identity", "[Elliptic functions]") {
	const double k = 0.17;
	for (auto& z : zs) {
		INFO("z=" << z)
		const auto [rotsn, rotcn, rotdn] = EllipticSNCNDN(1.0i * z, k);
		const auto [ysn, ycn, ydn] = EllipticSNCNDN(z, std::sqrt(1.0 - k * k));
		REQUIRE(rotsn == ApproxComplex(1.0i * ysn / ycn).margin(1e-15));
		REQUIRE(rotcn == ApproxComplex(1.0 / ycn).margin(1e-15));
		REQUIRE(rotdn == ApproxComplex(ydn / ycn).margin(1e-15));
	}
}

TEST_CASE("Elliptic SNCNDN descending Landen transformations identity", "[Elliptic functions]") {
	const double k = 0.17;
	const double kp = std::sqrt(1.0 - k * k);
	const double k1 = (1.0 - kp) / (1.0 + kp);
	for (auto& z : zs) {
		INFO("z=" << z)
		const auto [ysn, ycn, ydn] = EllipticSNCNDN(z, k);
		const auto [ysn1, ycn1, ydn1] = EllipticSNCNDN(z / (1.0 + k1), k1);
		REQUIRE(ysn == ApproxComplex((1.0 + k1) * ysn1 / (1.0 + k1 * ysn1 * ysn1)).margin(1e-15));
		REQUIRE(ycn == ApproxComplex(ycn1 * ydn1 / (1.0 + k1 * ysn1 * ysn1)).margin(1e-15));
		REQUIRE(ydn == ApproxComplex((ydn1 * ydn1 - 1.0 + k1) / (1.0 + k1 - ydn1 * ydn1)).margin(1e-15));
	}
}

TEST_CASE("Elliptic SNCNDN ascending Landen transformations identity", "[Elliptic functions]") {
	const double k = 0.17;
	const double k2 = 2.0 * std::sqrt(k) / (1.0 + k);
	const double k2p = (1.0 - k) / (1.0 + k);
	for (auto& z : zs) {
		INFO("z=" << z)
		const auto [ysn, ycn, ydn] = EllipticSNCNDN(z, k);
		const auto [ysn2, ycn2, ydn2] = EllipticSNCNDN(z / (1.0 + k2p), k2);
		REQUIRE(ysn == ApproxComplex((1.0 + k2p) * ysn2 * ycn2 / ydn2).margin(1e-15f));
		REQUIRE(ycn == ApproxComplex((1.0 + k2p) * (ydn2 * ydn2 - k2p) / (k2 * k2 * ydn2)).margin(1e-15f));
		REQUIRE(ydn == ApproxComplex((1.0 - k2p) * (ydn2 * ydn2 + k2p) / (k2 * k2 * ydn2)).margin(1e-15f));
	}
}