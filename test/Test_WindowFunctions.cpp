#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/DSP/WindowFunctions.hpp>

#include <Catch2/catch.hpp>
#include <algorithm>

using namespace enl;


TEST_CASE("Hamming window", "[AudioFramework:WindowFunctions]") {
	auto window = HammingWindow<float>(1024);

	REQUIRE(window.Length() == 1024);
	
	float max = *std::max_element(window.begin(), window.end());
	auto firstMinPosition = std::min_element(window.begin(), window.begin() + 512);
	auto lastMinPosition = std::min_element(window.begin() + 512, window.end());

	REQUIRE(Approx(max) == 1.0f);
	REQUIRE(firstMinPosition == window.begin());
	REQUIRE(lastMinPosition == window.end() - 1);
}


TEST_CASE("Hamming window complex", "[AudioFramework:WindowFunctions]") {
	auto window = HammingWindow<std::complex<float>>(1024);

	REQUIRE(window.Length() == 1024);

	
	REQUIRE(false);
}


TEST_CASE("Kaiser window", "[AudioFramework:WindowFunctions]") {
	auto window = KaiserWindow<float>(1024, 2.55f);

	REQUIRE(window.Length() == 1024);

	float max = *std::max_element(window.begin(), window.end());
	auto firstMinPosition = std::min_element(window.begin(), window.begin() + 512);
	auto lastMinPosition = std::min_element(window.begin() + 512, window.end());

	REQUIRE(Approx(max) == 1.0f);
	REQUIRE(firstMinPosition == window.begin());
	REQUIRE(lastMinPosition == window.end() - 1);	
}