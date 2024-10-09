#pragma once
#include <memory>


namespace dspbb {

enum class eSignalDomain {
	TIME,
	FREQUENCY,
	QUEFRENCY,
	DOMAINLESS,
};
static constexpr auto TIME_DOMAIN = eSignalDomain::TIME;
static constexpr auto FREQUENCY_DOMAIN = eSignalDomain::FREQUENCY;
static constexpr auto QUEFRENCY_DOMAIN = eSignalDomain::QUEFRENCY;
static constexpr auto DOMAINLESS = eSignalDomain::DOMAINLESS;


template <class T, eSignalDomain Domain>
class BasicSignal;


template <class T, eSignalDomain Domain>
class BasicSignalView;

} // namespace dspbb