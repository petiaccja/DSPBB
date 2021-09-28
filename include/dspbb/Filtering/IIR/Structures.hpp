#pragma once


#include "../../Primitives/Signal.hpp"


namespace dspbb {


template <class T>
class DirectForm1 {
public:
	Signal<T, eSignalDomain::DOMAINLESS> poleState;
	Signal<T, eSignalDomain::DOMAINLESS> zeroState;
};

template <class T>
class DirectForm2 {
public:
	Signal<T, eSignalDomain::DOMAINLESS> state;
};

} // namespace dspbb