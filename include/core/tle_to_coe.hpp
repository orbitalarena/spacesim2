#pragma once
#include "core/tle.hpp"
#include "physics/orbit.hpp"   // authoritative COE definition

// Fill existing COE struct from a TLE (km + radians)
void tle_to_coe(const TLE& t, double mu_km3_s2, COE& out);
