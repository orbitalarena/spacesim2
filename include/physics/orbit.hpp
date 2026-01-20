#pragma once
#include "physics/engine.hpp"
struct COE { double a,e,i,raan,argp,ta; }; // a km, angles rad
Body coe_to_body_eci(const COE&);
COE body_to_coe_eci(const Body&);
