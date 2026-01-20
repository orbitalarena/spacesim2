#pragma once
#include "physics/engine.hpp"
struct COE { double a,e,i,raan,argp,ta; };
Body coe_to_body(const COE&);
COE body_to_coe(const Body&);
