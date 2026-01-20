#pragma once
#include <array>
std::array<double,3> sun_eci_km(double jd);
double angle_deg_between(const std::array<double,3>& a,const std::array<double,3>& b);
