#pragma once
#include <array>
std::array<double,3> sun_eci(double jd);
double sun_angle_deg(const std::array<double,3>& sat,const std::array<double,3>& sun);
