#pragma once
#include <array>
std::array<double,3> vec_between(const std::array<double,3>& a,const std::array<double,3>& b);
double vec_mag(const std::array<double,3>&);
double vec_angle_deg(const std::array<double,3>&,const std::array<double,3>&);
