#pragma once
#include <array>
std::array<double,3> lla_to_ecef(double lat_deg,double lon_deg,double alt_km);
double local_solar_time_hours(double lon_deg,double jd);
double ecef_lon_deg(double x,double y,double z);
