#include "physics/sun.hpp"
#include <cmath>

std::array<double,3> sun_eci(double jd){
    double n=jd-2451545.0;
    double L=fmod(280.46+0.9856474*n,360.0);
    double g=fmod(357.528+0.9856003*n,360.0);
    double lambda=L+1.915*sin(g*M_PI/180)+0.02*sin(2*g*M_PI/180);
    return {cos(lambda*M_PI/180),sin(lambda*M_PI/180),0};
}

double sun_angle_deg(const std::array<double,3>& sat,const std::array<double,3>& sun){
    double dot=sat[0]*sun[0]+sat[1]*sun[1]+sat[2]*sun[2];
    double ms=sqrt(sat[0]*sat[0]+sat[1]*sat[1]+sat[2]*sat[2]);
    return acos(dot/ms)*180/M_PI;
}
