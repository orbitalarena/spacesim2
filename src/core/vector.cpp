#include "core/vector.hpp"
#include <cmath>

std::array<double,3> vec_between(const std::array<double,3>& a,const std::array<double,3>& b){
    return {b[0]-a[0],b[1]-a[1],b[2]-a[2]};
}

double vec_mag(const std::array<double,3>& v){
    return std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

double vec_angle_deg(const std::array<double,3>& a,const std::array<double,3>& b){
    double dot=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    double ma=vec_mag(a), mb=vec_mag(b);
    return std::acos(dot/(ma*mb))*180.0/M_PI;
}
