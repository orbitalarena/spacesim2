#include "physics/sun.hpp"
#include <cmath>

static constexpr double AU_KM=149597870.7;
static constexpr double DEG=M_PI/180.0;

static double wrap2pi(double x){
    x=fmod(x,2*M_PI);
    if(x<0) x+=2*M_PI;
    return x;
}

std::array<double,3> sun_eci_km(double jd){
    double T=(jd-2451545.0)/36525.0;

    double L0=wrap2pi((280.46646 + 36000.76983*T + 0.0003032*T*T)*DEG);
    double M =wrap2pi((357.52911 + 35999.05029*T - 0.0001537*T*T)*DEG);

    double C=( (1.914602 - 0.004817*T - 0.000014*T*T)*sin(M)
             + (0.019993 - 0.000101*T)*sin(2*M)
             + 0.000289*sin(3*M) )*DEG;

    double theta=L0 + C;
    double Omega=wrap2pi((125.04 - 1934.136*T)*DEG);
    double lambda=theta - (0.00569*DEG) - (0.00478*DEG)*sin(Omega);

    double eps0=(23.439291 - 0.0130042*T - 1.64e-7*T*T + 5.04e-7*T*T*T)*DEG;
    double eps=eps0 + (0.00256*DEG)*cos(Omega);

    double R_au = 1.000001018*(1-0.016708634*cos(M)-0.000042037*cos(2*M)-0.0000001236*cos(3*M));

    double x = R_au*cos(lambda);
    double y = R_au*cos(eps)*sin(lambda);
    double z = R_au*sin(eps)*sin(lambda);

    return {x*AU_KM,y*AU_KM,z*AU_KM};
}

double angle_deg_between(const std::array<double,3>& a,const std::array<double,3>& b){
    double dot=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    double ma=std::sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    double mb=std::sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
    double c=dot/(ma*mb);
    if(c>1) c=1;
    if(c<-1) c=-1;
    return std::acos(c)*180.0/M_PI;
}
