#include "orbit/lambert.hpp"
#include <cmath>
#include <algorithm>

static double C(double z){
    if(z>0) return (1-std::cos(std::sqrt(z)))/z;
    if(z<0) return (1-std::cosh(std::sqrt(-z)))/z;
    return 0.5;
}
static double S(double z){
    if(z>0) return (std::sqrt(z)-std::sin(std::sqrt(z)))/std::pow(z,1.5);
    if(z<0) return (std::sinh(std::sqrt(-z))-std::sqrt(-z))/std::pow(-z,1.5);
    return 1.0/6.0;
}

static LambertSolution solve(const std::array<double,3>& r1,
                             const std::array<double,3>& r2,
                             double tof,double mu){
    auto dot=[](auto&a,auto&b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];};
    auto nrm=[&](auto&a){return std::sqrt(dot(a,a));};

    double r1m=nrm(r1), r2m=nrm(r2);
    double dnu=std::acos(std::clamp(dot(r1,r2)/(r1m*r2m),-1.0,1.0));
    double A=std::sin(dnu)*std::sqrt(r1m*r2m/(1-std::cos(dnu)));

    double z=0;
    for(int i=0;i<40;i++){
        double y=r1m+r2m+A*(z*S(z)-1)/std::sqrt(C(z));
        double F=std::pow(y/C(z),1.5)*S(z)+A*std::sqrt(y)-std::sqrt(mu)*tof;
        z-=F/1000.0;
    }

    double y=r1m+r2m+A*(z*S(z)-1)/std::sqrt(C(z));
    double f=1-y/r1m;
    double g=A*std::sqrt(y/mu);
    double gdot=1-y/r2m;

    LambertSolution sol{};
    for(int i=0;i<3;i++){
        sol.v1_km_s[i]=(r2[i]-f*r1[i])/g;
        sol.v2_km_s[i]=(gdot*r2[i]-r1[i])/g;
    }
    return sol;
}

std::vector<LambertSolution>
lambert_intercept(const std::array<double,3>& r1,
                  const std::array<double,3>& r2,
                  double tof,double mu,int){
    return { solve(r1,r2,tof,mu) };
}

std::vector<LambertSolution>
lambert_rendezvous(const std::array<double,3>& r1,
                   const std::array<double,3>& v2,
                   const std::array<double,3>& r2,
                   double tof,double mu,int m){
    auto s=lambert_intercept(r1,r2,tof,mu,m);
    s[0].v2_km_s=v2;
    return s;
}

std::vector<LambertSolution>
lambert_nmc(const std::array<double,3>& r1,
            const std::array<double,3>&,
            const std::array<double,3>& r2,
            const std::array<double,3>& v2,
            double tof,const NMCParams& nmc,double mu){
    auto s=lambert_rendezvous(r1,v2,r2,tof,mu,0);
    double a=nmc.entry_angle_deg*M_PI/180.0;
    s[0].v2_km_s[0]+=nmc.semi_minor_km/1000*std::cos(a);
    s[0].v2_km_s[1]+=nmc.semi_minor_km/1000*std::sin(a);
    return s;
}
