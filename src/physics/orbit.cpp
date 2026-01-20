#include "physics/orbit.hpp"
#include <cmath>

static constexpr double MU=398600.4418;

static void rot3(double a,double& x,double& y){
    double c=std::cos(a), s=std::sin(a);
    double X=c*x - s*y;
    double Y=s*x + c*y;
    x=X; y=Y;
}
static void rot1(double a,double& y,double& z){
    double c=std::cos(a), s=std::sin(a);
    double Y=c*y - s*z;
    double Z=s*y + c*z;
    y=Y; z=Z;
}

Body coe_to_body_eci(const COE& c){
    double p=c.a*(1.0-c.e*c.e);
    double r=p/(1.0+c.e*std::cos(c.ta));
    double x=r*std::cos(c.ta);
    double y=r*std::sin(c.ta);
    double z=0.0;

    double h=std::sqrt(MU*p);
    double vx=-MU/h*std::sin(c.ta);
    double vy= MU/h*(c.e+std::cos(c.ta));
    double vz=0.0;

    rot3(c.argp,x,y);
    rot3(c.argp,vx,vy);

    rot1(c.i,y,z);
    rot1(c.i,vy,vz);

    rot3(c.raan,x,y);
    rot3(c.raan,vx,vy);

    Body b{};
    b.x=x; b.y=y; b.z=z;
    b.vx=vx; b.vy=vy; b.vz=vz;
    b.mass=1000.0;
    return b;
}

COE body_to_coe_eci(const Body& b){
    double rx=b.x, ry=b.y, rz=b.z;
    double vx=b.vx, vy=b.vy, vz=b.vz;

    double r=std::sqrt(rx*rx+ry*ry+rz*rz);
    double v=std::sqrt(vx*vx+vy*vy+vz*vz);

    double hx=ry*vz - rz*vy;
    double hy=rz*vx - rx*vz;
    double hz=rx*vy - ry*vx;
    double h=std::sqrt(hx*hx+hy*hy+hz*hz);

    double nx=-hy;
    double ny= hx;
    double n=std::sqrt(nx*nx+ny*ny);

    double ex=( (v*v - MU/r)*rx - (rx*vx+ry*vy+rz*vz)*vx )/MU;
    double ey=( (v*v - MU/r)*ry - (rx*vx+ry*vy+rz*vz)*vy )/MU;
    double ez=( (v*v - MU/r)*rz - (rx*vx+ry*vy+rz*vz)*vz )/MU;
    double e=std::sqrt(ex*ex+ey*ey+ez*ez);

    double i=std::acos(hz/h);

    double raan=0.0;
    if(n>1e-12){
        raan=std::acos(nx/n);
        if(ny<0) raan=2*M_PI-raan;
    }

    double argp=0.0;
    if(n>1e-12 && e>1e-12){
        argp=std::acos((nx*ex+ny*ey)/(n*e));
        if(ez<0) argp=2*M_PI-argp;
    }

    double ta=0.0;
    if(e>1e-12){
        ta=std::acos((ex*rx+ey*ry+ez*rz)/(e*r));
        if(rx*vx+ry*vy+rz*vz<0) ta=2*M_PI-ta;
    } else {
        ta=std::atan2(ry,rx);
        if(ta<0) ta+=2*M_PI;
    }

    double a=1.0/(2.0/r - (v*v)/MU);

    COE c{};
    c.a=a; c.e=e; c.i=i; c.raan=raan; c.argp=argp; c.ta=ta;
    return c;
}
