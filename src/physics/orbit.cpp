#include "physics/orbit.hpp"
#include <cmath>

static constexpr double MU=398600.4418;

Body coe_to_body(const COE& c){
    double r=c.a*(1-c.e*c.e)/(1+c.e*cos(c.ta));
    Body b{};
    b.x=r*cos(c.ta);
    b.y=r*sin(c.ta);
    b.z=0;
    double v=sqrt(MU/c.a);
    b.vx=-v*sin(c.ta);
    b.vy=v*cos(c.ta);
    b.vz=0;
    b.mass=1000;
    return b;
}

COE body_to_coe(const Body& b){
    COE c{};
    double r=sqrt(b.x*b.x+b.y*b.y+b.z*b.z);
    c.a=r;
    c.e=0;
    c.i=0;
    c.raan=0;
    c.argp=0;
    c.ta=atan2(b.y,b.x);
    return c;
}
