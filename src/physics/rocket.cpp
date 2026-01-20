#include "physics/rocket.hpp"
#include "physics/atmosphere.hpp"
#include <cmath>

static constexpr double G0=9.80665;

Rocket::Rocket(const std::vector<Stage>& st):stages(st),cur(0),sep(false){
    s={};
    s.m=0;
    for(auto&x:stages) s.m+=x.fuel+x.dry;
}

bool Rocket::step(double dt,double mu){
    sep=false;
    if(cur>=stages.size()) return false;
    auto& st=stages[cur];
    double r=sqrt(s.x*s.x+s.y*s.y+s.z*s.z);
    double alt=r-6371000.0;
    double rho=air_density(alt);
    double v=sqrt(s.vx*s.vx+s.vy*s.vy+s.vz*s.vz);
    double drag=0.5*rho*v*v*5.0;

    double mdot=st.thrust/(st.isp*G0);
    double dm=mdot*dt;
    if(dm>st.fuel) dm=st.fuel;
    st.fuel-=dm;
    s.m-=dm;

    double a_thrust=st.thrust/s.m;
    double ax=a_thrust*(s.x/r);
    double ay=a_thrust*(s.y/r);
    double az=a_thrust*(s.z/r);

    double ag=-mu/(r*r);
    ax+=ag*(s.x/r);
    ay+=ag*(s.y/r);
    az+=ag*(s.z/r);

    if(v>0){
        ax-=drag*(s.vx/v)/s.m;
        ay-=drag*(s.vy/v)/s.m;
        az-=drag*(s.vz/v)/s.m;
    }

    s.vx+=ax*dt;
    s.vy+=ay*dt;
    s.vz+=az*dt;
    s.x+=s.vx*dt;
    s.y+=s.vy*dt;
    s.z+=s.vz*dt;

    if(st.fuel<=0){
        s.m-=st.dry;
        cur++;
        sep=true;
    }
    return true;
}

RocketState Rocket::state() const{ return s; }
bool Rocket::stage_sep() const{ return sep; }
