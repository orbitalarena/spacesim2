#include "physics/rocket.hpp"
#include <cmath>
#include <algorithm>

static constexpr double G0 = 9.80665; // m/s^2

Rocket::Rocket(const std::vector<Stage>& stages): st(stages){
    if(!st.empty()){
        cur=0;
        fuel = st[0].fuel_kg;
        dry  = st[0].dry_kg;
        mass = dry + fuel;
    }
}

void Rocket::set_state(const RocketState& s){
    x=s.x; y=s.y; z=s.z;
    vx=s.vx; vy=s.vy; vz=s.vz;
    mass=s.mass;
}

RocketState Rocket::state() const{
    return {x,y,z,vx,vy,vz,mass};
}

double Rocket::mdot() const{
    if(cur>=st.size()) return 0.0;
    if(st[cur].isp_s<=0.0) return 0.0;
    return st[cur].thrust_n/(st[cur].isp_s*G0); // kg/s
}

bool Rocket::stage_sep(){
    bool r=sep;
    sep=false;
    return r;
}

bool Rocket::alive() const{
    return (cur < st.size()) && (mass > 0.0);
}

void Rocket::step(double dt,double mu_km3_s2,const std::array<double,3>& dir_in){
    // gravity (km, km/s, mu in km^3/s^2)
    double r = std::sqrt(x*x+y*y+z*z);
    double ax=0,ay=0,az=0;
    if(r>1e-9){
        double invr3 = 1.0/(r*r*r);
        ax -= mu_km3_s2*x*invr3;
        ay -= mu_km3_s2*y*invr3;
        az -= mu_km3_s2*z*invr3;
    }

    // thrust: convert m/s^2 -> km/s^2 by /1000
    if(cur < st.size()){
        auto& S = st[cur];
        if(fuel > 0.0 && S.thrust_n > 0.0 && mass > 1.0){
            double dx=dir_in[0], dy=dir_in[1], dz=dir_in[2];
            double dm = std::sqrt(dx*dx+dy*dy+dz*dz);
            if(dm < 1e-12){ dx=1; dy=0; dz=0; dm=1; }
            dx/=dm; dy/=dm; dz/=dm;

            double a_th_km = (S.thrust_n / mass) / 1000.0; // km/s^2
            ax += a_th_km*dx;
            ay += a_th_km*dy;
            az += a_th_km*dz;

            double use = mdot()*dt;
            if(use > fuel) use = fuel;
            fuel -= use;
            mass -= use;

            if(fuel <= 0.0){
                mass -= dry;
                sep=true;
                cur++;
                if(cur < st.size()){
                    fuel = st[cur].fuel_kg;
                    dry  = st[cur].dry_kg;
                    mass += (fuel + dry);
                }
            }
        }
    }

    // integrate (semi-implicit Euler) in km-space
    vx += ax*dt; vy += ay*dt; vz += az*dt;
    x  += vx*dt; y  += vy*dt; z  += vz*dt;
}
