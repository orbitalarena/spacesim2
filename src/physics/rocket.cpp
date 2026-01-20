#include "physics/rocket.hpp"
#include <cmath>

static constexpr double G0=9.80665;
static constexpr double R_E=6371000.0;

Rocket::Rocket(const std::vector<Stage>& st)
: stages(st), active_stage(0), sep_flag(false)
{
    mass=0.0;
    for(auto&s:st) mass+=s.fuel_mass+s.dry_mass;
    s={R_E,0,0,0,0,0};
}

void Rocket::step(double dt,double mu){
    sep_flag=false;

    double r=sqrt(s.x*s.x+s.y*s.y+s.z*s.z);
    double a_grav=-mu/(r*r);

    if(active_stage<stages.size()){
        auto& st=stages[active_stage];
        double mdot=st.thrust/(st.isp*G0);
        double dm=mdot*dt;

        if(st.fuel_mass>0){
            double use=std::min(st.fuel_mass,dm);
            st.fuel_mass-=use;
            mass-=use;

            double a_thrust=st.thrust/mass;
            s.vx+=(a_thrust+a_grav)*dt;
        } else {
            mass-=st.dry_mass;
            active_stage++;
            sep_flag=true;
        }
    } else {
        s.vx+=a_grav*dt;
    }

    s.x+=s.vx*dt;
}

RocketState Rocket::state() const { return s; }
bool Rocket::stage_sep() const { return sep_flag; }
