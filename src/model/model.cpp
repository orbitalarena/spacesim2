#include "model/model.hpp"
#include "sim/scenario.hpp"
#include "physics/orbit.hpp"
#include "physics/sun.hpp"
#include <fstream>
#include <array>
#include <cmath>

static std::array<double,3> v3(double x,double y,double z){ return {x,y,z}; }

static std::array<double,3> norm(const std::array<double,3>& a){
    double m=std::sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    return {a[0]/m,a[1]/m,a[2]/m};
}

static std::array<double,3> sub(const std::array<double,3>& a,const std::array<double,3>& b){
    return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};
}

void run_model(PhysicsEngine& e,double dt,double t_end){
    std::ofstream out("models/ace_geo.model");
    double jd=2451545.0;
    for(double t=0;t<=t_end;t+=dt){
        e.step(dt);
        jd+=dt/86400.0;

        if(std::fmod(t,3600.0)==0.0){
            auto& b=e.bodies[0];
            auto sat=v3(b.x,b.y,b.z);
            auto sun=sun_eci_km(jd);

            auto sun_from_sat = sub(sun,sat);
            auto zenith = sat;

            double zenith_angle = angle_deg_between(norm(zenith), norm(sun_from_sat));
            double nadir_angle  = angle_deg_between(norm({-zenith[0],-zenith[1],-zenith[2]}), norm(sun_from_sat));

            out<<"t "<<t<<" sun_zenith_angle_deg "<<zenith_angle<<" sun_nadir_angle_deg "<<nadir_angle<<"\n";
        }
    }

    auto& b=e.bodies[0];
    auto coe=body_to_coe_eci(b);

    out<<"final_pos_km "<<b.x<<" "<<b.y<<" "<<b.z<<"\n";
    out<<"final_vel_kmps "<<b.vx<<" "<<b.vy<<" "<<b.vz<<"\n";
    out<<"final_coe a_km "<<coe.a<<" e "<<coe.e<<" i_rad "<<coe.i<<" raan_rad "<<coe.raan<<" argp_rad "<<coe.argp<<" ta_rad "<<coe.ta<<"\n";
}
