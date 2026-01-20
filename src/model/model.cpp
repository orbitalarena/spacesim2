#include "model/model.hpp"
#include "physics/orbit.hpp"
#include "physics/sun.hpp"
#include <fstream>
#include <cmath>

void run_model(PhysicsEngine& e,double dt,double t_end){
    std::ofstream out("models/ace_geo.model");
    double jd=2451545.0;
    for(double t=0;t<=t_end;t+=dt){
        e.step(dt);
        jd+=dt/86400.0;
        auto& b=e.bodies[0];
        if(std::fmod(t,3600.0)==0){
            auto sun=sun_eci(jd);
            std::array<double,3> sat{b.x,b.y,b.z};
            out<<"t "<<t<<" sun_angle "<<sun_angle_deg(sat,sun)<<"\n";
        }
    }
    auto coe=body_to_coe(e.bodies[0]);
    auto& b=e.bodies[0];
    out<<"final_pos "<<b.x<<" "<<b.y<<" "<<b.z<<"\n";
    out<<"final_vel "<<b.vx<<" "<<b.vy<<" "<<b.vz<<"\n";
    out<<"final_coe "<<coe.a<<" "<<coe.e<<" "<<coe.i<<" "<<coe.raan<<" "<<coe.argp<<" "<<coe.ta<<"\n";
}
