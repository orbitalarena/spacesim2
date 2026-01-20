#include "model/model.hpp"
#include "core/environment.hpp"
#include "core/vector.hpp"
#include "physics/orbit.hpp"
#include <iostream>
#include <cmath>

static std::array<double,3> pos(const Body& b){ return {b.x,b.y,b.z}; }

void run_model(PhysicsEngine& e,double dt,double t_end,OutputWriter* ow){
    auto& env=Environment::instance();

    for(double t=0;t<=t_end;t+=dt){
        e.step(dt);

        if(ow && ow->enabled()) ow->tick(t,e);

        if(std::fmod(t,3600.0)==0.0){
            auto sat=pos(e.bodies[0]);
            for(const auto& kv : env.all()){
                const auto& name = kv.first;
                const auto& obj  = kv.second;
                auto v=vec_between(sat,pos(obj.body));
                std::cout<<"t "<<t<<" target "<<name
                         <<" dist_km "<<vec_mag(v)
                         <<" angle_deg "<<vec_angle_deg(sat,v)<<"\n";
            }
        }
    }

    auto coe=body_to_coe_eci(e.bodies[0]);
    std::cout<<"final_coe a_km "<<coe.a<<" e "<<coe.e<<" i_rad "<<coe.i
             <<" raan_rad "<<coe.raan<<" argp_rad "<<coe.argp<<" ta_rad "<<coe.ta<<"\n";
}
