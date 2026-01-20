#include "model/rocket_model.hpp"
#include "physics/orbit.hpp"
#include "core/vector.hpp"
#include <iostream>
#include <cmath>

static constexpr double MU_E=398600.4418e9;

void run_rocket_model(PhysicsEngine& e,double dt,double t_end){
    std::vector<Stage> stages{
        {7.6e6,263,395000,25600},
        {9.34e5,421,92670,4000},
        {9.34e5,450,15000,3500}
    };

    Rocket r(stages);
    auto& ace=e.bodies[0];

    for(double t=0;t<=t_end;t+=dt){
        r.step(dt,MU_E);
        auto s=r.state();

        double alt=sqrt(s.x*s.x+s.y*s.y+s.z*s.z)-6371000.0;
        if(fmod(t,60.0)==0.0)
            std::cout<<"t "<<t<<" rocket_alt_m "<<alt<<"\n";

        if(r.stage_sep())
            std::cout<<"t "<<t<<" rocket_stage_sep\n";

        auto d=vec_between({s.x,s.y,s.z},{ace.x,ace.y,ace.z});
        double range=vec_mag(d);

        if(range<1000.0){
            double rr=(d[0]*s.vx+d[1]*s.vy+d[2]*s.vz)/range;
            std::cout<<"ARRIVAL t "<<t
                     <<" range_m "<<range
                     <<" range_rate "<<rr<<"\n";
            break;
        }
    }
}
