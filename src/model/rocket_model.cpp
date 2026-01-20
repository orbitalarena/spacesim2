#include "model/rocket_model.hpp"
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

    double prev_range=-1.0;

    for(double t=0;t<=t_end;t+=dt){
        r.step(dt,MU_E);
        auto s=r.state();

        double alt=std::sqrt(s.x*s.x+s.y*s.y+s.z*s.z)-6371000.0;

        std::array<double,3> pR{s.x,s.y,s.z};
        std::array<double,3> vR{s.vx,s.vy,s.vz};
        std::array<double,3> pA{ace.x,ace.y,ace.z};
        std::array<double,3> vA{ace.vx,ace.vy,ace.vz};

        auto rel_p = vec_between(pR,pA);
        double range = vec_mag(rel_p);

        auto rel_v = vec_between(vR,vA);
        double rr_vec = (range>0) ? (rel_p[0]*rel_v[0]+rel_p[1]*rel_v[1]+rel_p[2]*rel_v[2])/range : 0.0;
        double rr_num = (prev_range<0) ? 0.0 : (range-prev_range)/dt;
        prev_range=range;

        if(std::fmod(t,60.0)==0.0){
            std::cout<<"t "<<t<<" rocket_alt_m "<<alt
                     <<" ace_range_m "<<range
                     <<" ace_rr_vec "<<rr_vec
                     <<" ace_rr_num "<<rr_num
                     <<"\n";
        }

        if(r.stage_sep())
            std::cout<<"t "<<t<<" rocket_stage_sep\n";

        if(range<1000.0){
            std::cout<<"ARRIVAL t "<<t<<" range_m "<<range
                     <<" ace_rr_vec "<<rr_vec
                     <<" ace_rr_num "<<rr_num
                     <<"\n";
            break;
        }
    }
}
