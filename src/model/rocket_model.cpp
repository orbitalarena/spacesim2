#include "physics/rocket.hpp"
#include "physics/orbit.hpp"
#include "core/vector.hpp"
#include <iostream>
#include <cmath>

void run_rocket(Rocket& r,const Body& ace,double mu,double dt,double t_end){
    double t=0;
    for(;t<=t_end;t+=dt){
        r.step(dt,mu);
        auto s=r.state();
        double alt=sqrt(s.x*s.x+s.y*s.y+s.z*s.z)-6371000.0;
        if(fmod(t,60)==0)
            std::cout<<"t "<<t<<" alt_m "<<alt<<"\n";
        if(r.stage_sep())
            std::cout<<"t "<<t<<" stage_sep\n";
        auto d=vec_between({s.x,s.y,s.z},{ace.x,ace.y,ace.z});
        double range=vec_mag(d);
        if(range<1000){
            double rr=(d[0]*s.vx+d[1]*s.vy+d[2]*s.vz)/range;
            std::cout<<"ARRIVAL t "<<t<<" range_rate "<<rr<<"\n";
            break;
        }
    }
}
