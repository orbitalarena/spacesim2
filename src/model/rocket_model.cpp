#include "model/rocket_model.hpp"
#include "physics/rocket.hpp"
#include "core/vector.hpp"
#include <iostream>
#include <cmath>
#include <limits>

static constexpr double MU_E_KM3_S2 = 398600.4418; // km^3/s^2
static constexpr double R_E_KM      = 6371.0;

static int find_idx(const PhysicsEngine& e,const std::string& name){
    for(size_t i=0;i<e.names.size();++i) if(e.names[i]==name) return (int)i;
    return -1;
}

static std::array<double,3> unit_to(const Body& from,const Body& to){
    double dx=to.x-from.x, dy=to.y-from.y, dz=to.z-from.z;
    double m=std::sqrt(dx*dx+dy*dy+dz*dz);
    if(m<1e-12) return {1,0,0};
    return {dx/m,dy/m,dz/m};
}

void run_rocket_model(PhysicsEngine& e,double dt,double t_end,OutputWriter* ow){
    int ace_i = find_idx(e,"Ace");
    int rok_i = find_idx(e,"Rocket");
    if(ace_i<0 || rok_i<0){
        std::cout<<"rocket_model_missing_entities\n";
        return;
    }

    // keep stage params (thrust N, isp s, masses kg)
    std::vector<Stage> stages{
        {7.6e6,263,395000,25600},
        {9.34e5,421,92670,4000},
        {9.34e5,450,15000,3500}
    };

    Rocket r(stages);

    // init rocket from scenario body (km, km/s)
    {
        auto& b = e.bodies[rok_i];
        RocketState s{};
        s.x=b.x; s.y=b.y; s.z=b.z;
        s.vx=b.vx; s.vy=b.vy; s.vz=b.vz;
        double m0 = stages[0].fuel_kg + stages[0].dry_kg;
        s.mass = (b.mass>1.0)? b.mass : m0;
        r.set_state(s);
        b.mass = s.mass;
    }

    double prev_range = std::numeric_limits<double>::quiet_NaN();

    for(double t=0;t<=t_end;t+=dt){
        // propagate other entities (Ace)
        e.step(dt);

        const auto& ace = e.bodies[ace_i];
        auto& rok = e.bodies[rok_i];

        auto dir = unit_to(rok,ace);

        r.step(dt, MU_E_KM3_S2, dir);
        auto s = r.state();

        // write back so --output sees rocket
        rok.x=s.x; rok.y=s.y; rok.z=s.z;
        rok.vx=s.vx; rok.vy=s.vy; rok.vz=s.vz;
        rok.mass=s.mass;

        if(ow && ow->enabled()) ow->tick(t,e);

        double rmag = std::sqrt(s.x*s.x+s.y*s.y+s.z*s.z);
        double alt_km = rmag - R_E_KM;

        if(std::fmod(t,60.0)==0.0){
            std::cout<<"t "<<t<<" rocket_alt_km "<<alt_km<<"\n";
        }

        if(r.stage_sep()){
            std::cout<<"t "<<t<<" rocket_stage_sep\n";
        }

        // range + range-rate to Ace (km, km/s)
        std::array<double,3> pR{s.x,s.y,s.z};
        std::array<double,3> vR{s.vx,s.vy,s.vz};
        std::array<double,3> pA{ace.x,ace.y,ace.z};
        std::array<double,3> vA{ace.vx,ace.vy,ace.vz};

        auto rel_p = vec_between(pR,pA);
        auto rel_v = vec_between(vR,vA);
        double range_km = vec_mag(rel_p);

        if(std::fmod(t,60.0)==0.0){
            double rr_vec = (rel_p[0]*rel_v[0] + rel_p[1]*rel_v[1] + rel_p[2]*rel_v[2]) / std::max(1e-12,range_km);
            double rr_num = std::isnan(prev_range) ? rr_vec : (range_km - prev_range)/dt;
            std::cout<<"t "<<t<<" rocket_to_ace_range_km "<<range_km<<" rr_vec_km_s "<<rr_vec<<" rr_num_km_s "<<rr_num<<"\n";
            prev_range = range_km;
        }

        if(range_km < 1.0){
            double rr_vec = (rel_p[0]*rel_v[0] + rel_p[1]*rel_v[1] + rel_p[2]*rel_v[2]) / std::max(1e-12,range_km);
            double rr_num = std::isnan(prev_range) ? rr_vec : (range_km - prev_range)/dt;
            std::cout<<"ARRIVAL t "<<t<<" range_km "<<range_km<<" rr_vec_km_s "<<rr_vec<<" rr_num_km_s "<<rr_num<<"\n";
            break;
        }

        if(!r.alive()){
            std::cout<<"ROCKET_DEAD t "<<t<<"\n";
            break;
        }
    }
}
