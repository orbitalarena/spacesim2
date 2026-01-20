#include "model/rocket_model.hpp"
#include "physics/rocket.hpp"
#include <iostream>
#include <cmath>
#include <array>
#include <algorithm>

static constexpr double MU_E_KM3_S2 = 398600.4418;
static constexpr double R_E_KM      = 6371.0;

static double mag3(const std::array<double,3>& v){
    return std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

static std::array<double,3> sub3(const std::array<double,3>& a,const std::array<double,3>& b){
    return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};
}

int main_dummy_to_keep_compiler_happy();

void run_rocket_model(PhysicsEngine& e,double dt,double t_end,OutputWriter* ow){
    if(e.bodies.size() < 2) return;

    auto& ace = e.bodies[0];
    auto& rb  = e.bodies[1];

    std::vector<Stage> stages{
        {7.6e6,263,395000,25600},
        {9.34e5,421,92670,4000},
        {9.34e5,450,15000,3500}
    };

    Rocket r(stages);

    RocketState s0{};
    s0.x = rb.x; s0.y = rb.y; s0.z = rb.z;
    s0.vx = rb.vx; s0.vy = rb.vy; s0.vz = rb.vz;
    s0.mass = stages[0].fuel_kg + stages[0].dry_kg;
    r.set_state(s0);

    double prev_range = -1.0;

    for(double t=0;t<=t_end;t+=dt){
        ace.x += ace.vx*dt;
        ace.y += ace.vy*dt;
        ace.z += ace.vz*dt;

        auto rs = r.state();
        std::array<double,3> pR{rs.x,rs.y,rs.z};
        std::array<double,3> vR{rs.vx,rs.vy,rs.vz};
        std::array<double,3> pA{ace.x,ace.y,ace.z};
        std::array<double,3> vA{ace.vx,ace.vy,ace.vz};

        auto relp = sub3(pA,pR);
        auto relv = sub3(vA,vR);
        double range_km = mag3(relp);

        double rr_vec = 0.0;
        if(range_km > 1e-9){
            rr_vec = (relp[0]*relv[0] + relp[1]*relv[1] + relp[2]*relv[2]) / range_km;
        }

        double closing = std::max(0.2, std::fabs(rr_vec));
        double tau = std::clamp(range_km / closing, 0.0, 7200.0);

        double target_xyz[3] = {
            ace.x + ace.vx*tau,
            ace.y + ace.vy*tau,
            ace.z + ace.vz*tau
        };

        r.step(dt, MU_E_KM3_S2, target_xyz);

        rs = r.state();
        rb.x=rs.x; rb.y=rs.y; rb.z=rs.z;
        rb.vx=rs.vx; rb.vy=rs.vy; rb.vz=rs.vz;
        rb.mass=rs.mass;

        if(ow && ow->enabled()) ow->tick(t,e);

        double rr = std::sqrt(rs.x*rs.x+rs.y*rs.y+rs.z*rs.z);
        double alt_km = rr - R_E_KM;

        pR = {rs.x,rs.y,rs.z};
        vR = {rs.vx,rs.vy,rs.vz};
        pA = {ace.x,ace.y,ace.z};
        vA = {ace.vx,ace.vy,ace.vz};

        relp = sub3(pA,pR);
        relv = sub3(vA,vR);
        range_km = mag3(relp);

        rr_vec = 0.0;
        if(range_km > 1e-9){
            rr_vec = (relp[0]*relv[0] + relp[1]*relv[1] + relp[2]*relv[2]) / range_km;
        }

        double rr_num = 0.0;
        if(prev_range >= 0.0) rr_num = (range_km - prev_range) / dt;
        prev_range = range_km;

        if(std::fmod(t,60.0)==0.0){
            std::cout<<"t "<<t<<" rocket_alt_km "<<alt_km<<"\n";
            std::cout<<"t "<<t<<" rocket_to_ace_range_km "<<range_km
                     <<" rr_vec_km_s "<<rr_vec
                     <<" rr_num_km_s "<<rr_num<<"\n";
        }

        if(r.stage_sep()){
            std::cout<<"t "<<t<<" rocket_stage_sep\n";
        }

        if(alt_km <= 0.0 && t > 10.0){
            std::cout<<"IMPACT t "<<t<<" rocket_alt_km "<<alt_km<<"\n";
            break;
        }

        if(range_km < 1.0){
            std::cout<<"ARRIVAL t "<<t
                     <<" range_km "<<range_km
                     <<" rr_vec_km_s "<<rr_vec
                     <<" rr_num_km_s "<<rr_num<<"\n";
            break;
        }
    }
}
