#include "model/model.hpp"
#include "core/environment.hpp"
#include "core/vector.hpp"
#include "core/geodesy.hpp"
#include "physics/orbit.hpp"
#include <iostream>
#include <cmath>

static std::array<double,3> pos(const Body& b){ return {b.x,b.y,b.z}; }
static std::array<double,3> norm(const std::array<double,3>& v){
    double m=vec_mag(v);
    if(m==0) return {0,0,0};
    return {v[0]/m,v[1]/m,v[2]/m};
}

void run_model(PhysicsEngine& e,double dt,double t_end){
    auto& env=Environment::instance();
    double jd=2451545.0;

    double last_hour_mark=0.0;

    for(double t=0.0;t<=t_end;t+=dt){
        e.step(dt);
        jd+=dt/86400.0;

        if(t - last_hour_mark >= 3600.0){
            last_hour_mark = std::floor(t/3600.0)*3600.0;

            auto ace=pos(e.bodies[0]);
            auto leo=pos(e.bodies[1]);

            auto dcB=env.get("WashingtonDC").body;
            auto laB=env.get("LosAngeles").body;

            auto vdc=vec_between(ace,pos(dcB));
            auto vla=vec_between(ace,pos(laB));

            auto radial = norm(ace);
            auto vdc_hat=norm(vdc);
            auto vla_hat=norm(vla);

            double dotdc=radial[0]*vdc_hat[0]+radial[1]*vdc_hat[1]+radial[2]*vdc_hat[2];
            if(dotdc>1) dotdc=1; if(dotdc<-1) dotdc=-1;

            double dotla=radial[0]*vla_hat[0]+radial[1]*vla_hat[1]+radial[2]*vla_hat[2];
            if(dotla>1) dotla=1; if(dotla<-1) dotla=-1;

            std::cout<<"t "<<last_hour_mark
                     <<" Ace_DC_deg "<<(std::acos(dotdc)*180.0/M_PI)<<" Ace_DC_km "<<vec_mag(vdc)
                     <<" Ace_LA_deg "<<(std::acos(dotla)*180.0/M_PI)<<" Ace_LA_km "<<vec_mag(vla)
                     <<"\n";

            double best=1e300; std::string best_city;
            double best_lon=0.0;

            for(auto& [name,obj]:env.all()){
                if(obj.body.mass!=0) continue;
                double d=vec_mag(vec_between(leo,pos(obj.body)));
                if(d<best){
                    best=d;
                    best_city=name;
                    best_lon=ecef_lon_deg(obj.body.x,obj.body.y,obj.body.z);
                }
            }

            double lst=local_solar_time_hours(best_lon,jd);
            std::cout<<"t "<<last_hour_mark
                     <<" LeoClosest "<<best_city
                     <<" dist_km "<<best
                     <<" LST "<<lst
                     <<"\n";
        }
    }

    auto coe=body_to_coe_eci(e.bodies[0]);
    std::cout<<"final_coe a_km "<<coe.a<<" e "<<coe.e<<" i_rad "<<coe.i
             <<" raan_rad "<<coe.raan<<" argp_rad "<<coe.argp<<" ta_rad "<<coe.ta
             <<"\n";
}
