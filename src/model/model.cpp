#include "model/model.hpp"
#include "core/environment.hpp"
#include "core/vector.hpp"
#include "core/geodesy.hpp"
#include <fstream>
#include <limits>
#include <cmath>

static std::array<double,3> pos(const Body& b){ return {b.x,b.y,b.z}; }

void run_model(PhysicsEngine& e,double dt,double t_end){
    std::ofstream out("models/ace_cities.model");
    auto& env=Environment::instance();
    double jd=2451545.0;

    for(double t=0;t<=t_end;t+=dt){
        e.step(dt);
        jd+=dt/86400.0;

        if(std::fmod(t,3600.0)==0){
            auto ace=pos(e.bodies[0]);
            auto leo=pos(e.bodies[1]);

            auto dc=env.get("WashingtonDC").body;
            auto la=env.get("LosAngeles").body;

            out<<"t "<<t
               <<" Ace_DC_deg "<<vec_angle_deg(ace,vec_between(ace,pos(dc)))
               <<" Ace_DC_km "<<vec_mag(vec_between(ace,pos(dc)))
               <<" Ace_LA_deg "<<vec_angle_deg(ace,vec_between(ace,pos(la)))
               <<" Ace_LA_km "<<vec_mag(vec_between(ace,pos(la)))<<"\n";

            double best=1e30; std::string best_city;
            for(auto& [name,obj]:env.all()){
                if(obj.body.mass!=0) continue;
                double d=vec_mag(vec_between(leo,pos(obj.body)));
                if(d<best){ best=d; best_city=name; }
            }
            double lst=local_solar_time_hours(0,jd);
            out<<"t "<<t<<" LeoClosest "<<best_city<<" dist_km "<<best<<" LST "<<lst<<"\n";
        }
    }
}
