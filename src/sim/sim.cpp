#include "sim/sim.hpp"
#include "core/environment.hpp"
#include "core/vector.hpp"
#include "core/geodesy.hpp"
#include <fstream>
#include <chrono>
#include <thread>
#include <cmath>

static std::array<double,3> pos(const Body& b){ return {b.x,b.y,b.z}; }
static std::array<double,3> norm(const std::array<double,3>& v){
    double m=vec_mag(v);
    if(m==0) return {0,0,0};
    return {v[0]/m,v[1]/m,v[2]/m};
}

static void report_hourly(std::ofstream& out, double t, double jd, PhysicsEngine& e){
    auto& env=Environment::instance();

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

    out<<"t "<<t
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
    out<<"t "<<t<<" LeoClosest "<<best_city<<" dist_km "<<best<<" LST "<<lst<<"\n";
}

void run_sim(PhysicsEngine& e, const ScenarioCfg& cfg, double speed, const std::string& out_path){
    using clock = std::chrono::steady_clock;
    std::ofstream out(out_path);
    out.setf(std::ios::fixed);
    out.precision(6);

    double dt=cfg.dt;
    double t_end=cfg.t_end;
    if(speed<=0) speed=1.0;

    double jd=2451545.0;
    double sim_t=0.0;

    auto t0 = clock::now();
    auto next_hb = t0;
    auto next_step = t0;

    while(sim_t<=t_end){
        e.step(dt);
        sim_t += dt;
        jd += dt/86400.0;

        if(std::fmod(sim_t,3600.0)==0.0){
            report_hourly(out, sim_t, jd, e);
            out.flush();
        }

        auto now = clock::now();
        if(now>=next_hb){
            double wall_s = std::chrono::duration<double>(now - t0).count();
            out<<"heartbeat wall_s "<<wall_s<<" sim_s "<<sim_t<<" speed "<<speed<<"\n";
            out.flush();
            next_hb = now + std::chrono::seconds(10);
        }

        next_step += std::chrono::duration_cast<clock::duration>(std::chrono::duration<double>(dt/speed));
        std::this_thread::sleep_until(next_step);
    }
}
