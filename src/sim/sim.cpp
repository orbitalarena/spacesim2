#include "sim/sim.hpp"
#include "core/environment.hpp"
#include "core/vector.hpp"
#include "core/geodesy.hpp"
#include <fstream>
#include <chrono>
#include <thread>
#include <cmath>
#include <string>

static std::array<double,3> pos(const Body& b){ return {b.x,b.y,b.z}; }
static std::array<double,3> norm(const std::array<double,3>& v){
    double m=vec_mag(v);
    if(m==0) return {0,0,0};
    return {v[0]/m,v[1]/m,v[2]/m};
}

static double read_speed_file(const std::string& path, double cur){
    std::ifstream f(path);
    double s=cur;
    if(f.good() && (f>>s)){
        if(s<0.01) s=0.01;
        if(s>1e6) s=1e6;
        return s;
    }
    return cur;
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

    const std::string speed_path="models/sim_speed.txt";
    {
        std::ifstream chk(speed_path);
        if(!chk.good()){
            std::ofstream init(speed_path);
            init<<"1\n";
        }
    }

    std::ofstream out(out_path);
    out.setf(std::ios::fixed);
    out.precision(6);

    double sim_dt=1.0;                 // 1 sim-second per tick at 1x
    double jd=2451545.0;
    double sim_t=0.0;

    auto start = clock::now();
    auto next_hb = start;
    auto next_tick = start;

    double last_hour_mark=0.0;

    for(;;){
        if(sim_t>=cfg.t_end) break;

        speed = read_speed_file(speed_path, speed);
        if(speed<=0) speed=1.0;

        e.step(sim_dt);
        sim_t += sim_dt;
        jd += sim_dt/86400.0;

        if(sim_t - last_hour_mark >= 3600.0){
            last_hour_mark = std::floor(sim_t/3600.0)*3600.0;
            report_hourly(out, last_hour_mark, jd, e);
            out.flush();
        }

        auto now = clock::now();
        if(now>=next_hb){
            double wall_s = std::chrono::duration<double>(now - start).count();
            out<<"heartbeat wall_s "<<wall_s<<" sim_s "<<sim_t<<" speed "<<speed<<"\n";
            out.flush();
            next_hb = now + std::chrono::seconds(10);
        }

        // sleep based on current speed
        double sleep_s = sim_dt / speed;
        next_tick = next_tick + std::chrono::duration_cast<clock::duration>(std::chrono::duration<double>(sleep_s));
        std::this_thread::sleep_until(next_tick);
    }
}
