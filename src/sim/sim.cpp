#include "sim/sim.hpp"
#include <chrono>
#include <thread>

void run_sim(PhysicsEngine& e,const ScenarioCfg& cfg,double speed,const char*,OutputWriter* ow){
    using clock=std::chrono::steady_clock;
    auto last=clock::now();
    double t=0.0;

    while(t<=cfg.t_end){
        auto now=clock::now();
        std::chrono::duration<double> wall=now-last;
        last=now;

        double dt_sim = wall.count()*speed;
        int n = (dt_sim>0) ? (int)(dt_sim/cfg.dt) : 0;
        if(n<1){ std::this_thread::sleep_for(std::chrono::milliseconds(1)); continue; }

        for(int i=0;i<n && t<=cfg.t_end;i++){
            e.step(cfg.dt);
            t+=cfg.dt;
            if(ow && ow->enabled()) ow->tick(t,e);
        }
    }
}
