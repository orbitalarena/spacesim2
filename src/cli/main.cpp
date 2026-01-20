#include "physics/engine.hpp"
#include "sim/scenario.hpp"
#include "model/model.hpp"
#include <string>

int main(int argc,char**argv){
    if(argc<3) return 1;
    std::string mode=argv[1];
    PhysicsEngine e;
    auto cfg=load_scenario(argv[2],e);
    if(mode=="--model") run_model(e,cfg.dt,cfg.t_end);
    if(mode=="--sim") for(double t=0;t<cfg.t_end;t+=cfg.dt) e.step(cfg.dt);
    return 0;
}
