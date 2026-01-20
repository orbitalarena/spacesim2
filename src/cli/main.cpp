#include "physics/engine.hpp"
#include "sim/scenario.hpp"
#include "sim/sim.hpp"
#include "model/model.hpp"
#include "core/solarsystem.hpp"
#include "core/cities.hpp"
#include <string>

static double arg_speed(int argc,char**argv,double def){
    for(int i=1;i<argc-1;i++){
        if(std::string(argv[i])=="--speed") return std::stod(argv[i+1]);
    }
    return def;
}

int main(int argc,char**argv){
    if(argc<3) return 1;

    load_solar_system();
    load_global_cities();

    PhysicsEngine e;
    auto cfg=load_scenario(argv[2],e);

    std::string mode=argv[1];
    if(mode=="--model"){
        run_model(e,cfg.dt,cfg.t_end);
        return 0;
    }
    if(mode=="--sim"){
        double speed=arg_speed(argc,argv,1.0);
        run_sim(e,cfg,speed,"models/ace_plus_leo.simulation");
        return 0;
    }
    return 1;
}
