#include "physics/engine.hpp"
#include "sim/scenario.hpp"
#include "sim/sim.hpp"
#include "model/model.hpp"
#include "model/rocket_model.hpp"
#include "core/solarsystem.hpp"
#include "core/cities.hpp"
#include <string>

int main(int argc,char**argv){
    if(argc<3) return 1;

    load_solar_system();
    load_global_cities();

    PhysicsEngine e;
    auto cfg=load_scenario(argv[2],e);

    std::string mode=argv[1];
    if(mode=="--model"){
        if(std::string(argv[2]).find("rocket")!=std::string::npos){
            run_rocket_model(e,cfg.dt,cfg.t_end);
        }else{
            run_model(e,cfg.dt,cfg.t_end);
        }
        return 0;
    }
    if(mode=="--sim"){
        run_sim(e,cfg,1.0,"models/sim.simulation");
        return 0;
    }
    return 0;
}
