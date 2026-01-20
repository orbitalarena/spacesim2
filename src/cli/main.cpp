#include "physics/engine.hpp"
#include "sim/scenario.hpp"
#include "model/model.hpp"
#include "core/solarsystem.hpp"
#include <string>

int main(int argc,char**argv){
    if(argc<3) return 1;
    load_solar_system();
    PhysicsEngine e;
    auto cfg=load_scenario(argv[2],e);
    if(std::string(argv[1])=="--model") run_model(e,cfg.dt,cfg.t_end);
    return 0;
}
