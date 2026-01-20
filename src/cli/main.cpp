#include "physics/engine.hpp"
#include "sim/scenario.hpp"
#include "model/model.hpp"
#include "core/solarsystem.hpp"
#include "core/cities.hpp"
#include <string>

int main(int argc,char**argv){
    load_solar_system();
    load_global_cities();
    PhysicsEngine e;
    auto cfg=load_scenario(argv[2],e);
    run_model(e,cfg.dt,cfg.t_end);
    return 0;
}
