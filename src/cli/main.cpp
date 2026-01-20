#include "physics/engine.hpp"
#include "sim/scenario.hpp"
#include "sim/sim.hpp"
#include "model/model.hpp"
#include "model/rocket_model.hpp"
#include "core/solarsystem.hpp"
#include "core/cities.hpp"
#include "core/output.hpp"
#include <string>
#include <cstring>
#include <cstdlib>

static const char* arg_after(int argc,char**argv,const char* key){
    for(int i=1;i+1<argc;i++) if(std::strcmp(argv[i],key)==0) return argv[i+1];
    return nullptr;
}

int main(int argc,char**argv){
    if(argc<3) return 1;

    const char* out_path = arg_after(argc,argv,"--output");
    const char* out_rate = arg_after(argc,argv,"--outputrate");
    if((out_path && !out_rate) || (!out_path && out_rate)) return 2;

    load_solar_system();
    load_global_cities();

    PhysicsEngine e;
    auto cfg=load_scenario(argv[2],e);

    OutputWriter ow;
    if(out_path && out_rate) ow.open(out_path,std::atof(out_rate));

    std::string mode=argv[1];

    if(mode=="--model"){
        std::string scen=argv[2];
        if(scen.find("rocket")!=std::string::npos){
            run_rocket_model(e,cfg.dt,cfg.t_end);
        }else{
            run_model(e,cfg.dt,cfg.t_end);
        }
        return 0;
    }

    if(mode=="--sim"){
        run_sim(e,cfg,1.0,"models/sim.simulation",(out_path&&out_rate)?&ow:nullptr);
        return 0;
    }

    return 0;
}
