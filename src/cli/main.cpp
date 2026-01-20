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
    auto cfg = load_scenario(argv[2],e);

    OutputWriter ow;
    if(out_path && out_rate) ow.open(out_path,std::atof(out_rate));

    std::string mode = argv[1];

    if(mode=="--model"){
        std::string path = argv[2];
        if(path.find("rocket") != std::string::npos){
            run_rocket_model(e,cfg.dt,cfg.t_end,(out_path?&ow:nullptr));
        }else{
            run_model(e,cfg.dt,cfg.t_end,(out_path?&ow:nullptr));
        }
        return 0;
    }

    if(mode=="--sim"){
        const char* sp = arg_after(argc,argv,"--speed");
        double speed = sp ? std::atof(sp) : 1.0;
        run_sim(e,cfg,speed,nullptr,(out_path?&ow:nullptr));
        return 0;
    }

    return 1;
}
