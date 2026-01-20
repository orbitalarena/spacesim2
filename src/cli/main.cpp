#include "physics/engine.hpp"
#include "sim/scenario.hpp"
#include "model/model.hpp"
#include <string>

int main(int argc,char**argv){
    if(argc<5) return 1;
    std::string mode=argv[1];
    PhysicsEngine e;
    load_scenario(argv[2],e);
    double dt=std::stod(argv[3]);
    double t_end=std::stod(argv[4]);
    if(mode=="--model") run_model(e,dt,t_end);
    if(mode=="--sim") for(double t=0;t<t_end;t+=dt) e.step(dt);
    return 0;
}
