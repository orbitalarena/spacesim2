#include "sim/scenario.hpp"
#include "physics/orbit.hpp"
#include <fstream>
#include <sstream>
#include <cmath>


struct Ent { std::string name; COE c{}; };

ScenarioCfg load_scenario(const std::string& path,PhysicsEngine& e){
    ScenarioCfg cfg{};
    std::ifstream in(path);
    std::string line;

    std::vector<Ent> ents;
    Ent cur{};
    bool in_coe=false;

    while(std::getline(in,line)){
        if(line.empty()) continue;
        std::istringstream ss(line);
        std::string k; ss>>k;
        if(k=="duration_seconds"){ ss>>cfg.t_end; }
        else if(k=="timestep_seconds"){ ss>>cfg.dt; }
        else if(k=="entity"){
            if(!cur.name.empty()) ents.push_back(cur);
            cur=Ent{}; ss>>cur.name;
            in_coe=false;
        }else if(k=="coe"){ in_coe=true; }
        else if(in_coe){
            if(k=="a_km") ss>>cur.c.a;
            else if(k=="e") ss>>cur.c.e;
            else if(k=="i_deg"){ ss>>cur.c.i; cur.c.i*=M_PI/180.0; }
            else if(k=="raan_deg"){ ss>>cur.c.raan; cur.c.raan*=M_PI/180.0; }
            else if(k=="argp_deg"){ ss>>cur.c.argp; cur.c.argp*=M_PI/180.0; }
            else if(k=="ta_deg"){ ss>>cur.c.ta; cur.c.ta*=M_PI/180.0; }
        }
    }
    if(!cur.name.empty()) ents.push_back(cur);

    for(auto& en:ents){
        e.add(coe_to_body_eci(en.c),en.name);
    }
    return cfg;
}
