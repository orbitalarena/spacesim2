#include "sim/scenario.hpp"
#include "physics/orbit.hpp"
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

static constexpr double DEG=M_PI/180.0;

struct NamedCOE { std::string name; COE c; };

ScenarioCfg load_scenario(const std::string& p, PhysicsEngine& e){
    std::ifstream f(p);
    std::string tok;
    ScenarioCfg cfg{};

    std::vector<NamedCOE> ents;
    NamedCOE* cur=nullptr;

    while(f>>tok){
        if(tok=="duration_seconds"){ f>>cfg.t_end; continue; }
        if(tok=="timestep_seconds"){ f>>cfg.dt; continue; }

        if(tok=="entity"){
            std::string name; f>>name;
            NamedCOE n{};
            n.name=name;
            n.c.a=42166.3; n.c.e=0; n.c.i=0; n.c.raan=0; n.c.argp=0; n.c.ta=0;
            ents.push_back(n);
            cur=&ents.back();
            continue;
        }

        if(!cur) continue;

        if(tok=="a_km"){ f>>cur->c.a; }
        else if(tok=="e"){ f>>cur->c.e; }
        else if(tok=="i_deg"){ f>>cur->c.i; cur->c.i*=DEG; }
        else if(tok=="raan_deg"){ f>>cur->c.raan; cur->c.raan*=DEG; }
        else if(tok=="argp_deg"){ f>>cur->c.argp; cur->c.argp*=DEG; }
        else if(tok=="ta_deg"){ f>>cur->c.ta; cur->c.ta*=DEG; }
    }

    for(size_t i=0;i<ents.size();++i){
        e.add(coe_to_body_eci(ents[i].c));
        cfg.body_index[ents[i].name]=i;
    }
    return cfg;
}
