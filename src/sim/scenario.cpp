#include "sim/scenario.hpp"
#include "physics/orbit.hpp"
#include <fstream>
#include <string>
#include <cmath>

static constexpr double DEG=M_PI/180.0;

ScenarioCfg load_scenario(const std::string& p, PhysicsEngine& e){
    std::ifstream f(p);
    std::string tok;
    ScenarioCfg cfg{};
    COE c{};
    c.a=42166.3; c.e=0; c.i=0; c.raan=0; c.argp=0; c.ta=0;
    while(f>>tok){
        if(tok=="duration_seconds") f>>cfg.t_end;
        else if(tok=="timestep_seconds") f>>cfg.dt;
        else if(tok=="a_km") f>>c.a;
        else if(tok=="e") f>>c.e;
        else if(tok=="i_deg"){ f>>c.i; c.i*=DEG; }
        else if(tok=="raan_deg"){ f>>c.raan; c.raan*=DEG; }
        else if(tok=="argp_deg"){ f>>c.argp; c.argp*=DEG; }
        else if(tok=="ta_deg"){ f>>c.ta; c.ta*=DEG; }
    }
    e.add(coe_to_body_eci(c));
    return cfg;
}
