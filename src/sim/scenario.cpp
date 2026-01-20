#include "sim/scenario.hpp"
#include "physics/orbit.hpp"
#include <fstream>
#include <string>

void load_scenario(const std::string& p, PhysicsEngine& e){
    std::ifstream f(p);
    std::string tok;
    COE c{};
    while(f>>tok){
        if(tok=="a_km") f>>c.a;
        else if(tok=="e") f>>c.e;
        else if(tok=="i_deg") f>>c.i;
        else if(tok=="raan_deg") f>>c.raan;
        else if(tok=="argp_deg") f>>c.argp;
        else if(tok=="ta_deg") f>>c.ta;
    }
    e.add(coe_to_body(c));
}
