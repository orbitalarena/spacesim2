#include "sim/scenario.hpp"
#include "physics/orbit.hpp"
#include "physics/rocket.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cctype>

static constexpr double DEG2RAD = 3.14159265358979323846/180.0;
static constexpr double R_E_KM  = 6371.0;

struct Ent{
    std::string name;
    std::string type;
    COE c{};
    double launch_lat=0, launch_lon=0, launch_az=90;
    std::vector<Stage> stages;
};

static bool is_blank(const std::string& s){
    for(char ch: s) if(!std::isspace((unsigned char)ch)) return false;
    return true;
}

ScenarioCfg load_scenario(const std::string& path, PhysicsEngine& e){
    ScenarioCfg cfg{};
    cfg.dt=1.0;
    cfg.t_end=86400.0;

    std::ifstream in(path);
    std::string line;

    std::vector<Ent> ents;
    Ent cur{};
    bool in_ent=false;

    while(std::getline(in,line)){
        if(is_blank(line)) continue;
        std::stringstream ss(line);
        std::string k; ss>>k;
        if(k=="name"){ continue; }
        else if(k=="duration_seconds"){ ss>>cfg.t_end; }
        else if(k=="timestep_seconds"){ ss>>cfg.dt; }
        else if(k=="entity"){
            if(in_ent) ents.push_back(cur);
            cur = Ent{};
            ss>>cur.name;
            in_ent=true;
        }else if(k=="type"){
            ss>>cur.type;
        }else if(k=="a_km"){ ss>>cur.c.a; }
        else if(k=="e"){ ss>>cur.c.e; }
        else if(k=="i_deg"){ ss>>cur.c.i; cur.c.i*=DEG2RAD; }
        else if(k=="raan_deg"){ ss>>cur.c.raan; cur.c.raan*=DEG2RAD; }
        else if(k=="argp_deg"){ ss>>cur.c.argp; cur.c.argp*=DEG2RAD; }
        else if(k=="ta_deg"){ ss>>cur.c.ta; cur.c.ta*=DEG2RAD; }
        else if(k=="launch_lat"){ ss>>cur.launch_lat; }
        else if(k=="launch_lon"){ ss>>cur.launch_lon; }
        else if(k=="launch_az_deg"){ ss>>cur.launch_az; }
        else if(k=="stage"){
            std::string tmp;
            int n=0;
            Stage st{};
            ss>>n;
            ss>>tmp>>st.thrust_n;
            ss>>tmp>>st.isp_s;
            ss>>tmp>>st.fuel_kg;
            ss>>tmp>>st.dry_kg;
            cur.stages.push_back(st);
        }
    }
    if(in_ent) ents.push_back(cur);

    for(const auto& en : ents){
        if(en.type=="satellite"){
            Body b = coe_to_body_eci(en.c);
            e.add(b,en.name);
        }else if(en.type=="rocket"){
            double lat= en.launch_lat * DEG2RAD;
            double lon= en.launch_lon * DEG2RAD;

            double cl=std::cos(lat), sl=std::sin(lat);
            double co=std::cos(lon), so=std::sin(lon);
            Body b{};
            b.x = R_E_KM * cl * co;
            b.y = R_E_KM * cl * so;
            b.z = R_E_KM * sl;

            b.vx=0; b.vy=0; b.vz=0;

            double m0=1.0;
            if(!en.stages.empty()){
                m0 = en.stages[0].fuel_kg + en.stages[0].dry_kg;
            }
            b.mass = m0;
            e.add(b,en.name);
        }
    }

    return cfg;
}
