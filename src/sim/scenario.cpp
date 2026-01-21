#include "sim/scenario.hpp"
#include "physics/orbit.hpp"
#include "physics/engine.hpp"
#include "core/tle.hpp"
#include "core/geodesy.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>

static constexpr double MU_E_KM3_S2 = 398600.4418; // km^3/s^2
static constexpr double DEG2RAD = 3.14159265358979323846/180.0;

// Earth rotation rate (rad/s). If we treat ECEF as ECI at t=0, adding this gives surface co-rotation velocity.
// Units: (rad/s) * km = km/s.
static constexpr double OMEGA_E_RAD_S = 7.2921159e-5;

ScenarioCfg load_scenario(const std::string& path, PhysicsEngine& e){
    ScenarioCfg cfg{};
    cfg.dt = 1.0;
    cfg.t_end = 3600.0;

    struct Ent{
        std::string name;
        std::string type;
        COE coe{};
        bool has_coe=false;

        // Rocket launch params (degrees in scenario file)
        double launch_lat_deg = 0.0;
        double launch_lon_deg = 0.0;
        double launch_az_deg  = 90.0;
        bool has_launch = false;
    };

    std::vector<Ent> ents;
    std::string tle_path;

    std::ifstream f(path);
    if(!f){
        std::cerr << "load_scenario: failed to open " << path << "\n";
        return cfg;
    }

    std::string line;
    Ent cur{};
    bool in_entity=false;

    auto flush_entity = [&](){
        if(in_entity && !cur.name.empty()){
            ents.push_back(cur);
        }
        cur = Ent{};
        in_entity=false;
    };

    while(std::getline(f,line)){
        if(line.empty()) continue;
        std::istringstream ss(line);
        std::string k;
        ss >> k;
        if(k.empty()) continue;

        if(k=="duration_seconds"){ ss >> cfg.t_end; }
        else if(k=="timestep_seconds"){ ss >> cfg.dt; }
        else if(k=="entity"){
            flush_entity();
            in_entity=true;
            ss >> cur.name;
        }else if(k=="type"){
            ss >> cur.type;
        }else if(k=="coe"){
            cur.has_coe=true;
        }else if(k=="a_km"){
            ss >> cur.coe.a;
        }else if(k=="e"){
            ss >> cur.coe.e;
        }else if(k=="i_deg"){
            ss >> cur.coe.i; cur.coe.i *= DEG2RAD;
        }else if(k=="raan_deg"){
            ss >> cur.coe.raan; cur.coe.raan *= DEG2RAD;
        }else if(k=="argp_deg"){
            ss >> cur.coe.argp; cur.coe.argp *= DEG2RAD;
        }else if(k=="ta_deg"){
            ss >> cur.coe.ta; cur.coe.ta *= DEG2RAD;
        }else if(k=="tle_file"){
            ss >> tle_path;
        }
        // Rocket launch keys
        else if(k=="launch_lat"){
            ss >> cur.launch_lat_deg;
            cur.has_launch = true;
        }else if(k=="launch_lon"){
            ss >> cur.launch_lon_deg;
            cur.has_launch = true;
        }else if(k=="launch_az_deg"){
            ss >> cur.launch_az_deg;
            cur.has_launch = true;
        }
        // NOTE: stages are currently handled inside rocket_model.cpp (hardcoded). We parse launch only for now.
    }
    flush_entity();

    // Add explicit scenario entities
    for(const auto& en : ents){
        if(en.type=="satellite" && en.has_coe){
            Body b = coe_to_body_eci(en.coe);
            e.add(b, en.name);
        }else if(en.type=="rocket"){
            // Spawn rocket at the given launch site on Earth's surface (alt_km=0)
            // Using simple spherical Earth ECEF. Treat ECEF==ECI at t=0.
            double alt_km = 0.0;
            if(!en.has_launch){
                std::cerr << "rocket entity '" << en.name << "' missing launch_lat/launch_lon; defaulting to 0,0\n";
            }
            auto r_ecef = lla_to_ecef(en.launch_lat_deg, en.launch_lon_deg, alt_km);

            Body b;
            b.x = r_ecef[0];
            b.y = r_ecef[1];
            b.z = r_ecef[2];

            // Initial velocity: Earth rotation (co-rotating ground)
            // omega along +Z, v = omega x r = (-w*y, w*x, 0)
            b.vx = -OMEGA_E_RAD_S * b.y;
            b.vy =  OMEGA_E_RAD_S * b.x;
            b.vz = 0.0;

            b.mass = 0.0; // rocket_model.cpp will set mass from stage totals if mass<=0
            e.add(b, en.name);
        }
    }

    // Load TLE satellites (each becomes an entity)
    if(!tle_path.empty()){
        std::vector<TLE> tles;
        if(load_tles_from_file(tle_path, tles)){
            for(const auto& t : tles){
                double x,y,z,vx,vy,vz;
                tle_mean_to_eci(t, MU_E_KM3_S2, x,y,z,vx,vy,vz);
                Body b;
                b.x=x; b.y=y; b.z=z;
                b.vx=vx; b.vy=vy; b.vz=vz;
                b.mass=0;
                e.add(b, t.name);
            }
        }
    }

    return cfg;
}
