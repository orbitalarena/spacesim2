#include "model/tle_report.hpp"
#include "physics/engine.hpp"
#include "core/geodesy.hpp"
#include <iostream>
#include <cmath>

void run_tle_hour_report(PhysicsEngine& e,
                         const ScenarioCfg& cfg,
                         OutputWriter* ow)
{
    // advance simulation normally
    double t = 0.0;
    while(t < cfg.t_end){
        e.step(cfg.dt);
        t += cfg.dt;
        if(ow && ow->enabled())
            ow->tick(t, e);
    }

    // Washington DC (lat, lon, alt km)
    auto dc = lla_to_ecef(38.9072, -77.0369, 0.0);

    std::cout << "\n--- Distance to Washington DC at t="
              << cfg.t_end << " ---\n";

    int ace_idx = -1;
    for(size_t i=0;i<e.names.size();i++){
        if(e.names[i] == "Ace"){
            ace_idx = (int)i;
            break;
        }
    }

    double best_range = 1e99;
    std::string best_name;

    for(size_t i=0;i<e.bodies.size();i++){
        const auto& b = e.bodies[i];

        double dx = b.x - dc[0];
        double dy = b.y - dc[1];
        double dz = b.z - dc[2];
        double r  = std::sqrt(dx*dx + dy*dy + dz*dz);

        std::cout << e.names[i] << "  " << r << " km\n";

        if((int)i != ace_idx){
            double ax = b.x - e.bodies[ace_idx].x;
            double ay = b.y - e.bodies[ace_idx].y;
            double az = b.z - e.bodies[ace_idx].z;
            double ar = std::sqrt(ax*ax + ay*ay + az*az);
            if(ar < best_range){
                best_range = ar;
                best_name  = e.names[i];
            }
        }
    }

    if(ace_idx >= 0){
        std::cout << "\nClosest to Ace: "
                  << best_name << "  "
                  << best_range << " km\n";
    }
}
