#include "model/tle_spawn.hpp"
#include "core/tle.hpp"
#include "core/tle_to_coe.hpp"
#include "physics/orbit.hpp"
#include "physics/engine.hpp"
#include <vector>

static constexpr double MU_E = 398600.4418; // km^3/s^2

void spawn_tle_entities(PhysicsEngine& e, const std::string& tle_path)
{
    std::vector<TLE> tles;
    if(!load_tles_from_file(tle_path, tles)) return;

    for(const auto& t : tles){
        COE c{};
        tle_to_coe(t, MU_E, c);

        // EXACT same pipeline as Ace
        Body b = coe_to_body_eci(c);
        e.bodies.push_back(b);
        e.names.push_back(t.name);
    }
}
