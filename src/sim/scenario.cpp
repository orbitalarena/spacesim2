#include "sim/scenario.hpp"
#include <fstream>

void load_scenario(const std::string& p, PhysicsEngine& e){
    std::ifstream f(p);
    Body b;
    while(f>>b.x>>b.y>>b.z>>b.vx>>b.vy>>b.vz>>b.mass) e.add(b);
}
