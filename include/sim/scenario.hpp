#pragma once
#include <string>
#include <unordered_map>
#include "physics/engine.hpp"

struct ScenarioCfg {
    double dt=60.0;
    double t_end=86400.0;
    std::unordered_map<std::string,size_t> body_index;
};

ScenarioCfg load_scenario(const std::string&, PhysicsEngine&);
