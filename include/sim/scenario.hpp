#pragma once
#include <string>
#include "physics/engine.hpp"

struct ScenarioCfg{
    double dt=1.0;
    double t_end=0.0;
};

ScenarioCfg load_scenario(const std::string& path,PhysicsEngine& e);
