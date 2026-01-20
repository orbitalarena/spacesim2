#pragma once
#include <string>
#include "physics/engine.hpp"
struct ScenarioCfg { double dt=60.0; double t_end=86400.0; };
ScenarioCfg load_scenario(const std::string&, PhysicsEngine&);
