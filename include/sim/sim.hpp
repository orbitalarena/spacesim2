#pragma once
#include <string>
#include "physics/engine.hpp"
#include "sim/scenario.hpp"
void run_sim(PhysicsEngine& e, const ScenarioCfg& cfg, double speed, const std::string& out_path);
