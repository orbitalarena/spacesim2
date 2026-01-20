#pragma once
#include "physics/engine.hpp"
#include "sim/scenario.hpp"
#include "core/output.hpp"
void run_sim(PhysicsEngine& e,const ScenarioCfg& cfg,double speed,const char* simfile,OutputWriter* ow);
