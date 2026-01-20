#pragma once
#include "physics/engine.hpp"
#include "sim/scenario.hpp"
#include "core/output.hpp"

void run_tle_hour_report(PhysicsEngine& e,
                         const ScenarioCfg& cfg,
                         OutputWriter* ow=nullptr);
