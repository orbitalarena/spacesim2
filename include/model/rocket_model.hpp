#pragma once
#include "physics/engine.hpp"
#include "core/output.hpp"

// New signature used by CLI:
void run_rocket_model(PhysicsEngine& e,double dt,double t_end,OutputWriter* ow);

// Back-compat wrapper (if anything else calls the old one):
inline void run_rocket_model(PhysicsEngine& e,double dt,double t_end){
    run_rocket_model(e,dt,t_end,nullptr);
}
