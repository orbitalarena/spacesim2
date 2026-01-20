#include "model/model.hpp"

void run_model(PhysicsEngine& e,double dt,double t_end){
    for(double t=0;t<t_end;t+=dt) e.step(dt);
}
