#include "physics/engine.hpp"
#include <cmath>
#include <algorithm>

static constexpr double MU_E_KM3_S2 = 398600.4418; // km^3/s^2

void PhysicsEngine::add(const Body& b,const std::string& name){
    bodies.push_back(b);
    names.push_back(name);
}

void PhysicsEngine::step(double dt){
    const size_t n=bodies.size();
    if(n==0) return;

    for(size_t i=0;i<n;i++){
        Body& b = bodies[i];

        const double r2 = b.x*b.x + b.y*b.y + b.z*b.z;
        const double r  = std::sqrt(std::max(1e-12, r2));
        const double inv_r3 = 1.0/(r*r*r);

        const double ax = -MU_E_KM3_S2 * b.x * inv_r3;
        const double ay = -MU_E_KM3_S2 * b.y * inv_r3;
        const double az = -MU_E_KM3_S2 * b.z * inv_r3;

        // semi-implicit Euler
        b.vx += ax*dt; b.vy += ay*dt; b.vz += az*dt;
        b.x  += b.vx*dt; b.y  += b.vy*dt; b.z  += b.vz*dt;
    }
}
