#include "physics/engine.hpp"
#include <cmath>
#include <atomic>
#include <thread>
#include <algorithm>
#include <vector>

// Units: km, km/s, s. Earth mu in km^3/s^2.
static constexpr double MU_E = 398600.4418;

void PhysicsEngine::add(const Body& b,const std::string& name){
    bodies.push_back(b);
    names.push_back(name);
}

void PhysicsEngine::step(double dt){
    const size_t n=bodies.size();
    if(n==0) return;

    std::vector<double> ax(n,0.0), ay(n,0.0), az(n,0.0);

    std::atomic<size_t> idx{0};
    auto worker_acc=[&](){
        for(;;){
            size_t i=idx.fetch_add(1);
            if(i>=n) break;

            const double x=bodies[i].x;
            const double y=bodies[i].y;
            const double z=bodies[i].z;

            const double r2 = x*x + y*y + z*z + 1e-12; // avoid 0-div
            const double r  = std::sqrt(r2);
            const double s  = -MU_E/(r2*r);           // -mu / r^3

            ax[i] = x*s;
            ay[i] = y*s;
            az[i] = z*s;
        }
    };

    unsigned t=std::max(1u,std::thread::hardware_concurrency());
    std::vector<std::thread> pool;
    pool.reserve(t);
    for(unsigned i=0;i<t;i++) pool.emplace_back(worker_acc);
    for(auto& th:pool) th.join();

    std::atomic<size_t> jdx{0};
    auto worker_step=[&](){
        for(;;){
            size_t i=jdx.fetch_add(1);
            if(i>=n) break;
            bodies[i].vx += ax[i]*dt;
            bodies[i].vy += ay[i]*dt;
            bodies[i].vz += az[i]*dt;
            bodies[i].x  += bodies[i].vx*dt;
            bodies[i].y  += bodies[i].vy*dt;
            bodies[i].z  += bodies[i].vz*dt;
        }
    };

    pool.clear();
    for(unsigned i=0;i<t;i++) pool.emplace_back(worker_step);
    for(auto& th:pool) th.join();
}
