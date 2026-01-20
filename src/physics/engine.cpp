#include "physics/engine.hpp"
#include <cmath>

static constexpr double MU=398600.4418; // km^3/s^2

void PhysicsEngine::add(const Body& b){ bodies.push_back(b); }

void PhysicsEngine::step(double dt){
    size_t n=bodies.size();
    if(n==0) return;

    std::vector<double> ax(n),ay(n),az(n);

    std::atomic<size_t> idx{0};
    auto worker_acc=[&](){
        for(;;){
            size_t i=idx.fetch_add(1);
            if(i>=n) break;
            auto& b=bodies[i];
            double r2=b.x*b.x+b.y*b.y+b.z*b.z;
            double r=std::sqrt(r2);
            double invr3=1.0/(r2*r);
            double k=-MU*invr3;
            ax[i]=k*b.x;
            ay[i]=k*b.y;
            az[i]=k*b.z;
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
            auto& b=bodies[i];
            b.vx+=ax[i]*dt;
            b.vy+=ay[i]*dt;
            b.vz+=az[i]*dt;
            b.x+=b.vx*dt;
            b.y+=b.vy*dt;
            b.z+=b.vz*dt;
        }
    };
    pool.clear();
    for(unsigned i=0;i<t;i++) pool.emplace_back(worker_step);
    for(auto& th:pool) th.join();
}
