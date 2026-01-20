#include "physics/engine.hpp"
#include <cmath>
#include <atomic>
#include <thread>
#include <algorithm>

static constexpr double MU_E=398600.4418e9;

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
            double axi=0,ayi=0,azi=0;
            for(size_t j=0;j<n;++j){
                if(i==j) continue;
                const double dx=bodies[j].x-bodies[i].x;
                const double dy=bodies[j].y-bodies[i].y;
                const double dz=bodies[j].z-bodies[i].z;
                const double r2=dx*dx+dy*dy+dz*dz+1e-9;
                const double r=std::sqrt(r2);
                const double mu = (j==0) ? MU_E : 0.0;
                const double s = mu/(r2*r);
                axi += dx*s;
                ayi += dy*s;
                azi += dz*s;
            }
            ax[i]=axi; ay[i]=ayi; az[i]=azi;
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
