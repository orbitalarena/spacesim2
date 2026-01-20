#include "physics/engine.hpp"

void PhysicsEngine::add(const Body& b){ bodies.push_back(b); }

void PhysicsEngine::step(double dt){
    size_t n=bodies.size();
    std::atomic<size_t> idx{0};
    auto worker=[&](){
        for(;;){
            size_t i=idx.fetch_add(1);
            if(i>=n) break;
            bodies[i].x+=bodies[i].vx*dt;
            bodies[i].y+=bodies[i].vy*dt;
            bodies[i].z+=bodies[i].vz*dt;
        }
    };
    unsigned t=std::thread::hardware_concurrency();
    std::vector<std::thread> pool;
    for(unsigned i=0;i<t;i++) pool.emplace_back(worker);
    for(auto& th:pool) th.join();
}
