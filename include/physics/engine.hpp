#pragma once
#include <vector>
#include <thread>
#include <atomic>

struct Body {
    double x,y,z;
    double vx,vy,vz;
    double mass;
};

class PhysicsEngine {
public:
    void step(double dt);
    void add(const Body& b);
    std::vector<Body> bodies;
};
