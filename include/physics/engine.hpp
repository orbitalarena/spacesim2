#pragma once
#include <vector>
#include <string>

struct Body{
    double x=0,y=0,z=0;
    double vx=0,vy=0,vz=0;
    double mass=0;
};

class PhysicsEngine{
public:
    std::vector<Body> bodies;
    std::vector<std::string> names;

    void add(const Body& b,const std::string& name);
    void step(double dt);
};
