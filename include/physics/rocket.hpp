#pragma once
#include <vector>
struct Stage{
    double thrust;
    double isp;
    double fuel;
    double dry;
};
struct RocketState{
    double m,x,y,z,vx,vy,vz;
};
class Rocket{
public:
    Rocket(const std::vector<Stage>&);
    bool step(double dt,double mu);
    RocketState state() const;
    bool stage_sep() const;
private:
    std::vector<Stage> stages;
    size_t cur;
    RocketState s;
    bool sep;
};
