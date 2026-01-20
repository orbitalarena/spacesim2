#pragma once
#include <vector>

struct Stage {
    double thrust;
    double isp;
    double fuel_mass;
    double dry_mass;
};

struct RocketState {
    double x,y,z;
    double vx,vy,vz;
};

class Rocket {
public:
    explicit Rocket(const std::vector<Stage>& stages);
    void step(double dt,double mu);
    RocketState state() const;
    bool stage_sep() const;
private:
    std::vector<Stage> stages;
    size_t active_stage;
    double mass;
    RocketState s;
    bool sep_flag;
};
