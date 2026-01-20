#pragma once
#include <vector>
#include <array>

struct Stage{
    double thrust_n;      // N
    double isp_s;         // s
    double fuel_kg;
    double dry_kg;
};

struct RocketState{
    double x,y,z;
    double vx,vy,vz;
    double mass;
};

class Rocket{
public:
    explicit Rocket(const std::vector<Stage>& stages);

    void set_state(const RocketState& s);
    RocketState state() const;

    // mu is in m^3/s^2, dir_eci is desired thrust direction (will be normalized)
    void step(double dt,double mu_m3_s2,const std::array<double,3>& dir_eci);

    bool stage_sep();       // true on the tick a stage is dropped
    bool alive() const;     // has mass and stages

private:
    std::vector<Stage> st;
    size_t cur=0;
    bool sep=false;

    double x=0,y=0,z=0;
    double vx=0,vy=0,vz=0;

    double fuel=0;
    double dry=0;
    double mass=0;

    double mdot() const; // kg/s
};
