#pragma once
#include <vector>

struct Stage{
    double thrust_n=0.0;
    double isp_s=0.0;
    double fuel_kg=0.0;
    double dry_kg=0.0;
};

struct RocketState{
    double x=0,y=0,z=0;     // km
    double vx=0,vy=0,vz=0;  // km/s
    double mass=0;          // kg
};

class Rocket{
public:
    Rocket()=default;
    Rocket(const std::vector<Stage>& stages);

    void set_state(const RocketState& s);
    RocketState state() const;

    // km-space propagation; target_xyz_km may be nullptr to fly "straight up" profile only
    void step(double dt_s, double mu_km3_s2, const double* target_xyz_km=nullptr);

    bool stage_sep() const { return sep; }
    int  stage_index() const { return (int)cur; }
    bool has_thrust() const { return powered; }

private:
    double mdot() const;

    std::vector<Stage> st;
    size_t cur=0;

    double x=0,y=0,z=0;
    double vx=0,vy=0,vz=0;

    double fuel=0.0;
    double dry=0.0;
    double mass=1.0;

    bool sep=false;
    bool powered=false;
};
