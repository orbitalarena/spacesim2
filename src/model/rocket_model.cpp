#include "model/rocket_model.hpp"
#include "physics/rocket.hpp"
#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>

static constexpr double MU_E_KM3_S2 = 398600.4418; // km^3/s^2
static constexpr double R_E_KM = 6371.0;

static inline void step_body_central(Body& b, double dt){
    const double r2 = b.x*b.x + b.y*b.y + b.z*b.z;
    const double r  = std::sqrt(std::max(1e-12, r2));
    const double inv_r3 = 1.0/(r*r*r);
    const double ax = -MU_E_KM3_S2 * b.x * inv_r3;
    const double ay = -MU_E_KM3_S2 * b.y * inv_r3;
    const double az = -MU_E_KM3_S2 * b.z * inv_r3;
    b.vx += ax*dt; b.vy += ay*dt; b.vz += az*dt;
    b.x  += b.vx*dt; b.y  += b.vy*dt; b.z  += b.vz*dt;
}

struct PlanResult{
    double cost = std::numeric_limits<double>::infinity();
    double miss_km = std::numeric_limits<double>::infinity();
    double tca_s = 0.0;
    double thrust_scale = 0.25;
    double lead_tau_s = 4000.0;
};

static double simulate_coarse_cost(const PhysicsEngine& e0,
                                  double dt,
                                  double t_end,
                                  double thrust_scale,
                                  double lead_tau_s,
                                  double* out_miss_km,
                                  double* out_tca_s)
{
    std::vector<Stage> stages{
        {7.6e6*thrust_scale,263,395000,25600},
        {9.34e5*thrust_scale,421,92670,4000},
        {9.34e5*thrust_scale,450,15000,3500}
    };

    Rocket r(stages);

    Body ace = e0.bodies[0];
    const Body r0 = e0.bodies[1];

    RocketState rs{};
    rs.x=r0.x; rs.y=r0.y; rs.z=r0.z;
    rs.vx=r0.vx; rs.vy=r0.vy; rs.vz=r0.vz;
    rs.mass = (r0.mass>0 ? r0.mass : (stages[0].dry_kg+stages[0].fuel_kg));
    r.set_state(rs);

    double best_range = std::numeric_limits<double>::infinity();
    double best_t = 0.0;

    for(double t=0;t<=t_end;t+=dt){
        if(t>0.0) step_body_central(ace, dt);

        const double tp[3] = { ace.x + ace.vx*lead_tau_s,
                               ace.y + ace.vy*lead_tau_s,
                               ace.z + ace.vz*lead_tau_s };
        const double tv[3] = { ace.vx, ace.vy, ace.vz };

        r.step(dt, MU_E_KM3_S2, tp, tv, lead_tau_s);

        const auto s = r.state();
        const double dx = ace.x - s.x;
        const double dy = ace.y - s.y;
        const double dz = ace.z - s.z;
        const double range = std::sqrt(dx*dx + dy*dy + dz*dz);

        if(range < best_range){
            best_range = range;
            best_t = t;
        }
    }

    if(out_miss_km) *out_miss_km = best_range;
    if(out_tca_s)   *out_tca_s   = best_t;

    // Cost: prioritize miss, but prefer ~4 hour closest approach
    const double target_tca = 4.0*3600.0;
    const double t_pen = std::abs(best_t - target_tca);
    const double cost = best_range + 0.15*(t_pen/60.0);
    return cost;
}

static inline void orbital_diags(const RocketState& s, double& v_km_s, double& eps, double& a_km, double& e, double& rp_km, double& ra_km){
    const double rx=s.x, ry=s.y, rz=s.z;
    const double vx=s.vx, vy=s.vy, vz=s.vz;

    const double r = std::sqrt(rx*rx+ry*ry+rz*rz);
    const double v2 = vx*vx+vy*vy+vz*vz;
    v_km_s = std::sqrt(v2);

    eps = 0.5*v2 - MU_E_KM3_S2/r; // km^2/s^2

    // h = r x v
    const double hx = ry*vz - rz*vy;
    const double hy = rz*vx - rx*vz;
    const double hz = rx*vy - ry*vx;
    const double h  = std::sqrt(hx*hx+hy*hy+hz*hz);

    // evec = (v x h)/mu - r/r
    const double vxh_x = vy*hz - vz*hy;
    const double vxh_y = vz*hx - vx*hz;
    const double vxh_z = vx*hy - vy*hx;

    const double ex = vxh_x/MU_E_KM3_S2 - rx/r;
    const double ey = vxh_y/MU_E_KM3_S2 - ry/r;
    const double ez = vxh_z/MU_E_KM3_S2 - rz/r;
    e = std::sqrt(ex*ex+ey*ey+ez*ez);

    if(std::abs(eps) < 1e-12){
        a_km = std::numeric_limits<double>::infinity();
        rp_km = 0.0;
        ra_km = std::numeric_limits<double>::infinity();
        return;
    }

    a_km = -MU_E_KM3_S2/(2.0*eps);

    // For bound orbits (eps<0), rp/ra meaningful; for eps>0, rp meaningful, ra infinite (escape)
    if(eps < 0.0){
        rp_km = a_km*(1.0 - e);
        ra_km = a_km*(1.0 + e);
    }else{
        // hyperbolic / escape
        rp_km = a_km*(1.0 - e); // note: a is negative on hyperbola with this convention
        ra_km = std::numeric_limits<double>::infinity();
    }
}

static void final_run(const PhysicsEngine& e0, double dt, double t_end, double thrust_scale, double lead_tau_s, OutputWriter* ow){
    std::vector<Stage> stages{
        {7.6e6*thrust_scale,263,395000,25600},
        {9.34e5*thrust_scale,421,92670,4000},
        {9.34e5*thrust_scale,450,15000,3500}
    };

    Rocket r(stages);

    Body ace = e0.bodies[0];
    const Body r0 = e0.bodies[1];

    RocketState rs{};
    rs.x=r0.x; rs.y=r0.y; rs.z=r0.z;
    rs.vx=r0.vx; rs.vy=r0.vy; rs.vz=r0.vz;
    rs.mass = (r0.mass>0 ? r0.mass : (stages[0].dry_kg+stages[0].fuel_kg));
    r.set_state(rs);

    // Target: suborbital apogee ~ GEO radius (Earth radius + 35786 km)
    const double r_target_km = R_E_KM + 35786.0;
    bool cutoff_announced=false;

    double prev_range = -1.0;
    double last_print_t = -1e9;

    double best_range = std::numeric_limits<double>::infinity();
    double best_t = 0.0;

    for(double t=0;t<=t_end;t+=dt){
        if(t>0.0) step_body_central(ace, dt);

        const double tp[3] = { ace.x + ace.vx*lead_tau_s,
                               ace.y + ace.vy*lead_tau_s,
                               ace.z + ace.vz*lead_tau_s };
        const double tv[3] = { ace.vx, ace.vy, ace.vz };

        r.step(dt, MU_E_KM3_S2, tp, tv, lead_tau_s);

        auto s = r.state();

        // Suborbital cutoff: once we can reach GEO-radius apogee (two-body estimate), stop thrust and coast
        double v_km_s=0, eps=0, a_km=0, ecc=0, rp_km=0, ra_km=0;
        orbital_diags(s, v_km_s, eps, a_km, ecc, rp_km, ra_km);

        const double rmag = std::sqrt(s.x*s.x + s.y*s.y + s.z*s.z);
        const double alt_km = rmag - R_E_KM;

        if(!cutoff_announced && r.has_thrust() && alt_km > 150.0 && std::isfinite(ra_km) && ra_km >= r_target_km){
            r.force_cutoff();
            cutoff_announced=true;
            std::cout<<"t "<<t<<" ROCKET_CUTOFF reason apogee_reached ra_km "<<ra_km<<" target_ra_km "<<r_target_km<<"\n";
        }

        // Metrics vs Ace
        const double dx = ace.x - s.x;
        const double dy = ace.y - s.y;
        const double dz = ace.z - s.z;
        const double range = std::sqrt(dx*dx + dy*dy + dz*dz);

        if(range < best_range){
            best_range = range;
            best_t = t;
        }

        if(ow && ow->enabled()){
            PhysicsEngine tmp;
            tmp.bodies.resize(2);
            tmp.names.resize(2);
            tmp.bodies[0] = ace;
            tmp.names[0]  = "Ace";
            tmp.bodies[1] = Body{ s.x,s.y,s.z, s.vx,s.vy,s.vz, s.mass };
            tmp.names[1]  = "Rocket";
            ow->tick(t, tmp);
        }

        if((t - last_print_t) >= 60.0 - 1e-9){
            last_print_t = t;

            const double rvx = s.vx - ace.vx;
            const double rvy = s.vy - ace.vy;
            const double rvz = s.vz - ace.vz;

            const double rr_vec = (dx*rvx + dy*rvy + dz*rvz) / std::max(1e-9, range);

            double rr_num = 0.0;
            if(prev_range >= 0.0) rr_num = (range - prev_range) / 60.0;
            prev_range = range;

            std::cout<<"t "<<t<<" rocket_alt_km "<<alt_km
                     <<" rocket_v_km_s "<<v_km_s
                     <<" eps_km2_s2 "<<eps
                     <<" a_km "<<a_km
                     <<" e "<<ecc
                     <<" ra_km "<<ra_km
                     <<" rp_km "<<rp_km<<"\n";

            std::cout<<"t "<<t<<" rocket_to_ace_range_km "<<range
                     <<" rr_vec_km_s "<<rr_vec
                     <<" rr_num_km_s "<<rr_num
                     <<" lead_tau_s "<<lead_tau_s
                     <<" thrust_scale "<<thrust_scale
                     <<" powered "<<(r.has_thrust()?1:0)<<"\n";
        }
    }

    std::cout<<"CLOSEST_APPROACH t "<<best_t<<" range_km "<<best_range<<"\n";
}

void run_rocket_model(PhysicsEngine& e, double dt, double t_end, OutputWriter* ow){
    const double t_search = std::min(t_end, 6.0*3600.0);
    const double dt_search = std::max(2.0, dt);

    const double thrust_grid[] = {0.10,0.12,0.14,0.15,0.16,0.18,0.20,0.22,0.25,0.28,0.30};
    const double lead_grid[]   = {2500,3000,3500,4000,4500,5000,5500};

    PlanResult best;

    for(double ts : thrust_grid){
        for(double lt : lead_grid){
            double miss=0.0, tca=0.0;
            double cost = simulate_coarse_cost(e, dt_search, t_search, ts, lt, &miss, &tca);
            if(cost < best.cost){
                best.cost = cost;
                best.miss_km = miss;
                best.tca_s = tca;
                best.thrust_scale = ts;
                best.lead_tau_s = lt;
            }
        }
    }

    std::cout<<"AUTO_PLAN best_thrust_scale "<<best.thrust_scale
             <<" best_lead_tau_s "<<best.lead_tau_s
             <<" best_miss_km "<<best.miss_km
             <<" best_tca_s "<<best.tca_s<<"\n";

    final_run(e, dt, t_end, best.thrust_scale, best.lead_tau_s, ow);
}
