#include "model/rocket_model.hpp"
#include "physics/rocket.hpp"
#include <iostream>
#include <cmath>
#include <array>
#include <algorithm>

static constexpr double MU_E_KM3_S2 = 398600.4418;
static constexpr double R_E_KM      = 6371.0;

static double dot3(const std::array<double,3>& a,const std::array<double,3>& b){
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
static std::array<double,3> add3(const std::array<double,3>& a,const std::array<double,3>& b){
    return {a[0]+b[0],a[1]+b[1],a[2]+b[2]};
}
static std::array<double,3> sub3(const std::array<double,3>& a,const std::array<double,3>& b){
    return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};
}
static std::array<double,3> mul3(const std::array<double,3>& a,double s){
    return {a[0]*s,a[1]*s,a[2]*s};
}
static double mag3(const std::array<double,3>& v){
    return std::sqrt(dot3(v,v));
}

static double stumpC(double z){
    if(std::fabs(z) < 1e-8) return 0.5;
    if(z > 0){
        double sz = std::sqrt(z);
        return (1.0 - std::cos(sz)) / z;
    }else{
        double sz = std::sqrt(-z);
        return (1.0 - std::cosh(sz)) / z;
    }
}
static double stumpS(double z){
    if(std::fabs(z) < 1e-8) return 1.0/6.0;
    if(z > 0){
        double sz = std::sqrt(z);
        return (sz - std::sin(sz)) / (sz*sz*sz);
    }else{
        double sz = std::sqrt(-z);
        return (std::sinh(sz) - sz) / (sz*sz*sz);
    }
}

// Propagate (r0,v0) under 2-body for dt using universal variable formulation.
// Returns r,v at dt.
static void propagate_2body_universal(
    const std::array<double,3>& r0,
    const std::array<double,3>& v0,
    double dt,
    std::array<double,3>& r,
    std::array<double,3>& v
){
    double r0mag = mag3(r0);
    double v0mag2 = dot3(v0,v0);
    double vr0 = dot3(r0,v0) / r0mag;

    double alpha = 2.0/r0mag - v0mag2/MU_E_KM3_S2; // 1/a

    // Initial guess for chi
    double chi = 0.0;
    if(std::fabs(alpha) > 1e-10){
        chi = std::sqrt(MU_E_KM3_S2) * std::fabs(alpha) * dt;
    }else{
        chi = std::sqrt(MU_E_KM3_S2) * dt / r0mag;
    }
    chi = std::clamp(chi, 1e-6, 1e6);

    // Newton iterations
    for(int it=0; it<30; ++it){
        double chi2 = chi*chi;
        double z = alpha*chi2;
        double C = stumpC(z);
        double S = stumpS(z);

        double F = (r0mag*vr0/std::sqrt(MU_E_KM3_S2))*chi2*C
                 + (1.0 - alpha*r0mag)*chi*chi2*S
                 + r0mag*chi
                 - std::sqrt(MU_E_KM3_S2)*dt;

        double dF = (r0mag*vr0/std::sqrt(MU_E_KM3_S2))*chi*(1.0 - z*S)
                  + (1.0 - alpha*r0mag)*chi2*C
                  + r0mag;

        double dchi = -F / dF;
        chi += dchi;
        if(std::fabs(dchi) < 1e-10) break;
        if(chi < 1e-8) chi = 1e-8;
    }

    double chi2 = chi*chi;
    double z = alpha*chi2;
    double C = stumpC(z);
    double S = stumpS(z);

    double f = 1.0 - (chi2/r0mag)*C;
    double g = dt - (chi*chi2/std::sqrt(MU_E_KM3_S2))*S;

    r = add3(mul3(r0,f), mul3(v0,g));
    double rmag = mag3(r);

    double fdot = (std::sqrt(MU_E_KM3_S2)/(rmag*r0mag)) * (z*S - 1.0) * chi;
    double gdot = 1.0 - (chi2/rmag)*C;

    v = add3(mul3(r0,fdot), mul3(v0,gdot));
}

void run_rocket_model(PhysicsEngine& e,double dt,double t_end,OutputWriter* ow){
    if(e.bodies.size() < 2) return;

    auto& aceB = e.bodies[0];
    auto& rb   = e.bodies[1];

    // snapshot Ace initial state for clean 2-body prediction
    const std::array<double,3> ace_r0{aceB.x, aceB.y, aceB.z};
    const std::array<double,3> ace_v0{aceB.vx, aceB.vy, aceB.vz};

    std::vector<Stage> stages{
        {7.6e6,263,395000,25600},
        {9.34e5,421,92670,4000},
        {9.34e5,450,15000,3500}
    };

    Rocket r(stages);
    RocketState s0{};
    s0.x = rb.x; s0.y = rb.y; s0.z = rb.z;
    s0.vx = rb.vx; s0.vy = rb.vy; s0.vz = rb.vz;
    s0.mass = stages[0].fuel_kg + stages[0].dry_kg;
    r.set_state(s0);

    double prev_range = -1.0;
    double best_range = 1e300;
    double best_t = 0.0;

    for(double t=0; t<=t_end; t+=dt){
        // Predict Ace at current time t under 2-body (no manual linear drift)
        std::array<double,3> ace_r, ace_v;
        propagate_2body_universal(ace_r0, ace_v0, t, ace_r, ace_v);

        // Current rocket state
        auto rs = r.state();
        std::array<double,3> pR{rs.x, rs.y, rs.z};
        std::array<double,3> vR{rs.vx, rs.vy, rs.vz};

        // Relative now
        auto relp = sub3(ace_r, pR);
        auto relv = sub3(ace_v, vR);
        double range_km = mag3(relp);

        // Lead time = time-to-closest-approach under linear relative motion
        double relv2 = dot3(relv, relv);
        double tau = 0.0;
        if(relv2 > 1e-12){
            tau = -dot3(relp, relv) / relv2;
        }
        tau = std::clamp(tau, 0.0, 6.0*3600.0); // up to 6h lead

        // Predict Ace at t+tau (target point)
        std::array<double,3> ace_rt, ace_vt;
        propagate_2body_universal(ace_r0, ace_v0, t+tau, ace_rt, ace_vt);

        double target_xyz[3] = {ace_rt[0], ace_rt[1], ace_rt[2]};
        r.step(dt, MU_E_KM3_S2, target_xyz);

        // Write back into engine for CSV output
        rs = r.state();
        rb.x=rs.x; rb.y=rs.y; rb.z=rs.z;
        rb.vx=rs.vx; rb.vy=rs.vy; rb.vz=rs.vz;
        rb.mass=rs.mass;

        aceB.x=ace_r[0]; aceB.y=ace_r[1]; aceB.z=ace_r[2];
        aceB.vx=ace_v[0]; aceB.vy=ace_v[1]; aceB.vz=ace_v[2];

        if(ow && ow->enabled()) ow->tick(t, e);

        // Reporting (per-minute)
        double rmag = std::sqrt(rs.x*rs.x + rs.y*rs.y + rs.z*rs.z);
        double alt_km = rmag - R_E_KM;

        // Recompute range/rr to current Ace for reporting
        pR = {rs.x, rs.y, rs.z};
        vR = {rs.vx, rs.vy, rs.vz};
        relp = sub3(ace_r, pR);
        relv = sub3(ace_v, vR);
        range_km = mag3(relp);

        double rr_vec = 0.0;
        if(range_km > 1e-9) rr_vec = dot3(relp, relv) / range_km;

        double rr_num = 0.0;
        if(prev_range >= 0.0) rr_num = (range_km - prev_range) / dt;
        prev_range = range_km;

        if(range_km < best_range){
            best_range = range_km;
            best_t = t;
        }

        if(std::fmod(t,60.0)==0.0){
            std::cout<<"t "<<t<<" rocket_alt_km "<<alt_km<<"\n";
            std::cout<<"t "<<t<<" rocket_to_ace_range_km "<<range_km
                     <<" rr_vec_km_s "<<rr_vec
                     <<" rr_num_km_s "<<rr_num
                     <<" lead_tau_s "<<tau
                     <<"\n";
        }

        if(r.stage_sep()){
            std::cout<<"t "<<t<<" rocket_stage_sep\n";
        }

        if(range_km < 1.0){
            std::cout<<"ARRIVAL t "<<t
                     <<" range_km "<<range_km
                     <<" rr_vec_km_s "<<rr_vec
                     <<" rr_num_km_s "<<rr_num
                     <<"\n";
            break;
        }

        if(alt_km <= 0.0 && t > 10.0){
            std::cout<<"IMPACT t "<<t<<" rocket_alt_km "<<alt_km<<"\n";
            break;
        }
    }

    std::cout<<"CLOSEST_APPROACH t "<<best_t<<" range_km "<<best_range<<"\n";
}
