#include "model/lambert_demo.hpp"
#include "orbit/lambert.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <algorithm>

static constexpr double MU_E_KM3_S2 = 398600.4418;

static inline std::array<double,3> r_of(const Body& b){ return {b.x,b.y,b.z}; }
static inline std::array<double,3> v_of(const Body& b){ return {b.vx,b.vy,b.vz}; }

static inline std::array<double,3> sub3(const std::array<double,3>& a,const std::array<double,3>& b){
    return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};
}
static inline double dot3(const std::array<double,3>& a,const std::array<double,3>& b){
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
static inline std::array<double,3> cross3(const std::array<double,3>& a,const std::array<double,3>& b){
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
static inline double norm3(const std::array<double,3>& a){ return std::sqrt(std::max(0.0, dot3(a,a))); }

static inline std::array<double,3> unit3_or(const std::array<double,3>& a, const std::array<double,3>& fallback){
    double n = norm3(a);
    if(n < 1e-12) return fallback;
    return {a[0]/n,a[1]/n,a[2]/n};
}

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

static inline void dv_to_ric_m_s(const std::array<double,3>& dv_km_s,
                                const std::array<double,3>& r_km,
                                const std::array<double,3>& v_km_s)
{
    auto eR = unit3_or(r_km, {1,0,0});

    double vr = dot3(v_km_s, eR);
    std::array<double,3> vt = { v_km_s[0]-vr*eR[0], v_km_s[1]-vr*eR[1], v_km_s[2]-vr*eR[2] };

    std::array<double,3> arb = (std::abs(eR[0]) < 0.9) ? std::array<double,3>{1,0,0} : std::array<double,3>{0,1,0};
    auto eI = unit3_or(vt, unit3_or(cross3(arb, eR), {0,1,0}));
    auto eC = unit3_or(cross3(eR, eI), {0,0,1});

    double dR = dot3(dv_km_s, eR) * 1000.0;
    double dI = dot3(dv_km_s, eI) * 1000.0;
    double dC = dot3(dv_km_s, eC) * 1000.0;

    std::cout << "DV_RIC_m_s " << dR << " " << dI << " " << dC << "\n";
}

static Body propagate_to_tof(const Body& b0, double tof_s){
    Body b = b0;
    const double dt = 10.0;
    double t = 0.0;
    while(t + dt < tof_s){
        step_body_central(b, dt);
        t += dt;
    }
    step_body_central(b, std::max(0.0, tof_s - t));
    return b;
}

static LambertSolution pick_min_dv1(const std::vector<LambertSolution>& sols,
                                   const std::array<double,3>& v_now_km_s)
{
    if(sols.empty()){
        return LambertSolution{{NAN,NAN,NAN},{NAN,NAN,NAN}};
    }
    double best = std::numeric_limits<double>::infinity();
    size_t best_i = 0;
    for(size_t i=0;i<sols.size();++i){
        auto dv = sub3(sols[i].v1_km_s, v_now_km_s);
        double m = norm3(dv);
        if(std::isfinite(m) && m < best){
            best = m;
            best_i = i;
        }
    }
    return sols[best_i];
}

static void run_case_intercept(Body ace0, Body chat0, PhysicsEngine& e){
    const double tof_s = 12.0*3600.0;     // 12 hours
    const int    max_rev = 3;

    Body chat_f = propagate_to_tof(chat0, tof_s);

    auto sols = lambert_intercept(r_of(ace0), r_of(chat_f), tof_s, MU_E_KM3_S2, max_rev);
    auto best = pick_min_dv1(sols, v_of(ace0));

    auto dv   = sub3(best.v1_km_s, v_of(ace0));
    std::cout << "=== INTERCEPT (TOF 43200s) ===\n";
    std::cout << "DV_mag_m_s " << (norm3(dv)*1000.0) << "\n";
    dv_to_ric_m_s(dv, r_of(ace0), v_of(ace0));

    // apply burn to ace
    ace0.vx = best.v1_km_s[0]; ace0.vy = best.v1_km_s[1]; ace0.vz = best.v1_km_s[2];

    // propagate 6 hours, printing every minute
    e.bodies.clear(); e.names.clear();
    e.add(chat0, "Chat");
    e.add(ace0,  "Ace");

    std::cout << "t_min range_km range_rate_km_s\n";
    double prev_range = -1.0;

    for(int m=0; m<=720; ++m){
        if(m>0) e.step(60.0);

        const Body& chat = e.bodies[0];
        const Body& ace  = e.bodies[1];

        const double dx = chat.x - ace.x;
        const double dy = chat.y - ace.y;
        const double dz = chat.z - ace.z;
        const double range = std::sqrt(dx*dx+dy*dy+dz*dz);

        double rr = 0.0;
        if(prev_range >= 0.0) rr = (range - prev_range) / 60.0;
        prev_range = range;

        std::cout << m << " " << range << " " << rr << "\n";
    }
}

static void run_case_rendezvous(Body ace0, Body chat0){
    const double tof_s = 12.0*3600.0;
    const int    max_rev = 3;

    Body chat_f = propagate_to_tof(chat0, tof_s);

    auto sols = lambert_rendezvous(r_of(ace0), v_of(chat_f), r_of(chat_f), tof_s, MU_E_KM3_S2, max_rev);
    auto best = pick_min_dv1(sols, v_of(ace0));

    auto dv   = sub3(best.v1_km_s, v_of(ace0));

    std::cout << "=== RENDEZVOUS (burn only, TOF 43200s) ===\n";
    std::cout << "DV_mag_m_s " << (norm3(dv)*1000.0) << "\n";
    dv_to_ric_m_s(dv, r_of(ace0), v_of(ace0));
}

static void run_case_nmc(Body ace0, Body chat0){
    const double tof_s = 12.0*3600.0;

    Body chat_f = propagate_to_tof(chat0, tof_s);

    NMCParams nmc{5.0, 45.0};
    auto sols = lambert_nmc(r_of(ace0), v_of(ace0), r_of(chat_f), v_of(chat_f), tof_s, nmc, MU_E_KM3_S2);
    auto best = pick_min_dv1(sols, v_of(ace0));

    auto dv   = sub3(best.v1_km_s, v_of(ace0));

    std::cout << "=== NMC (burn only, TOF 43200s) ===\n";
    std::cout << "DV_mag_m_s " << (norm3(dv)*1000.0) << "\n";
    dv_to_ric_m_s(dv, r_of(ace0), v_of(ace0));
}

void run_lambert_demo(PhysicsEngine& e)
{
    if(e.bodies.size() < 2){
        std::cerr << "lambert_demo: need 2 bodies (Chat, Ace)\n";
        return;
    }

    const Body chat0 = e.bodies[0];
    const Body ace0  = e.bodies[1];

    run_case_intercept(ace0, chat0, e);
    run_case_rendezvous(ace0, chat0);
    run_case_nmc(ace0, chat0);
}
