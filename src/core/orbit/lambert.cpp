#include "orbit/lambert.hpp"
#include <cmath>
#include <algorithm>

static inline double stumpff_C(double z){
    if(z >  1e-12) return (1.0 - std::cos(std::sqrt(z))) / z;
    if(z < -1e-12) return (1.0 - std::cosh(std::sqrt(-z))) / z;
    return 0.5;
}

static inline double stumpff_S(double z){
    if(z >  1e-12) return (std::sqrt(z) - std::sin(std::sqrt(z))) / std::pow(z, 1.5);
    if(z < -1e-12) return (std::sinh(std::sqrt(-z)) - std::sqrt(-z)) / std::pow(-z, 1.5);
    return 1.0/6.0;
}

static inline double dot3(const std::array<double,3>& a, const std::array<double,3>& b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline double norm3(const std::array<double,3>& a){
    return std::sqrt(std::max(0.0, dot3(a,a)));
}

static LambertSolution lambert_uv_single(const std::array<double,3>& r1,
                                         const std::array<double,3>& r2,
                                         double tof_s,
                                         double mu_km3_s2,
                                         int k_rev)
{
    const double r1m = norm3(r1);
    const double r2m = norm3(r2);

    double cos_dnu = dot3(r1,r2) / (r1m*r2m);
    cos_dnu = std::clamp(cos_dnu, -1.0, 1.0);

    const double TWO_PI = 2.0 * M_PI;
    const double dnu = std::acos(cos_dnu) + TWO_PI * (double)k_rev;

    const double sin_dnu = std::sin(dnu);
    const double one_minus_cos = 1.0 - std::cos(dnu);

    LambertSolution sol{};
    if(std::abs(one_minus_cos) < 1e-12){
        sol.v1_km_s = {0,0,0};
        sol.v2_km_s = {0,0,0};
        return sol;
    }

    const double A = sin_dnu * std::sqrt((r1m*r2m)/one_minus_cos);

    // Solve for z (universal variable) using a robust bracket + relaxed Newton
    double z = 0.0;
    double z_lo = -4.0*M_PI*M_PI;
    double z_hi =  4.0*M_PI*M_PI;

    auto time_of_flight = [&](double zz)->double{
        const double C = stumpff_C(zz);
        const double S = stumpff_S(zz);

        const double sqrtC = std::sqrt(std::max(1e-16, C));
        const double y = r1m + r2m + A*(zz*S - 1.0)/sqrtC;
        if(y < 0.0) return std::numeric_limits<double>::infinity();

        const double x = std::sqrt(y / C);
        return (x*x*x*S + A*std::sqrt(y)) / std::sqrt(mu_km3_s2);
    };

    // bracket target tof
    for(int i=0;i<60;i++){
        const double t_lo = time_of_flight(z_lo);
        const double t_hi = time_of_flight(z_hi);
        if(std::isfinite(t_lo) && std::isfinite(t_hi) && t_lo <= tof_s && tof_s <= t_hi) break;
        z_lo *= 1.4;
        z_hi *= 1.4;
    }

    // binary + relaxed Newton steps
    for(int i=0;i<80;i++){
        const double t = time_of_flight(z);
        if(!std::isfinite(t)){
            z = 0.5*(z + z_hi);
            continue;
        }
        const double err = t - tof_s;
        if(std::abs(err) < 1e-6) break;

        if(err > 0) z_hi = z; else z_lo = z;

        // relaxed step toward bracket mid
        const double mid = 0.5*(z_lo + z_hi);
        z = 0.7*z + 0.3*mid;
    }

    const double C = stumpff_C(z);
    const double S = stumpff_S(z);
    const double sqrtC = std::sqrt(std::max(1e-16, C));
    const double y = r1m + r2m + A*(z*S - 1.0)/sqrtC;

    const double f = 1.0 - y/r1m;
    const double g = A * std::sqrt(y/mu_km3_s2);
    const double gdot = 1.0 - y/r2m;

    for(int i=0;i<3;i++){
        sol.v1_km_s[i] = (r2[i] - f*r1[i]) / g;
        sol.v2_km_s[i] = (gdot*r2[i] - r1[i]) / g;
    }
    return sol;
}

std::vector<LambertSolution>
lambert_intercept(const std::array<double,3>& r1,
                  const std::array<double,3>& r2,
                  double tof_s,
                  double mu_km3_s2,
                  int max_rev)
{
    std::vector<LambertSolution> out;
    max_rev = std::max(0, max_rev);
    for(int k=0;k<=max_rev;k++){
        out.push_back(lambert_uv_single(r1,r2,tof_s,mu_km3_s2,k));
    }
    return out;
}

std::vector<LambertSolution>
lambert_rendezvous(const std::array<double,3>& r1,
                   const std::array<double,3>& v_target_km_s,
                   const std::array<double,3>& r2,
                   double tof_s,
                   double mu_km3_s2,
                   int max_rev)
{
    auto out = lambert_intercept(r1,r2,tof_s,mu_km3_s2,max_rev);
    for(auto& s : out){
        s.v2_km_s = v_target_km_s; // terminal velocity match burn
    }
    return out;
}

std::vector<LambertSolution>
lambert_nmc(const std::array<double,3>& r1,
            const std::array<double,3>& v1,
            const std::array<double,3>& r2_target,
            const std::array<double,3>& v2_target,
            double tof_s,
            const NMCParams& nmc,
            double mu_km3_s2)
{
    auto dot3 = [](const std::array<double,3>& a, const std::array<double,3>& b){
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    };
    auto cross3 = [](const std::array<double,3>& a, const std::array<double,3>& b){
        return std::array<double,3>{
            a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]
        };
    };
    auto norm3 = [&](const std::array<double,3>& a){
        const double n2 = dot3(a,a);
        return std::sqrt(std::max(0.0, n2));
    };
    auto scale3 = [](const std::array<double,3>& a, double k){
        return std::array<double,3>{a[0]*k, a[1]*k, a[2]*k};
    };
    auto add3 = [](const std::array<double,3>& a, const std::array<double,3>& b){
        return std::array<double,3>{a[0]+b[0], a[1]+b[1], a[2]+b[2]};
    };
    auto unit3_or = [&](const std::array<double,3>& a, const std::array<double,3>& fb){
        const double n = norm3(a);
        if(n < 1e-12) return fb;
        return scale3(a, 1.0/n);
    };

    // Target LVLH basis at arrival (R/I/C)
    const auto eR = unit3_or(r2_target, {1,0,0});
    const auto h  = cross3(r2_target, v2_target);
    const auto eC = unit3_or(h, {0,0,1});
    const auto eI = unit3_or(cross3(eC, eR), {0,1,0});

    // Angular rate vector of LVLH frame: omega = h / r^2 (rad/s)
    const double rmag = std::max(1e-12, norm3(r2_target));
    const double r2   = rmag*rmag;
    const auto omega  = scale3(h, 1.0/r2);
    const double n    = norm3(omega); // mean-ish rate for CW-like NMC sizing

    // NMC: "2x1 ellipse" in (R,I) plane, size parameter is semi-minor axis b (closest approach)
    // Choose a=2b (semi-major in-track), param theta (entry angle).
    const double b_km = std::max(1e-9, nmc.semi_minor_km);
    const double a_km = 2.0*b_km;
    const double th   = nmc.entry_angle_deg * (3.14159265358979323846/180.0);

    // Relative position in LVLH (km)
    // x_R = b cos(th), y_I = a sin(th), z_C = 0
    const double xR = b_km * std::cos(th);
    const double yI = a_km * std::sin(th);
    const double zC = 0.0;

    // Relative velocity in LVLH rotating frame (km/s) for the parametric ellipse
    // xdot = -b n sin(th), ydot =  a n cos(th)
    const double xRdot = -b_km * n * std::sin(th);
    const double yIdot =  a_km * n * std::cos(th);
    const double zCdot = 0.0;

    // Convert LVLH-relative position to ECI
    const auto r_rel_eci = add3(add3(scale3(eR, xR), scale3(eI, yI)), scale3(eC, zC));
    const auto r2_des    = add3(r2_target, r_rel_eci);

    // Convert LVLH-relative velocity to ECI:
    // v_rel_eci = (xdot*eR + ydot*eI + zdot*eC) + omega x r_rel_eci
    const auto v_rel_basis = add3(add3(scale3(eR, xRdot), scale3(eI, yIdot)), scale3(eC, zCdot));
    const auto omega_x_r   = cross3(omega, r_rel_eci);
    const auto v_rel_eci   = add3(v_rel_basis, omega_x_r);
    const auto v2_des      = add3(v2_target, v_rel_eci);

    // Solve as "rendezvous to state": arrive at r2_des with v2_des at tof.
    // NOTE: we ignore multi-rev here (NMC sizing is already tied to local motion); caller can sweep if needed.
    return lambert_rendezvous(r1, v2_des, r2_des, tof_s, mu_km3_s2, /*max_rev=*/0);
}

