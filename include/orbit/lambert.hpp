#pragma once
#include <array>
#include <vector>

struct LambertSolution{
    std::array<double,3> v1_km_s;
    std::array<double,3> v2_km_s;
};

struct NMCParams{
    double semi_minor_km;   // closest approach (km)
    double entry_angle_deg; // degrees
};

std::vector<LambertSolution>
lambert_intercept(const std::array<double,3>& r1,
                  const std::array<double,3>& r2,
                  double tof_s,
                  double mu_km3_s2,
                  int max_rev);

std::vector<LambertSolution>
lambert_rendezvous(const std::array<double,3>& r1,
                   const std::array<double,3>& v_target_km_s,
                   const std::array<double,3>& r2,
                   double tof_s,
                   double mu_km3_s2,
                   int max_rev);

std::vector<LambertSolution>
lambert_nmc(const std::array<double,3>& r1,
            const std::array<double,3>& v1,
            const std::array<double,3>& r2,
            const std::array<double,3>& v2,
            double tof_s,
            const NMCParams& nmc,
            double mu_km3_s2);
