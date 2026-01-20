#pragma once
#include <string>
#include <vector>

struct TLE {
    std::string name;
    std::string l1;
    std::string l2;

    double inc_rad = 0.0;
    double raan_rad = 0.0;
    double ecc = 0.0;
    double argp_rad = 0.0;
    double M_rad = 0.0;
    double n_rev_per_day = 0.0;
};

bool load_tles_from_file(const std::string& path, std::vector<TLE>& out);

void tle_mean_to_eci(const TLE& t, double mu_km3_s2,
                     double& x_km, double& y_km, double& z_km,
                     double& vx_km_s, double& vy_km_s, double& vz_km_s);
