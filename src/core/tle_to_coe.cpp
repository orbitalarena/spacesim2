#include "core/tle_to_coe.hpp"
#include <cmath>

static constexpr double TWO_PI = 2.0 * M_PI;
static constexpr double SECONDS_PER_DAY = 86400.0;

void tle_to_coe(const TLE& t, double mu, COE& c)
{
    // Mean motion [rad/s]
    double n = t.n_rev_per_day * TWO_PI / SECONDS_PER_DAY;

    // Semi-major axis [km]
    c.a = std::cbrt(mu / (n * n));

    c.e    = t.ecc;
    c.i    = t.inc_rad;
    c.raan = t.raan_rad;
    c.argp = t.argp_rad;

    // Mean anomaly ≈ true anomaly at epoch
    c.ta   = t.M_rad;
}
