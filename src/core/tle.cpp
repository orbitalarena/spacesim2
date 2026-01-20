#include "core/tle.hpp"
#include <fstream>
#include <sstream>
#include <cmath>

static constexpr double DEG2RAD = M_PI / 180.0;
static constexpr double TWO_PI = 2.0 * M_PI;
static constexpr double SECONDS_PER_DAY = 86400.0;

static double parse_double(const std::string& s){
    return std::stod(s);
}

bool load_tles_from_file(const std::string& path, std::vector<TLE>& out){
    std::ifstream f(path);
    if(!f) return false;

    std::string line;
    while(true){
        std::string name,l1,l2;
        if(!std::getline(f,name)) break;
        if(!std::getline(f,l1)) break;
        if(!std::getline(f,l2)) break;

        if(l1.size() < 69 || l2.size() < 69) continue;

        TLE t;
        t.name = name;
        t.l1 = l1;
        t.l2 = l2;

        t.inc_rad  = parse_double(l2.substr(8,8))  * DEG2RAD;
        t.raan_rad = parse_double(l2.substr(17,8)) * DEG2RAD;
        t.ecc      = parse_double("0."+l2.substr(26,7));
        t.argp_rad = parse_double(l2.substr(34,8)) * DEG2RAD;
        t.M_rad    = parse_double(l2.substr(43,8)) * DEG2RAD;
        t.n_rev_per_day = parse_double(l2.substr(52,11));

        out.push_back(t);
    }
    return true;
}

static double solve_kepler(double M, double e){
    double E = M;
    for(int i=0;i<15;i++){
        double f = E - e*std::sin(E) - M;
        double fp = 1.0 - e*std::cos(E);
        E -= f/fp;
    }
    return E;
}

void tle_mean_to_eci(const TLE& t, double mu,
                     double& x, double& y, double& z,
                     double& vx, double& vy, double& vz)
{
    // mean motion [rad/s]
    double n = t.n_rev_per_day * TWO_PI / SECONDS_PER_DAY;

    // semi-major axis [km]
    double a = std::cbrt(mu / (n*n));

    double E = solve_kepler(t.M_rad, t.ecc);

    double cosE = std::cos(E);
    double sinE = std::sin(E);

    double r = a * (1.0 - t.ecc * cosE);

    double x_p = a * (cosE - t.ecc);
    double y_p = a * std::sqrt(1.0 - t.ecc*t.ecc) * sinE;

    double vx_p = -a * n * sinE / (1.0 - t.ecc*cosE);
    double vy_p =  a * n * std::sqrt(1.0 - t.ecc*t.ecc) * cosE
                   / (1.0 - t.ecc*cosE);

    double cosO = std::cos(t.raan_rad);
    double sinO = std::sin(t.raan_rad);
    double cosi = std::cos(t.inc_rad);
    double sini = std::sin(t.inc_rad);
    double cosw = std::cos(t.argp_rad);
    double sinw = std::sin(t.argp_rad);

    double R11 = cosO*cosw - sinO*sinw*cosi;
    double R12 = -cosO*sinw - sinO*cosw*cosi;
    double R21 = sinO*cosw + cosO*sinw*cosi;
    double R22 = -sinO*sinw + cosO*cosw*cosi;
    double R31 = sinw*sini;
    double R32 = cosw*sini;

    x  = R11*x_p + R12*y_p;
    y  = R21*x_p + R22*y_p;
    z  = R31*x_p + R32*y_p;

    vx = R11*vx_p + R12*vy_p;
    vy = R21*vx_p + R22*vy_p;
    vz = R31*vx_p + R32*vy_p;
}
