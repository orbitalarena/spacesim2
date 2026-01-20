#include "core/geodesy.hpp"
#include <cmath>

static constexpr double R_E=6378.137;
static constexpr double DEG=M_PI/180.0;

std::array<double,3> lla_to_ecef(double lat,double lon,double alt){
    lat*=DEG; lon*=DEG;
    double r=R_E+alt;
    return {
        r*cos(lat)*cos(lon),
        r*cos(lat)*sin(lon),
        r*sin(lat)
    };
}

double local_solar_time_hours(double lon_deg,double jd){
    double T=(jd-2451545.0);
    double gst=fmod(18.697374558+24.06570982441908*T,24.0);
    double lst=fmod(gst+lon_deg/15.0,24.0);
    if(lst<0) lst+=24;
    return lst;
}

double ecef_lon_deg(double x,double y,double){
    double lon=atan2(y,x)*180.0/M_PI;
    return lon;
}
