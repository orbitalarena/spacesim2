#include "core/cities.hpp"
#include "core/environment.hpp"
#include "core/geodesy.hpp"
#include "core/city.hpp"

static void add(const char* n,double lat,double lon){
    auto p=lla_to_ecef(lat,lon,0);
    EnvEntity e{};
    e.name=n;
    e.body.x=p[0]; e.body.y=p[1]; e.body.z=p[2];
    e.body.vx=0; e.body.vy=0; e.body.vz=0;
    e.body.mass=0;
    Environment::instance().add(e);
}

void load_global_cities(){
    add("LosAngeles",34.05,-118.25);
    add("WashingtonDC",38.907,-77.037);
    add("NewYork",40.7128,-74.0060);
    add("Tokyo",35.6895,139.6917);
    add("London",51.5074,-0.1278);
    add("Paris",48.8566,2.3522);
    add("Beijing",39.9042,116.4074);
    add("Shanghai",31.2304,121.4737);
    add("Mumbai",19.0760,72.8777);
    add("SaoPaulo",-23.55,-46.633);
}
