#include "core/solarsystem.hpp"
#include "core/environment.hpp"

static EnvEntity make(const char* n,double x,double y,double z,double m){
    EnvEntity e{};
    e.name=n;
    e.body.x=x; e.body.y=y; e.body.z=z;
    e.body.vx=0; e.body.vy=0; e.body.vz=0;
    e.body.mass=m;
    return e;
}

void load_solar_system(){
    auto& env=Environment::instance();

    env.add(make("Earth",0,0,0,5.972e24));
    env.add(make("Sun",-149597870,0,0,1.9885e30));
    env.add(make("Moon",384400,0,0,7.342e22));

    env.add(make("Mars",227939200,0,0,6.39e23));
    env.add(make("Jupiter",778340821,0,0,1.898e27));
    env.add(make("Io",778340821+421700,0,0,8.93e22));
    env.add(make("Europa",778340821+671034,0,0,4.8e22));
    env.add(make("Ganymede",778340821+1070412,0,0,1.48e23));
    env.add(make("Callisto",778340821+1882709,0,0,1.08e23));
    env.add(make("Saturn",1426666422,0,0,5.683e26));
    env.add(make("Uranus",2870658186,0,0,8.681e25));
    env.add(make("Neptune",4498396441,0,0,1.024e26));
    env.add(make("Pluto",5906380000,0,0,1.309e22));
}
