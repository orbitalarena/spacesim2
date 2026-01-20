#include "core/cities.hpp"
#include "core/environment.hpp"
#include "core/geodesy.hpp"

static void add(const char* name,double lat,double lon){
    EnvEntity e{};
    e.name=name;
    auto p=lla_to_ecef(lat,lon,0.0);
    e.body.x=p[0]; e.body.y=p[1]; e.body.z=p[2];
    e.body.vx=0; e.body.vy=0; e.body.vz=0;
    e.body.mass=0;
    Environment::instance().add(e);
}

void load_global_cities(){
    add("WashingtonDC",38.9072,-77.0369);
    add("LosAngeles",34.0522,-118.2437);
    add("NewYork",40.7128,-74.0060);
    add("London",51.5074,-0.1278);
    add("Tokyo",35.6895,139.6917);

    add("MexicoCity",19.4326,-99.1332);
    add("Delhi",28.7041,77.1025);
    add("Dhaka",23.8103,90.4125);
    add("Cairo",30.0444,31.2357);
    add("Karachi",24.8607,67.0011);
    add("BuenosAires",-34.6037,-58.3816);
    add("Istanbul",41.0082,28.9784);
    add("Kolkata",22.5726,88.3639);
    add("Lagos",6.5244,3.3792);
    add("Manila",14.5995,120.9842);
    add("RioDeJaneiro",-22.9068,-43.1729);
    add("Tianjin",39.3434,117.3616);
    add("Kinshasa",-4.4419,15.2663);
    add("Guangzhou",23.1291,113.2644);
    add("Moscow",55.7558,37.6173);
    add("Shenzhen",22.5431,114.0579);
    add("Lahore",31.5204,74.3587);
    add("Bangalore",12.9716,77.5946);
    add("Paris",48.8566,2.3522);
    add("Bogota",4.7110,-74.0721);
    add("Jakarta",-6.2088,106.8456);
    add("Chennai",13.0827,80.2707);
    add("Lima",-12.0464,-77.0428);
    add("Bangkok",13.7563,100.5018);
    add("Seoul",37.5665,126.9780);
    add("Nagoya",35.1815,136.9066);
    add("Hyderabad",17.3850,78.4867);
    add("Tehran",35.6892,51.3890);
    add("Chicago",41.8781,-87.6298);
    add("Chengdu",30.5728,104.0668);
    add("Nanjing",32.0603,118.7969);
    add("Wuhan",30.5928,114.3055);
    add("HoChiMinhCity",10.8231,106.6297);
    add("Luanda",-8.8390,13.2894);
    add("Ahmedabad",23.0225,72.5714);
    add("KualaLumpur",3.1390,101.6869);
    add("XiAn",34.3416,108.9398);
    add("HongKong",22.3193,114.1694);
    add("Dongguan",23.0207,113.7518);
    add("Hangzhou",30.2741,120.1551);
    add("Foshan",23.0215,113.1214);
    add("Shenyang",41.8057,123.4315);
    add("Riyadh",24.7136,46.6753);
    add("Baghdad",33.3152,44.3661);
}
