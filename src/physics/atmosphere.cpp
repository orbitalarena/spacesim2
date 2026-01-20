#include "physics/atmosphere.hpp"
#include <cmath>
double air_density(double h){
    if(h<0) h=0;
    if(h>100000) return 0.0;
    return 1.225*exp(-h/8500.0);
}
