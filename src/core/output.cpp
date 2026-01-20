#include "core/output.hpp"

void OutputWriter::open(const std::string& path,double rate_s){
    out.open(path);
    rate=rate_s;
    next_t=0.0;
    is_on=true;
    out<<"Time,EntityID,EntityName,X,Y,Z,VX,VY,VZ\n";
}

void OutputWriter::tick(double t,const PhysicsEngine& e){
    if(!is_on) return;
    if(t+1e-9 < next_t) return;
    next_t = t + rate;

    for(size_t i=0;i<e.bodies.size();++i){
        const auto& b=e.bodies[i];
        const char* name = (i<e.names.size()) ? e.names[i].c_str() : "";
        out<<t<<","<<i<<","<<name<<","<<b.x<<","<<b.y<<","<<b.z<<","<<b.vx<<","<<b.vy<<","<<b.vz<<"\n";
    }
}
