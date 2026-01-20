#pragma once
#include <fstream>
#include <string>
#include "physics/engine.hpp"

class OutputWriter {
public:
    OutputWriter()=default;
    void open(const std::string& path,double rate_s);
    void tick(double t,const PhysicsEngine& e);
    bool enabled() const { return is_on; }
private:
    std::ofstream out;
    double rate=0.0;
    double next_t=0.0;
    bool is_on=false;
};
