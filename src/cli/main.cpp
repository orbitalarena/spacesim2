#include "physics/engine.hpp"
#include "sim/scenario.hpp"
#include "sim/expand_scenario.hpp"
#include "model/model.hpp"
#include "model/rocket_model.hpp"
#include "model/tle_report.hpp"
#include "core/output.hpp"
#include "model/lambert_demo.hpp"
#include <iostream>
#include <string>

static void usage(){
    std::cerr << "usage:\n";
    std::cerr << "  spacesim2 --model <scenario> [--output <file> --outputrate <sec>]\n";
    std::cerr << "  spacesim2 --expand <in.scenario> <out.scenario>\n";
}

int main(int argc, char** argv){
    if(argc < 2){
        usage();
        return 2;
    }

    std::string mode = argv[1];

    if(mode == "--expand"){
        if(argc < 4){
            usage();
            return 2;
        }
        std::string in_path = argv[2];
        std::string out_path = argv[3];
        if(!expand_scenario_with_tles(in_path, out_path)){
            std::cerr << "expand failed: " << in_path << " -> " << out_path << "\n";
            return 1;
        }
        std::cout << "expanded scenario written: " << out_path << "\n";
        return 0;
    }

    if(mode != "--model"){
        usage();
        return 2;
    }

    if(argc < 3){
        usage();
        return 2;
    }

    std::string scenario_path = argv[2];

    std::string out_file;
    double out_rate = 0.0;

    for(int i=3;i<argc;i++){
        std::string a = argv[i];
        if(a=="--output" && i+1<argc) out_file = argv[++i];
        else if(a=="--outputrate" && i+1<argc) out_rate = std::stod(argv[++i]);
    }

    PhysicsEngine e;
    ScenarioCfg cfg = load_scenario(scenario_path, e);


    // Lambert demo
    if(scenario_path.find("lambert_demo") != std::string::npos){
        run_lambert_demo(e);
        return 0;
    }
    OutputWriter ow;
    OutputWriter* ow_ptr = nullptr;
    if(!out_file.empty() && out_rate > 0.0){
        ow.open(out_file, out_rate);
        ow_ptr = &ow;
    }

    if(scenario_path.find("tle_hour") != std::string::npos){
        run_tle_hour_report(e, cfg, ow_ptr);
        return 0;
    }
    if(scenario_path.find("rocket") != std::string::npos){
        run_rocket_model(e, cfg.dt, cfg.t_end, ow_ptr);
        return 0;
    }
    run_model(e, cfg.dt, cfg.t_end, ow_ptr);
    return 0;
}
