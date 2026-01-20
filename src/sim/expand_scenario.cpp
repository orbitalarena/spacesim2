#include "sim/expand_scenario.hpp"

#include "core/tle.hpp"
#include "core/tle_to_coe.hpp"
#include "physics/orbit.hpp"

#include <fstream>
#include <sstream>
#include <unordered_set>
#include <vector>
#include <cctype>
#include <cmath>

static constexpr double MU_E_KM3_S2 = 398600.4418;
static constexpr double RAD2DEG = 180.0 / M_PI;

static std::string trim(const std::string& s){
    size_t a = 0;
    while(a < s.size() && std::isspace((unsigned char)s[a])) a++;
    size_t b = s.size();
    while(b > a && std::isspace((unsigned char)s[b-1])) b--;
    return s.substr(a, b-a);
}

static bool starts_with(const std::string& s, const std::string& p){
    return s.size() >= p.size() && s.compare(0, p.size(), p) == 0;
}

static std::string sanitize_name(const std::string& name){
    std::string out;
    out.reserve(name.size());
    for(char ch : name){
        unsigned char c = (unsigned char)ch;
        if(std::isalnum(c)) out.push_back((char)c);
        else if(std::isspace(c) || ch=='-' || ch=='.') out.push_back('_');
        // drop everything else
    }
    if(out.empty()) out = "TLE_SAT";
    // collapse multiple underscores
    std::string collapsed;
    collapsed.reserve(out.size());
    bool last_us = false;
    for(char ch : out){
        if(ch=='_'){
            if(!last_us) collapsed.push_back(ch);
            last_us = true;
        } else {
            collapsed.push_back(ch);
            last_us = false;
        }
    }
    // trim underscores
    while(!collapsed.empty() && collapsed.front()=='_') collapsed.erase(collapsed.begin());
    while(!collapsed.empty() && collapsed.back()=='_') collapsed.pop_back();
    if(collapsed.empty()) collapsed = "TLE_SAT";
    return collapsed;
}

static std::string make_unique(std::unordered_set<std::string>& used, const std::string& base){
    if(!used.count(base)){
        used.insert(base);
        return base;
    }
    for(int i=2;i<1000000;i++){
        std::string c = base + "_" + std::to_string(i);
        if(!used.count(c)){
            used.insert(c);
            return c;
        }
    }
    // extremely unlikely
    return base + "_X";
}

bool expand_scenario_with_tles(const std::string& in_path,
                               const std::string& out_path)
{
    std::ifstream in(in_path);
    if(!in) return false;

    std::vector<std::string> kept_lines;
    kept_lines.reserve(2048);

    std::vector<std::string> tle_paths;

    std::string line;
    while(std::getline(in, line)){
        std::string t = trim(line);
        if(t.empty()){
            kept_lines.push_back(line);
            continue;
        }
        // strip comments starting with '#'
        auto hash = t.find('#');
        if(hash != std::string::npos){
            t = trim(t.substr(0, hash));
        }
        if(starts_with(t, "tle_file")){
            std::istringstream iss(t);
            std::string kw, path;
            iss >> kw >> path;
            if(!path.empty()){
                tle_paths.push_back(path);
            }
            // do NOT keep this line in the output (we expand it)
            continue;
        }
        kept_lines.push_back(line);
    }

    // load all TLEs referenced
    std::vector<TLE> all;
    for(const auto& p : tle_paths){
        std::vector<TLE> tles;
        if(load_tles_from_file(p, tles)){
            all.insert(all.end(), tles.begin(), tles.end());
        }
    }

    std::ofstream out(out_path);
    if(!out) return false;

    // write original scenario content without tle_file lines
    for(const auto& l : kept_lines) out << l << "\n";

    if(!all.empty()){
        out << "\n# --- Expanded TLE entities (" << all.size() << ") ---\n";
    }

    std::unordered_set<std::string> used_names;

    // try to seed used_names from existing "entity <name>" lines to avoid collisions
    for(const auto& l : kept_lines){
        std::string t = trim(l);
        if(starts_with(t, "entity ")){
            std::istringstream iss(t);
            std::string kw, nm;
            iss >> kw >> nm;
            if(!nm.empty()) used_names.insert(nm);
        }
    }

    // append each TLE as an entity with COE fields in degrees
    for(const auto& t : all){
        COE c{};
        tle_to_coe(t, MU_E_KM3_S2, c);

        std::string base = sanitize_name(t.name);
        std::string nm = make_unique(used_names, base);

        out << "\nentity " << nm << "\n";
        out << "type satellite\n";
        out << "coe\n";
        out << "a_km "    << c.a << "\n";
        out << "e "       << c.e << "\n";
        out << "i_deg "   << (c.i * RAD2DEG) << "\n";
        out << "raan_deg "<< (c.raan * RAD2DEG) << "\n";
        out << "argp_deg "<< (c.argp * RAD2DEG) << "\n";
        out << "ta_deg "  << (c.ta * RAD2DEG) << "\n";
    }

    return true;
}
