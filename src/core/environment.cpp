#include "core/environment.hpp"

Environment& Environment::instance(){
    static Environment e;
    return e;
}

void Environment::add(const EnvEntity& e){
    entities[e.name]=e;
}

const EnvEntity& Environment::get(const std::string& n) const{
    return entities.at(n);
}

const std::unordered_map<std::string,EnvEntity>& Environment::all() const{
    return entities;
}
