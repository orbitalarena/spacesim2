#pragma once
#include <string>
#include <unordered_map>
#include "physics/engine.hpp"

struct EnvEntity {
    std::string name;
    Body body;
};

class Environment {
public:
    static Environment& instance();
    void add(const EnvEntity&);
    const EnvEntity& get(const std::string&) const;
    const std::unordered_map<std::string,EnvEntity>& all() const;
private:
    std::unordered_map<std::string,EnvEntity> entities;
};
