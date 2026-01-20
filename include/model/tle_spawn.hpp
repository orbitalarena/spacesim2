#pragma once
#include <string>

class PhysicsEngine;

// Spawn TLE satellites into the engine (Ace already exists)
void spawn_tle_entities(PhysicsEngine& e, const std::string& tle_path);
