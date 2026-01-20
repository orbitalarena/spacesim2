#pragma once
#include <string>

// Reads an input .scenario, removes any `tle_file <path>` directives,
// and appends generated COE-based entities for each 3-line TLE entry.
// Returns true on success.
bool expand_scenario_with_tles(const std::string& in_path,
                               const std::string& out_path);
