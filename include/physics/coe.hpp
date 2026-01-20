#pragma once

// This project already defines COE + coe_to_body_eci in the
// orbital mechanics layer used by the scenario loader.
// Expose them here so new code can include physics/coe.hpp
// without changing existing includes.

#include "physics/orbit.hpp"
