# ===============================
# README
# ===============================
cat > README.md <<'EOF'
# SpaceSim2

SpaceSim2 is a clean-slate rebuild of SpaceSim, designed to be minimal, fast, patchable, and physics-first.
The project is C++-based, CMake/Ninja driven, and structured to support real-time simulation and ultra-fast
model execution.

---

## Phase 0 – Reset and Foundation

- Full restart as **spacesim2**
- Hard baseline in C++
- Git-first workflow, patch-friendly
- Deterministic folder structure
- No UI assumptions

Core structure:
- `core/` – global environment, reference frames, shared utilities
- `physics/` – integrators, force models, orbital mechanics
- `sim/` – scenario parsing and real-time execution
- `model/` – max-speed batch execution
- `cli/` – command-line entry point
- `scenarios/` – user-defined scenario inputs
- `models/` – generated simulation outputs

---

## Phase 1 – Physics Engine

- Multi-threaded physics engine
- Two-body Newtonian gravity (Earth-centered)
- State propagation via velocity + acceleration
- Designed for extension to N-body and non-gravitational forces

---

## Phase 2 – Orbital Mechanics

- Classical Orbital Elements (COE)
- COE → ECI conversion
- ECI → COE recovery
- Units standardized to km / s / radians
- Validated GEO stability over 24-hour runs

---

## Phase 3 – Scenario System

- `.scenario` files define:
  - Duration
  - Timestep
  - Initial orbital elements
- Scenario loader converts COEs into physics entities
- Scenario controls both sim-mode and model-mode execution

---

## Phase 4 – Global Environment (Always-Loaded)

A permanent solar-system environment exists independently of scenarios.

Always available entities:
- Earth (origin / reference frame)
- Sun
- Moon
- Mars
- Jupiter
- Io, Europa, Ganymede, Callisto
- Saturn
- Uranus
- Neptune
- Pluto

Each entity has:
- Global XYZ position
- Mass
- Velocity placeholder (future ephemeris support)

These entities:
- Do not require scenario definitions
- Can be queried by any simulation
- Can optionally participate in gravity calculations later

---

## Phase 5 – Vector Geometry

- Global XYZ reference frame
- Vector construction between any two entities
- Magnitude, normalization, and angle calculations
- Stable angle math with clamping (no NaNs)

---

## Phase 6 – Model Execution

- `--model` runs at maximum possible speed
- Hourly reporting during execution
- Outputs:
  - Distance from Ace → every global body
  - Angle relative to Ace local radial vector
  - Final orbital elements at end of run

Output written to:
- `models/*.model`

---

## Current State

- GEO satellite remains stable at ~42166 km
- Earth-centered frame is correct
- Sun-relative geometry behaves correctly over 24 hours
- COEs remain consistent
- Architecture is extensible to:
  - N-body gravity
  - SRP
  - Atmospheric drag
  - Maneuvers
  - Real-time UI
  - Visualization backends

---

## Build

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build

##Run Example
./build/src/cli/spacesim2 --model scenarios/ace_geo.scenario

SpaceSim2 is now at a stable physics baseline.
Next steps: force models, ephemeris propagation, maneuvering, and UI.
EOF


```bash
# ===============================
# COMMIT + PUSH
# ===============================
git add .
git commit -m "spacesim2: README documenting full architecture and progress"
git branch -M main
git push -u origin main
