# N-Body Gravity Engine

A high-performance N-body gravitational simulator for solar system dynamics, written in C++23. Integrates 35 bodies (the Sun, 8 planets, Pluto, and 25 natural satellites) over multi-century timescales using a 4th-order Yoshida symplectic integrator. Validated against NASA JPL Horizons DE441 ephemeris data.

**[Technical Paper](docs/N_Body_Technical_Paper.pdf)** — Full derivations, validation methodology, and error analysis.

---

## Key Results

**249-year simulation (1950–2199), 35 bodies, Δt = 360 s, Yoshida 4th-order:**

| Metric | Value |
|---|---|
| Mean max relative error (ex Sun) | 0.079% |
| Worst-case body | Callisto (0.411%) |
| Planet error range | 0.0007% – 0.417% |
| Energy conservation ΔE/E | 1.5 × 10⁻¹² |
| Angular momentum ΔL/L | 4.3 × 10⁻¹³ |

The Sun's maximum relative error of 1.34% is a normalization artifact from the barycentric coordinate frame; its absolute position error is only 238 km, the smallest of any body in the simulation.

Residual errors are attributed to physics model differences (Newtonian gravity with 35 bodies vs. JPL DE441's ~300 bodies with post-Newtonian relativity, solar oblateness J₂, asteroid perturbations, and tidal effects), not numerical integration error. A timestep convergence study confirms that Δt ≤ 900 s lies on the floating-point precision floor.

---

## Quick Start

### Prerequisites

- C++23 compiler (GCC / Clang / MSVC)
- CMake 3.25+
- OpenMP
- Python 3.8+ with `numpy` and `matplotlib`
- Internet connection (for JPL Horizons API)

### Build & Run

```bash
# 1. Fetch JPL reference data and generate initial conditions
python src/jpl_compare.py fetch --moons

# 2. Build
cmake -B build
cmake --build build

# 3. Run simulation
./build/main

# 4. Validate against JPL Horizons
python src/jpl_compare.py compare

# 5. Visualize
python src/visualize.py
```

Step 1 queries the NASA JPL Horizons API for all 35 bodies, generates `src/Body.hpp` (initial conditions), and saves the reference ephemeris to `tests/`. This requires an internet connection and takes ~30 seconds.

Steps 2–5 work offline.

Alternatively, use the runner script which handles configuration, building, and validation in one command:

```bash
./scripts/run.sh                          # full pipeline with defaults
./scripts/run.sh --dt 360 --years 100     # custom timestep and duration
./scripts/run.sh --no-moons               # planets only (10 bodies)
./scripts/run.sh --visualize              # open viewer after sim
```

---

## Configuration

All simulation parameters live in `src/Config.hpp`:

```cpp
inline static constexpr double dt{ 900.0 };               // Integration timestep (seconds)
inline static constexpr std::size_t num_years{ 249 };     // Simulation duration
inline static constexpr std::size_t output_hours{ 487 };  // Output interval (hours)
```

Change these values, rebuild, and everything adapts automatically. The Python scripts read `num_years` and `output_hours` directly from this file. A `static_assert` enforces that `output_hours × 3600` is divisible by `dt`.

**Important:** After changing `num_years` or `output_hours`, re-run `python src/jpl_compare.py fetch --moons` to regenerate the JPL reference data at the new cadence before validating.

---

## Validation

`jpl_compare.py compare` reports per-body metrics in two tables:

- **Relative position error** — max and RMS relative error (%) for each body, with all-body and Sun-excluded means
- **Absolute position error** — max, RMS, and mean absolute error in km

Results are also exported to `tests/comparison_results.json`.

```bash
# Compare all bodies
python src/jpl_compare.py compare

# Compare specific bodies
python src/jpl_compare.py compare --bodies Mercury,Venus,Earth
```

---

## Visualization

### Matplotlib (interactive 3D viewer)

```bash
python src/visualize.py
python src/visualize.py --speed 4        # start at 4× speed
```

| Key | Action |
|---|---|
| Space | Play / Pause |
| Right / Left | Step forward / backward |
| Up / Down | Speed up (2×) / slow down (0.5×) |
| R | Reset to frame 0 |
| T | Toggle trails |
| Q | Quit |

### Rerun (dashboard with diagnostics)

```bash
pip install rerun-sdk
python src/render.py                     # open Rerun viewer
python src/render.py --save sim.rrd      # save to file
```

Displays 3D orbits alongside energy conservation, momentum drift, and heliocentric distance plots.

---

## Project Structure

```
├── src/
│   ├── main.cpp                # Entry point
│   ├── Config.hpp              # Single-source configuration
│   ├── Body.hpp                # Initial conditions (auto-generated)
│   ├── Force/                  # Gravity: O(N²), branchless SIMD kernel
│   ├── Integrator/             # Yoshida 4th-order + Velocity Verlet
│   ├── Particle/               # SoA particle data (single contiguous allocation)
│   ├── Simulation/             # Time-stepping loop, conservation diagnostics
│   ├── Output/                 # Binary output format
│   ├── jpl_compare.py          # JPL fetch + validation pipeline
│   ├── visualize.py            # Interactive 3D orbit viewer (Matplotlib)
│   └── render.py               # Rerun dashboard with diagnostics
├── scripts/
│   ├── run.sh                  # Full pipeline runner (Linux/macOS)
│   └── run.ps1                 # Full pipeline runner (Windows)
├── tests/                      # Generated validation data (gitignored)
├── docs/
│   └── N_Body_Technical_Paper.pdf
├── CMakeLists.txt
└── README.md
```

---

## Implementation Notes

**Yoshida 4th-order symplectic integrator.** Four position drifts interleaved with three force evaluations per timestep, achieving O(Δt⁴) local truncation error. The negative intermediate coefficient introduces a backward sub-step essential for error cancellation. Energy errors remain bounded and oscillatory over arbitrarily long integrations.

**Structure-of-Arrays memory layout.** All particle data occupies a single contiguous allocation (~3.6 KB for 35 bodies), fitting entirely in L1 cache. `__restrict__`-qualified pointers enable SIMD auto-vectorization.

**Branchless force kernel.** Self-interaction is eliminated with a floating-point mask rather than a conditional branch, preserving SIMD vectorization. Newton's third law symmetry is intentionally not exploited; the doubled FLOP count is traded for regular memory access patterns and freedom from race conditions under OpenMP.

**OpenMP parallelization.** Thread parallelism activates at N ≥ 500. SIMD vectorization of the inner loop is always active. At N = 1,050 (synthetic scaling benchmark), OpenMP yields 3.96× speedup on 6 cores.

---

## License

MIT