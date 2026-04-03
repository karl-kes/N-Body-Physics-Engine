#!/usr/bin/env bash
set -euo pipefail

# N-Body Simulation Runner
# Configures, builds, and runs the full simulation pipeline.
#
# Usage:
#   ./run.sh                          # defaults: dt=900, years=249, output_hours=487, with moons
#   ./run.sh --dt 360 --years 100     # custom timestep and duration
#   ./run.sh --no-moons               # planets only (10 bodies)
#   ./run.sh --skip-fetch             # reuse existing JPL data
#   ./run.sh --visualize              # open Matplotlib viewer after sim
#   ./run.sh --render                 # open Rerun dashboard after sim

# Defaults:
DT="900.0"
YEARS="249"
OUTPUT_HOURS="487"
MOONS="--moons"
START_DATE="1950-01-01"
SKIP_FETCH=false
VISUALIZE=false
RENDER=false

# Parse arguments:
while [[ $# -gt 0 ]]; do
    case $1 in
        --dt)           DT="$2";           shift 2 ;;
        --years)        YEARS="$2";        shift 2 ;;
        --output-hours) OUTPUT_HOURS="$2"; shift 2 ;;
        --start)        START_DATE="$2";   shift 2 ;;
        --no-moons)     MOONS="";          shift ;;
        --skip-fetch)   SKIP_FETCH=true;   shift ;;
        --visualize)    VISUALIZE=true;    shift ;;
        --render)       RENDER=true;       shift ;;
        -h|--help)
            echo "Usage: ./run.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --dt VALUE           Integration timestep in seconds (default: 900.0)"
            echo "  --years VALUE        Simulation duration in years (default: 249)"
            echo "  --output-hours VALUE Output interval in hours (default: 487)"
            echo "  --start DATE         Start date YYYY-MM-DD (default: 1950-01-01)"
            echo "  --no-moons           Planets only (10 bodies instead of 35)"
            echo "  --skip-fetch         Skip JPL Horizons fetch (reuse existing data)"
            echo "  --visualize          Open Matplotlib viewer after simulation"
            echo "  --render             Open Rerun dashboard after simulation"
            echo "  -h, --help           Show this help message"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

CONFIG="src/Config.hpp"

# Validate:
REMAINDER=$(python3 -c "
dt = $DT
oh = $OUTPUT_HOURS
r = int(oh * 3600) % int(dt)
print(r)
")

if [[ "$REMAINDER" != "0" ]]; then
    echo "ERROR: output_hours * 3600 ($OUTPUT_HOURS * 3600) must be divisible by dt ($DT)"
    echo "       Remainder: $REMAINDER"
    exit 1
fi

# Patch Config.hpp:
echo "Configuring..."
echo "  dt:           $DT s"
echo "  years:        $YEARS"
echo "  output_hours: $OUTPUT_HOURS"
echo "  start:        $START_DATE"
echo "  moons:        $([ -n "$MOONS" ] && echo 'yes' || echo 'no')"

sed -i.bak -E "s/(inline static constexpr double dt\{)[^}]*/\1 ${DT} /" "$CONFIG"
sed -i.bak -E "s/(inline static constexpr std::size_t num_years\{)[^}]*/\1 ${YEARS} /" "$CONFIG"
sed -i.bak -E "s/(inline static constexpr std::size_t output_hours\{)[^}]*/\1 ${OUTPUT_HOURS} /" "$CONFIG"
rm -f "${CONFIG}.bak"

echo "  Config.hpp updated."

# Fetch JPL data:
if [ "$SKIP_FETCH" = false ]; then
    echo ""
    echo "Fetching JPL Horizons data..."
    python3 src/jpl_compare.py fetch --start "$START_DATE" $MOONS
else
    echo ""
    echo "Skipping JPL fetch (--skip-fetch)"
fi

# Build:
echo ""
echo "Building..."
cmake -B build -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1
cmake --build build --config Release

# Run:
echo ""
echo "Running simulation..."
./build/main

# Validate:
echo ""
echo "Validating against JPL Horizons..."
python3 src/jpl_compare.py compare

# Visualize:
if [ "$VISUALIZE" = true ]; then
    echo ""
    python3 src/visualize.py
fi

if [ "$RENDER" = true ]; then
    echo ""
    python3 src/render.py
fi

echo ""
echo "Done."
