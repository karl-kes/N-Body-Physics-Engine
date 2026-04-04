#!/usr/bin/env bash
set -euo pipefail

# Unit Test & Benchmark Runner
# Builds and runs the test suite and optional scaling benchmark.
#
# Usage:
#   ./test.sh                # build and run unit tests
#   ./test.sh --benchmark    # also run the scaling benchmark
#   ./test.sh --bench-only   # skip tests, run benchmark only
#   ./test.sh --steps 200    # pass --steps to benchmark
#   ./test.sh --max-n 8192   # pass --max-n to benchmark

RUN_TESTS=true
RUN_BENCH=false
BENCH_ARGS=()

while [[ $# -gt 0 ]]; do
    case $1 in
        --benchmark)    RUN_BENCH=true;  shift ;;
        --bench-only)   RUN_TESTS=false; RUN_BENCH=true; shift ;;
        --steps)        BENCH_ARGS+=("--steps" "$2"); shift 2 ;;
        --max-n)        BENCH_ARGS+=("--max-n" "$2"); shift 2 ;;
        -h|--help)
            echo "Usage: ./test.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --benchmark    Run scaling benchmark after unit tests"
            echo "  --bench-only   Skip unit tests, run benchmark only"
            echo "  --steps N      Steps per trial for benchmark (default: 100)"
            echo "  --max-n N      Maximum N for benchmark sweep (default: 4096)"
            echo "  -h, --help     Show this help message"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Build
echo "Building..."
cmake -B build -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1

if [ "$RUN_TESTS" = true ]; then
    cmake --build build --target tests --config Release
fi

if [ "$RUN_BENCH" = true ]; then
    cmake --build build --target benchmark --config Release
fi

# Run tests
if [ "$RUN_TESTS" = true ]; then
    echo ""
    ./build/tests
    TEST_EXIT=$?

    if [ $TEST_EXIT -ne 0 ]; then
        echo "Tests failed."
        exit $TEST_EXIT
    fi
fi

# Run benchmark
if [ "$RUN_BENCH" = true ]; then
    echo ""
    ./build/benchmark "${BENCH_ARGS[@]}"
fi

echo ""
echo "Done."
