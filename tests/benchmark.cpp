// tests/benchmark/benchmark.cpp
//
// Integration-step scaling benchmark: serial vs OpenMP across a range of N.
// Times full Yoshida integration steps (3 force evaluations + drift/kick updates).
// Reports wall time, estimated GFLOP/s, and speedup for each configuration.
//
// Build:  cmake --build build --target benchmark
// Run:    ./build/benchmark
//         ./build/benchmark --max-n 16384 --trials 5

#include "../src/Particle/Particle.hpp"
#include "../src/Force/Force.hpp"
#include "../src/Integrator/Integrator.hpp"
#include "../src/Config.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <vector>
#include <memory>
#include <random>
#include <string>
#include <algorithm>

#include <omp.h>

// Populate N bodies with random state. Enforces a minimum pairwise separation
// of `min_sep` meters to avoid near-singular force evaluations that would
// produce non-representative dynamics during the timed integration.
static void populate_random( Particles &p, std::size_t const N, double const min_sep = 1e9,
                             unsigned const seed = 42 ) {
    std::mt19937_64 rng{ seed };
    std::uniform_real_distribution<double> pos_dist{ -1e12, 1e12 };
    std::uniform_real_distribution<double> vel_dist{ -1e4, 1e4 };
    std::uniform_real_distribution<double> mass_dist{ 1e20, 1e30 };

    for ( std::size_t i{}; i < N; ++i ) {
        // Rejection-sample positions until minimum separation is satisfied.
        // For the densities used here, rejection is rare (< 1% of attempts).
        bool accepted{ false };
        while ( !accepted ) {
            double const px{ pos_dist( rng ) };
            double const py{ pos_dist( rng ) };
            double const pz{ pos_dist( rng ) };

            accepted = true;
            for ( std::size_t j{}; j < i; ++j ) {
                double const dx{ px - p.pos_x()[j] };
                double const dy{ py - p.pos_y()[j] };
                double const dz{ pz - p.pos_z()[j] };
                if ( dx*dx + dy*dy + dz*dz < min_sep * min_sep ) {
                    accepted = false;
                    break;
                }
            }

            if ( accepted ) {
                p.pos_x()[i] = px;
                p.pos_y()[i] = py;
                p.pos_z()[i] = pz;
            }
        }

        p.vel_x()[i] = vel_dist( rng );
        p.vel_y()[i] = vel_dist( rng );
        p.vel_z()[i] = vel_dist( rng );
        p.acc_x()[i] = 0.0;
        p.acc_y()[i] = 0.0;
        p.acc_z()[i] = 0.0;
        p.mass()[i]  = mass_dist( rng );
    }
}

struct BenchResult {
    std::size_t N;
    std::size_t steps;
    double serial_ms;
    double omp_ms;
    double serial_gflops;
    double omp_gflops;
    double speedup;
};

// Run `steps` Yoshida integration steps and return wall time in ms.
// Thread count is set before calling this function.
static double run_timed( std::size_t const N, std::size_t const steps ) {
    Particles p{ N };
    populate_random( p, N );

    std::vector<std::unique_ptr<Force>> forces;
    forces.push_back( std::make_unique<Gravity>() );
    Yoshida integ{ 900.0 };

    // Warmup: 3 steps to stabilize thread pool and caches.
    for ( std::size_t s{}; s < 3; ++s ) {
        integ.integrate( p, forces );
    }

    auto const start{ std::chrono::high_resolution_clock::now() };

    for ( std::size_t s{}; s < steps; ++s ) {
        integ.integrate( p, forces );
    }

    auto const end{ std::chrono::high_resolution_clock::now() };

    return std::chrono::duration<double, std::milli>( end - start ).count();
}

// Run multiple trials, return the median wall time.
static double run_median( std::size_t const N, std::size_t const steps,
                          int const num_threads, std::size_t const trials ) {
    omp_set_num_threads( num_threads );

    std::vector<double> timings;
    timings.reserve( trials );

    for ( std::size_t t{}; t < trials; ++t ) {
        timings.push_back( run_timed( N, steps ) );
    }

    std::sort( timings.begin(), timings.end() );
    return timings[trials / 2];
}

// Estimated FLOP count per Yoshida step.
// 3 force evaluations x N x N pairwise interactions x ~27 FLOPs per pair
// (sub, mul, add for dx/dy/dz, R_sq, 1/sqrt, mul chain, mask, accumulate)
// plus drift/kick updates: ~84N FLOPs per step.
// These are modeled estimates, not measured instruction counts.
static double estimated_flops_per_step( std::size_t const N ) {
    return 3.0 * N * N * 27.0 + 84.0 * N;
}

// Calibrate step count so serial runtime is approximately `target_ms`.
static std::size_t calibrate_steps( std::size_t const N, double const target_ms ) {
    omp_set_num_threads( 1 );
    double const calibration_ms{ run_timed( N, 5 ) };
    double const ms_per_step{ calibration_ms / 5.0 };
    std::size_t const estimated{ static_cast<std::size_t>( target_ms / ms_per_step ) };
    return std::clamp( estimated, static_cast<std::size_t>( 3 ), static_cast<std::size_t>( 500 ) );
}

int main( int argc, char* argv[] ) {
    std::size_t max_n{ 8192 };
    std::size_t num_trials{ 3 };
    double target_ms{ 2000.0 };

    // Capture max threads before any omp_set_num_threads calls.
    int const max_threads{ omp_get_max_threads() };

    for ( int i{ 1 }; i < argc; ++i ) {
        std::string const arg{ argv[i] };
        if ( arg == "--max-n" && i + 1 < argc ) { max_n = std::stoull( argv[++i] ); }
        else if ( arg == "--trials" && i + 1 < argc ) { num_trials = std::stoull( argv[++i] ); }
        else if ( arg == "--target-ms" && i + 1 < argc ) { target_ms = std::stod( argv[++i] ); }
        else if ( arg == "-h" || arg == "--help" ) {
            std::cout << "Usage: benchmark [--max-n N] [--trials N] [--target-ms MS]\n";
            std::cout << "  --max-n N       Maximum N for sweep (default: 8192)\n";
            std::cout << "  --trials N      Trials per config, reports median (default: 3)\n";
            std::cout << "  --target-ms MS  Target serial runtime per trial in ms (default: 2000)\n";
            return 0;
        }
    }

    std::cout << "\n<--- N-Body Scaling Benchmark --->\n";
    std::cout << "  Integrator:     Yoshida 4th-order (dt = 900 s)\n";
    std::cout << "  OMP threads:    " << max_threads << "\n";
    std::cout << "  Trials/config:  " << num_trials << " (median)\n";
    std::cout << "  Target serial:  " << std::fixed << std::setprecision( 0 ) << target_ms << " ms per trial\n";
    std::cout << "  OMP threshold:  N >= " << config::OMP_THRESHOLD << "\n\n";

    // Start at OMP_THRESHOLD (rounded up to next power of 2) so the parallel
    // code path in Force.cpp is always active.
    std::vector<std::size_t> N_values;
    std::size_t n_start{ 1 };
    while ( n_start < config::OMP_THRESHOLD ) n_start *= 2;
    for ( std::size_t n{ n_start }; n <= max_n; n *= 2 ) {
        N_values.push_back( n );
    }

    std::vector<BenchResult> results;

    std::cout << std::left
              << std::setw( 8 )  << "N"
              << std::setw( 8 )  << "Steps"
              << std::setw( 14 ) << "Serial (ms)"
              << std::setw( 14 ) << "OMP (ms)"
              << std::setw( 10 ) << "Speedup"
              << std::setw( 16 ) << "~Serial GFLOP/s"
              << std::setw( 16 ) << "~OMP GFLOP/s"
              << "\n";
    std::cout << std::string( 86, '=' ) << "\n";

    for ( std::size_t const N : N_values ) {
        std::size_t const steps{ calibrate_steps( N, target_ms ) };

        double const serial_ms{ run_median( N, steps, 1, num_trials ) };
        double const omp_ms{ run_median( N, steps, max_threads, num_trials ) };

        double const total_flops{ estimated_flops_per_step( N ) * steps };
        double const serial_gflops{ total_flops / ( serial_ms * 1e6 ) };
        double const omp_gflops{ total_flops / ( omp_ms * 1e6 ) };
        double const speedup{ serial_ms / omp_ms };

        results.push_back( { N, steps, serial_ms, omp_ms,
                             serial_gflops, omp_gflops, speedup } );

        std::cout << std::left << std::fixed
                  << std::setw( 8 )  << N
                  << std::setw( 8 )  << steps
                  << std::setw( 14 ) << std::setprecision( 1 ) << serial_ms
                  << std::setw( 14 ) << std::setprecision( 1 ) << omp_ms
                  << std::setw( 10 ) << std::setprecision( 2 ) << speedup
                  << std::setw( 16 ) << std::setprecision( 2 ) << serial_gflops
                  << std::setw( 16 ) << std::setprecision( 2 ) << omp_gflops
                  << "\n" << std::flush;
    }

    std::cout << std::string( 86, '=' ) << "\n\n";

    // Restore thread count
    omp_set_num_threads( max_threads );

    // CSV output for plotting
    std::cout << "CSV (for plotting):\n";
    std::cout << "N,steps,serial_ms,omp_ms,speedup,serial_gflops,omp_gflops\n";
    for ( auto const &r : results ) {
        std::cout << r.N << "," << r.steps << ","
                  << std::fixed << std::setprecision( 2 )
                  << r.serial_ms << "," << r.omp_ms << ","
                  << std::setprecision( 3 ) << r.speedup << ","
                  << std::setprecision( 3 ) << r.serial_gflops << ","
                  << r.omp_gflops << "\n";
    }

    return 0;
}