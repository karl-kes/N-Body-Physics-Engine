#include "Force/Force.hpp"
#include "Force/BarnesHut.hpp"
#include "Integrator/Integrator.hpp"
#include "Particle/Particle.hpp"
#include "Simulation/Simulation.hpp"
#include "Body.hpp"

#include <iostream>
#include <iomanip>
#include <cstddef>
#include <cmath>
#include <string>
#include <string_view>

int main( int argc, char* argv[] ) {
    std::string_view force_kind{ "direct" };
    double theta{ 0.5 };

    for ( int i{ 1 }; i < argc; ++i ) {
        std::string_view const arg{ argv[i] };
        if ( arg == "--force" && i + 1 < argc ) { force_kind = argv[++i]; }
        else if ( arg == "--theta" && i + 1 < argc ) { theta = std::stod( argv[++i] ); }
        else if ( arg == "-h" || arg == "--help" ) {
            std::cout << "Usage: main [--force {direct|bh}] [--theta T]\n"
                      << "  --force direct  Direct O(N^2) summation (default)\n"
                      << "  --force bh      Barnes-Hut O(N log N) approximation\n"
                      << "  --theta T       Opening angle for BH (default 0.5)\n";
            return 0;
        }
    }

    if ( force_kind != "direct" && force_kind != "bh" ) {
        std::cerr << "error: --force must be 'direct' or 'bh'\n";
        return 1;
    }

    static constexpr std::size_t num_bodies{ sizeof( bodies ) / sizeof( bodies[0] ) };

    std::vector<std::string> names{};
    names.reserve( num_bodies );
    for ( std::size_t i{}; i < num_bodies; ++i ) {
        names.emplace_back( bodies[i].name );
    }

    Simulation sim{
        num_bodies,
        config::total_steps,
        config::output_interval,
        std::move( names ),
        "tests/sim_output.bin"
    };

    if ( force_kind == "bh" ) {
        sim.add_force( std::make_unique<Gravity_BarnesHut>( theta ) );
    } else {
        sim.add_force( std::make_unique<Gravity>() );
    }
    sim.set_integrator( std::make_unique<Yoshida>( config::dt ) );
    initialize_bodies( sim.particles(), num_bodies );

    sim.run();

    return 0;
}
