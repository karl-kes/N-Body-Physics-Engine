#include "Force/Force.hpp"
#include "Integrator/Integrator.hpp"
#include "Particle/Particle.hpp"
#include "Simulation/Simulation.hpp"
#include "Body.hpp"

#include <iostream>
#include <iomanip>
#include <cstddef>
#include <cmath>
#include <string>

int main() {
    /*

    Initialize:
    python src/jpl_compare.py fetch --moons

    Compile:
    rm -r -Force build
    cmake -B build -G "MinGW Makefiles"
    cmake --build build
    ./build/main.exe

    Validate:
    python src/jpl_compare.py compare

    Visualize:
    python src/visualize.py

    */

    // Determine number of bodies:
    static constexpr std::size_t num_bodies{ sizeof(bodies) / sizeof(bodies[0]) };

    // Prepare bodies:
    std::vector<std::string> names{};
    names.reserve( num_bodies );
    for ( std::size_t i{}; i < num_bodies; ++i ) {
        names.emplace_back( bodies[i].name );
    }

    // Initialize Simulation:
    Simulation sim{ 
        num_bodies,
        config::total_steps,
        config::output_interval,
        std::move( names ),
        "tests/sim_output.bin"
    };
    sim.add_force( std::make_unique<Gravity>() );
    sim.set_integrator( std::make_unique<Yoshida>( config::dt ) );
    initialize_bodies( sim.particles(), num_bodies );

    // Run simulation:
    sim.initial_output();
    sim.run();

    return 0;
}