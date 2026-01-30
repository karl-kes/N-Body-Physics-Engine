#include "Simulation.hpp"

Simulation::Simulation( std::size_t num_particles, std::size_t steps, std::size_t output_interval )
: particles_{ num_particles }
, forces_{}
, integrator_{ nullptr }
, num_steps_{ steps }
, output_interval_{ output_interval }
{ }
