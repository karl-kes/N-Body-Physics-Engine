#include "Simulation.hpp"

Simulation::Simulation( std::size_t const num_particles,
                        std::size_t const steps, 
                        std::size_t const output_interval,
                        std::vector<std::string> names,
                        std::string output_path )
: particles_{ num_particles }
, forces_{}
, integrator_{ nullptr }
, num_bodies_{ num_particles }
, num_steps_{ steps }
, output_interval_{ output_interval }
, body_names_{ std::move( names ) }
, output_path_{ std::move( output_path ) }
{ }

void Simulation::run() {
    double const initial_energy{ total_energy() };
    double min_energy{ initial_energy };
    double max_energy{ initial_energy };

    double const initial_ang_momentum{ total_ang_momentum() };
    double min_ang_momentum{ initial_ang_momentum };
    double max_ang_momentum{ initial_ang_momentum };

    double const initial_lin_momentum{ total_lin_momentum() };
    double min_lin_momentum{ initial_lin_momentum };
    double max_lin_momentum{ initial_lin_momentum };

    Binary_Output bin{ output_path_, body_names_, num_bodies() };
    bin.write( particles(), num_bodies(), 0, 0.0 );

    auto const start_time{ std::chrono::high_resolution_clock::now() };

    for ( std::size_t curr_step{1}; curr_step <= steps(); ++curr_step ) {
        integrator()->integrate( particles(), forces() );

        if ( curr_step % ( 10*output_interval() ) == 0 ) {
            // Energy Tracking:
            double const E{ total_energy() };
            max_energy = std::max( E, max_energy );
            min_energy = std::min( E, min_energy );

            // Ang. Momentum Tracking:
            double const L{ total_ang_momentum() };
            max_ang_momentum = std::max( max_ang_momentum, L );
            min_ang_momentum = std::min( min_ang_momentum, L );

            // Lin. Momentum Tracking:
            double const P{ total_lin_momentum() };
            max_lin_momentum = std::max( max_lin_momentum, P );
            min_lin_momentum = std::min( min_lin_momentum, P );

            // Progress Bar:
            print_progress( curr_step, steps() );
        }

        if ( curr_step % output_interval() == 0 ) {
            bin.write( particles(), num_bodies(), curr_step, curr_step * config::dt );
        }
    }
    std::cout << "\rProgress: 100%" << std::flush;

    auto const end_time{ std::chrono::high_resolution_clock::now() };
    auto const duration{ std::chrono::duration_cast<std::chrono::milliseconds>( end_time - start_time ) };

    double const energy_drift{ std::abs( 100.0 * ( max_energy - min_energy ) / initial_energy ) };
    double const ang_momentum_drift{ std::abs( 100.0 * ( max_ang_momentum - min_ang_momentum ) / initial_ang_momentum ) };
    double const lin_momentum_drift{ std::abs( 100.0 * ( max_lin_momentum - min_lin_momentum ) / initial_lin_momentum ) };

    std::cout << "\nMax Energy Drift: " << std::scientific << std::setprecision( 6 ) << energy_drift << "%" << std::endl;
    std::cout << "Max Angular Momentum Drift: " << std::scientific << std::setprecision( 6 ) << ang_momentum_drift << "%" << std::endl;
    std::cout << "Max Linear Momentum Drift: " << std::scientific << std::setprecision( 6 ) << lin_momentum_drift << "%" << std::endl;
    std::cout << "Duration of Simulation: " << duration.count() << " ms" << std::endl;
}

void Simulation::add_force( std::unique_ptr<Force> force ) {
    forces().emplace_back( std::move( force ) );
}

void Simulation::set_integrator( std::unique_ptr<Integrator> sim_integrator ) {
    integrator() = std::move( sim_integrator );
}

double Simulation::total_energy() const {
    std::size_t const N{ particles().num_particles() };

    double const* RESTRICT px{ particles().pos_x() };
    double const* RESTRICT py{ particles().pos_y() };
    double const* RESTRICT pz{ particles().pos_z() };
    double const* RESTRICT vx{ particles().vel_x() };
    double const* RESTRICT vy{ particles().vel_y() };
    double const* RESTRICT vz{ particles().vel_z() };
    double const* RESTRICT mass{ particles().mass() };

    constexpr double eps_sq{ config::EPS * config::EPS };
    constexpr double G{ config::G };
    constexpr double OMP_THRESHOLD{ config::OMP_THRESHOLD };

    double kinetic_energy{};
    auto kinetic_kernel = [=]( std::size_t i ) -> double {
        double const v_sq{ vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i] };
        return 0.5 * mass[i] * v_sq;
    };

    double potential_energy{};
    auto potential_kernel = [=]( std::size_t i ) -> double {
        double const pxi{ px[i] }, pyi{ py[i] }, pzi{ pz[i] };
        double const mi{ mass[i] };
        double row_pot{ 0.0 };

        #pragma omp simd reduction( +:row_pot )
        for ( std::size_t j = i + 1; j < N; ++j ) {
            double const dx{ px[j] - pxi };
            double const dy{ py[j] - pyi };
            double const dz{ pz[j] - pzi };

            double const R_sq{ dx*dx + dy*dy + dz*dz + eps_sq };
            double const inv_R{ 1.0 / std::sqrt( R_sq ) };

            row_pot -= mass[j] * inv_R;
        }
        
        return G * mi * row_pot;
    };

    if ( N >= OMP_THRESHOLD ) {
        #pragma omp parallel for simd reduction( +:kinetic_energy ) schedule( static )
        for ( std::size_t i = 0; i < N; ++i ) {
            kinetic_energy += kinetic_kernel(i); 
        }
    } else {
        #pragma omp simd reduction( +:kinetic_energy )
        for ( std::size_t i = 0; i < N; ++i ) {
            kinetic_energy += kinetic_kernel(i); 
        }
    }

    if ( N >= OMP_THRESHOLD ) {
        #pragma omp parallel for reduction( +:potential_energy ) schedule( guided )
        for ( std::size_t i = 0; i < N; ++i ) {
            potential_energy += potential_kernel(i);
        }
    } else {
        for ( std::size_t i = 0; i < N; ++i ) {
            potential_energy += potential_kernel(i);
        }
    }

    return kinetic_energy + potential_energy;
}

double Simulation::total_ang_momentum() const {
    std::size_t const N{ particles().num_particles() };
    double const* RESTRICT m{ particles().mass() };

    double const* RESTRICT px{ particles().pos_x() };
    double const* RESTRICT py{ particles().pos_y() };
    double const* RESTRICT pz{ particles().pos_z() };

    double const* RESTRICT vx{ particles().vel_x() };
    double const* RESTRICT vy{ particles().vel_y() };
    double const* RESTRICT vz{ particles().vel_z() };

    double Lx{}, Ly{}, Lz{};
    #pragma omp simd reduction( +:Lx, Ly, Lz )
    for ( std::size_t i = 0; i < N; ++i ) {
        Lx += m[i] * ( py[i]*vz[i] - pz[i]*vy[i] );
        Ly += m[i] * ( pz[i]*vx[i] - px[i]*vz[i] );
        Lz += m[i] * ( px[i]*vy[i] - py[i]*vx[i] );
    }

    return std::sqrt( Lx*Lx + Ly*Ly + Lz*Lz );
}

double Simulation::total_lin_momentum() const {
    std::size_t const N{ particles().num_particles() };

    double const* RESTRICT m{ particles().mass() };

    double const* RESTRICT vx{ particles().vel_x() };
    double const* RESTRICT vy{ particles().vel_y() };
    double const* RESTRICT vz{ particles().vel_z() };

    double Px{}, Py{}, Pz{};
    #pragma omp simd reduction( +:Px, Py, Pz )
    for ( std::size_t i = 0; i < N; ++i ) {
        Px += m[i] * vx[i];
        Py += m[i] * vy[i];
        Pz += m[i] * vz[i];
    }

    return std::sqrt( Px*Px + Py*Py + Pz*Pz );
}

void Simulation::initial_output() {
    std::cout << "\n<--- Solar System Simulation --->" << std::endl;
    std::cout << "Bodies: " << num_bodies() << std::endl;
    std::cout << "Integrator: " << integrator()->name() << std::endl;
    std::cout << "Dt: " << config::dt << " seconds" << std::endl;
    std::cout << "Duration: " << config::num_years << " years" << std::endl;
    std::cout << "Parallelization: " << ( ( config::OMP_THRESHOLD <= num_bodies() ) ? "Enabled" : "Disabled" ) << std::endl;
    std::cout << std::endl;
}