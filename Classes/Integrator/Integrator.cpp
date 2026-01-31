#include "Integrator.hpp"

Velocity_Verlet::Velocity_Verlet( double dt )
: Integrator{ dt }
{ }

void Velocity_Verlet::integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const {
    std::size_t const N{ particles.num_particles() };

    for ( std::size_t i{}; i < N; ++i ) {
        auto &a{ particles };

        a.pos_x()[i] += dt() * ( a.vel_x()[i] + 0.5 * a.acc_x()[i] * dt() );
        a.pos_y()[i] += dt() * ( a.vel_y()[i] + 0.5 * a.acc_y()[i] * dt() );
        a.pos_z()[i] += dt() * ( a.vel_z()[i] + 0.5 * a.acc_z()[i] * dt() );

        a.old_acc_x()[i] = a.acc_x()[i];
        a.old_acc_y()[i] = a.acc_y()[i];
        a.old_acc_z()[i] = a.acc_z()[i];

        a.acc_x()[i] = 0.0;
        a.acc_y()[i] = 0.0;
        a.acc_z()[i] = 0.0;
    }

    for ( auto const &force : forces ) {
        force->apply( particles );
    }

    for ( std::size_t i{}; i < N; ++i ) {
        auto &a{ particles };

        a.vel_x()[i] += 0.5 * ( a.old_acc_x()[i] + a.acc_x()[i] ) * dt();
        a.vel_y()[i] += 0.5 * ( a.old_acc_y()[i] + a.acc_y()[i] ) * dt();
        a.vel_z()[i] += 0.5 * ( a.old_acc_z()[i] + a.acc_z()[i] ) * dt();
    }
}

Yoshida::Yoshida( double dt )
: Integrator{ dt }
, w_0_{ -cbrt_2() / ( 2.0 - cbrt_2() ) }
, w_1_{ 1.0 / ( 2.0 - cbrt_2() ) }
, c_1_{ w_1() / 2.0 }
, c_2_{ ( w_0() + w_1() ) / 2.0 }
, c_3_{ c_2() }
, c_4_{ c_1() }
, d_1_{ w_1() }
, d_2_{ w_0() }
, d_3_{ w_1() }
{ }

void Yoshida::integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const {
    std::size_t const N{ particles.num_particles() };

    auto calculate_pos = [&]( double c ){
        for ( std::size_t i{}; i < N; ++i ) {
            particles.pos_x()[i] += c * dt() * particles.vel_x()[i];
            particles.pos_y()[i] += c * dt() * particles.vel_y()[i];
            particles.pos_z()[i] += c * dt() * particles.vel_z()[i];
        }
    };

    auto apply_force = [&](){
        for ( std::size_t i{}; i < N; ++i ) {
            particles.acc_x()[i] = 0.0;
            particles.acc_y()[i] = 0.0;
            particles.acc_z()[i] = 0.0;
        }

        for ( auto const &force : forces ) {
            force->apply( particles );
        };
    };

    auto calculate_vel = [&]( double d ){
        for ( std::size_t i{}; i < N; ++i ) {
            particles.vel_x()[i] += d * dt() * particles.acc_x()[i];
            particles.vel_y()[i] += d * dt() * particles.acc_y()[i];
            particles.vel_z()[i] += d * dt() * particles.acc_z()[i];
        }
    };

    // Step 1:
    calculate_pos( c_1() );
    apply_force();
    calculate_vel( d_1() );

    // Step 2:
    calculate_pos( c_2() );
    apply_force();
    calculate_vel( d_2() );

    // Step 3:
    calculate_pos( c_3() );
    apply_force();
    calculate_vel( d_3() );

    // Step 4:
    calculate_pos( c_4() );
}