#include "Integrator.hpp"

Velocity_Verlet::Velocity_Verlet( double dt )
: dt_{ dt }
{ }

void Velocity_Verlet::integrate( Particles &particles, std::vector<std::unique_ptr<Force_Law>> const &forces ) const {
    std::size_t const N{ sizeof( particles.mass() ) / particles.mass()[0] };

    std::vector<double> old_acc_x{}, old_acc_y{}, old_acc_z{};
    old_acc_x.resize( N );
    old_acc_y.resize( N );
    old_acc_z.resize( N );

    for ( std::size_t i{}; i < N; ++i ) {
        auto &a{ particles };

        a.pos_x()[i] += dt() * ( a.vel_x()[i] + 0.5 * a.acc_x()[i] * dt() );
        a.pos_y()[i] += dt() * ( a.vel_y()[i] + 0.5 * a.acc_y()[i] * dt() );
        a.pos_z()[i] += dt() * ( a.vel_z()[i] + 0.5 * a.acc_z()[i] * dt() );

        old_acc_x[i] = a.acc_x()[i];
        old_acc_y[i] = a.acc_y()[i];
        old_acc_z[i] = a.acc_z()[i];

        a.acc_x()[i] = 0.0;
        a.acc_y()[i] = 0.0;
        a.acc_z()[i] = 0.0;
    }

    for ( auto const &force : forces ) {
        force->apply( particles );
    }

    for ( std::size_t i{}; i < N; ++i ) {
        auto &a{ particles };

        a.vel_x()[i] += 0.5 * ( old_acc_x[i] + a.acc_x()[i] ) * dt();
        a.vel_y()[i] += 0.5 * ( old_acc_y[i] + a.acc_y()[i] ) * dt();
        a.vel_z()[i] += 0.5 * ( old_acc_z[i] + a.acc_z()[i] ) * dt();
    }
}