#include "Force.hpp"

Gravity::Gravity()
{ }

void Gravity::apply( Particles &particles ) const {
    std::size_t const N{ particles.num_particles() };

    for ( std::size_t i{}; i < N; ++i ) {
        auto &a{ particles };
        for ( std::size_t j{i + 1}; j < N; ++j ) {
            auto &b{ particles };

            double const dx{ b.pos_x()[j] - a.pos_x()[i] };
            double const dy{ b.pos_y()[j] - a.pos_y()[i] };
            double const dz{ b.pos_z()[j] - a.pos_z()[i] };

            double const R_sq{ dx*dx + dy*dy + dz*dz + constant::EPS*constant::EPS };
            double const inv_R_cb{ 1.0 / ( std::sqrt( R_sq ) * R_sq ) };

            double const factor_a{ constant::G * b.mass()[j] * inv_R_cb };
            double const factor_b{ constant::G * a.mass()[i] * inv_R_cb };


            a.acc_x()[i] += factor_a * dx;
            a.acc_y()[i] += factor_a * dy;
            a.acc_z()[i] += factor_a * dz;

            b.acc_x()[j] -= factor_b * dx;
            b.acc_y()[j] -= factor_b * dy;
            b.acc_z()[j] -= factor_b * dz;
        }
    }
}