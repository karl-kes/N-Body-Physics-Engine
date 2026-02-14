#pragma once

#include "../Particle/Particle.hpp"
#include "../Config.hpp"

#include <vector>
#include <cmath>
#include <cstddef>
#include <omp.h>

#if defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define RESTRICT __restrict
#else
    #define RESTRICT
#endif

class Force {
public:
    virtual ~Force() = default;
    virtual void apply( Particles &particles ) const = 0;
};

class Gravity : public Force {
public:
    Gravity();
    static void compute_forces(
        double const pxi, double const pyi, double const pzi, double const mi,
        double const pxj, double const pyj, double const pzj, double const mj,
        double &a_xi, double &a_yi, double &a_zi,
        double const G, double const eps_sq, double const mask
    );

    void apply( Particles &particles ) const override;
};