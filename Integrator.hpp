#pragma once

#include "Particle.hpp"
#include "Force.hpp"
#include <vector>
#include <cstddef>
#include <memory>

class Integrator {
public:
    virtual ~Integrator() = default;
    virtual void integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const = 0;
};

class Velocity_Verlet : public Integrator {
private:
    double dt_; // Seconds

public:
    Velocity_Verlet( double dt = 1.0 );

    void integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const override;
    [[nodiscard]] double dt() const { return dt_; }
};