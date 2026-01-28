#pragma once

#include "Particle.hpp"
#include <vector>
#include <cmath>
#include <cstddef>

class Force_Law {
public:
    virtual ~Force_Law() = default;
    virtual void apply( Particles &particles ) const = 0;
};

class Gravity : public Force_Law {
private:
    double G_;

public:
    Gravity();

    void apply( Particles &particles ) const override;
    [[nodiscard]] double G() const { return G_; }
};