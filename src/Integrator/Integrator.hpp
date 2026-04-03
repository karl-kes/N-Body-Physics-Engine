#pragma once

#include "../Particle/Particle.hpp"
#include "../Force/Force.hpp"

#include <vector>
#include <cstddef>
#include <memory>
#include <string>

#if defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define RESTRICT __restrict
#else
    #define RESTRICT
#endif

class Integrator {
protected:
    double dt_;
    std::string name_;
    
public:
    Integrator( double dt, std::string const &name ) 
    : dt_{ dt }
    , name_{ name }
    { }
    
    virtual ~Integrator() = default;
    virtual void integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const = 0;
    [[nodiscard]] double dt() const { return dt_; }
    [[nodiscard]] const std::string &name() const { return name_; }
};

class Velocity_Verlet : public Integrator {
public:
    Velocity_Verlet( double dt = 360.0 );
    void integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const override;
};

class Yoshida : public Integrator {
private:
    // Yoshida 4th-order coefficients.
    // w0 and w1 are pure constants: w1 = 1/(2 - 2^(1/3)), w0 = -2^(1/3)/(2 - 2^(1/3)).
    // c and d coefficients scale by dt, so they are computed at construction.
    static constexpr double cbrt_2_{ 1.2599210498948732 };
    static constexpr double w_0_{ -cbrt_2_ / ( 2.0 - cbrt_2_ ) };
    static constexpr double w_1_{ 1.0 / ( 2.0 - cbrt_2_ ) };

    double c_1_, c_2_, c_3_, c_4_;
    double d_1_, d_2_, d_3_;

public:
    Yoshida( double const dt = 900.0 );

    // Sequence per timestep: r(c1) -> v(d1) -> r(c2) -> v(d2) -> r(c3) -> v(d3) -> r(c4)
    // 4 position drifts interleaved with 3 force evaluations.
    void integrate( Particles &particles, std::vector<std::unique_ptr<Force>> const &forces ) const override;

    [[nodiscard]] static constexpr double cbrt_2() { return cbrt_2_; }
    [[nodiscard]] static constexpr double w_0() { return w_0_; }
    [[nodiscard]] static constexpr double w_1() { return w_1_; }

    [[nodiscard]] double c_1() const { return c_1_; }
    [[nodiscard]] double c_2() const { return c_2_; }
    [[nodiscard]] double c_3() const { return c_3_; }
    [[nodiscard]] double c_4() const { return c_4_; }

    [[nodiscard]] double d_1() const { return d_1_; }
    [[nodiscard]] double d_2() const { return d_2_; }
    [[nodiscard]] double d_3() const { return d_3_; }
};