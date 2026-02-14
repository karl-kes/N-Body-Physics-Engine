#pragma once

#include "Body.hpp"
#include "Particle/Particle.hpp"

#include <fstream>
#include <iomanip>

class CSV_Output {
private:
    std::ofstream file_;

public:
    explicit CSV_Output( std::string const &path )
    : file_{ path } {
        file_ << "step,time_s,name,x_m,y_m,z_m,vx_ms,vy_ms,vz_ms\n";
    }

    void write(
        Particles const &particles,
        Body const *bodies,
        std::size_t const num,
        std::size_t const step,
        double const time ) {

        for (std::size_t i = 0; i < num; ++i) {
            file_ << step << "," << std::fixed << std::setprecision(2) << time << ","
               << bodies[i].name << "," << std::scientific << std::setprecision(12)
               << particles.pos_x(i) << "," << particles.pos_y(i) << "," << particles.pos_z(i) << ","
               << particles.vel_x(i) << "," << particles.vel_y(i) << "," << particles.vel_z(i) << "\n";
        }
        file_.flush();
    }
};