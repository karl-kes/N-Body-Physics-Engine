#pragma once

#include "Body.hpp"
#include "Classes/Particle/Particle.hpp"

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
        Particles const &p,
        Body const *b,
        std::size_t const n,
        std::size_t const step,
        double const t ) {

        for (std::size_t i = 0; i < n; ++i) {
            file_ << step << "," << std::fixed << std::setprecision(2) << t << ","
               << b[i].name << "," << std::scientific << std::setprecision(12)
               << p.pos_x(i) << "," << p.pos_y(i) << "," << p.pos_z(i) << ","
               << p.vel_x(i) << "," << p.vel_y(i) << "," << p.vel_z(i) << "\n";
        }
        file_.flush();
    }
};