#include "Particle.hpp"

Particle::Particle(
    double mass,
    double init_pos_x, double init_pos_y, double init_pos_z,
    double init_vel_x, double init_vel_y, double init_vel_z,
    double init_acc_x, double init_acc_y, double init_acc_z )
    : mass_{ mass } {
    pos_[0] = init_pos_x;
    pos_[1] = init_pos_y;
    pos_[2] = init_pos_z;

    vel_[0] = init_vel_x;
    vel_[1] = init_vel_y;
    vel_[2] = init_vel_z;

    acc_[0] = init_acc_x;
    acc_[1] = init_acc_y;
    acc_[2] = init_acc_z;
}