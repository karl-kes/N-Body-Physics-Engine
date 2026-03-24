#pragma once

#include <memory>

class Particles {
private:
    // Number of particles
    std::size_t N_;

    // Readable memory indexing:
    enum ArrayIndex : std::size_t {
        POS_X,
        POS_Y,
        POS_Z,
        VEL_X,
        VEL_Y,
        VEL_Z,
        ACC_X,
        ACC_Y,
        ACC_Z,
        OLD_ACC_X,
        OLD_ACC_Y,
        OLD_ACC_Z,
        MASS,
        NUM_SUB_ARRAYS
    };

    // Memory monolith:
    std::unique_ptr<double[]> mem_block_;

public:
    // Constructor:
    explicit Particles( std::size_t const num_particles )
    : N_{ num_particles }
    , mem_block_{ std::make_unique<double[]>( NUM_SUB_ARRAYS * num_particles ) }
    { }

    // Getters & Setters:
    [[nodiscard]] std::size_t num_particles() const { return N_; }

    // Mutable raw pointers:
    [[nodiscard]] double* pos_x() { return mem_block_.get() + POS_X * N_; }
    [[nodiscard]] double* pos_y() { return mem_block_.get() + POS_Y * N_; }
    [[nodiscard]] double* pos_z() { return mem_block_.get() + POS_Z * N_; }
    [[nodiscard]] double* vel_x() { return mem_block_.get() + VEL_X * N_; }
    [[nodiscard]] double* vel_y() { return mem_block_.get() + VEL_Y * N_; }
    [[nodiscard]] double* vel_z() { return mem_block_.get() + VEL_Z * N_; }
    [[nodiscard]] double* acc_x() { return mem_block_.get() + ACC_X * N_; }
    [[nodiscard]] double* acc_y() { return mem_block_.get() + ACC_Y * N_; }
    [[nodiscard]] double* acc_z() { return mem_block_.get() + ACC_Z * N_; }
    [[nodiscard]] double* old_acc_x() { return mem_block_.get() + OLD_ACC_X * N_; }
    [[nodiscard]] double* old_acc_y() { return mem_block_.get() + OLD_ACC_Y * N_; }
    [[nodiscard]] double* old_acc_z() { return mem_block_.get() + OLD_ACC_Z * N_; }
    [[nodiscard]] double* mass() { return mem_block_.get() + MASS * N_; }

    // Const raw pointers:
    [[nodiscard]] const double* pos_x() const { return mem_block_.get() + POS_X * N_; }
    [[nodiscard]] const double* pos_y() const { return mem_block_.get() + POS_Y * N_; }
    [[nodiscard]] const double* pos_z() const { return mem_block_.get() + POS_Z * N_; }
    [[nodiscard]] const double* vel_x() const { return mem_block_.get() + VEL_X * N_; }
    [[nodiscard]] const double* vel_y() const { return mem_block_.get() + VEL_Y * N_; }
    [[nodiscard]] const double* vel_z() const { return mem_block_.get() + VEL_Z * N_; }
    [[nodiscard]] const double* acc_x() const { return mem_block_.get() + ACC_X * N_; }
    [[nodiscard]] const double* acc_y() const { return mem_block_.get() + ACC_Y * N_; }
    [[nodiscard]] const double* acc_z() const { return mem_block_.get() + ACC_Z * N_; }
    [[nodiscard]] const double* old_acc_x() const { return mem_block_.get() + OLD_ACC_X * N_; }
    [[nodiscard]] const double* old_acc_y() const { return mem_block_.get() + OLD_ACC_Y * N_; }
    [[nodiscard]] const double* old_acc_z() const { return mem_block_.get() + OLD_ACC_Z * N_; }
    [[nodiscard]] const double* mass() const { return mem_block_.get() + MASS * N_; }
};