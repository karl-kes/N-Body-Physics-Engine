#pragma once

#include <memory>

class Particles {
private:
    std::size_t N_;
    
    // SoA field indices into the contiguous memory block.
    // All 13 arrays are packed end-to-end: [pos_x|pos_y|...|mass], each of length N.
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

    // Single contiguous allocation: 13 * N doubles.
    // Eliminates per-array allocation overhead and guarantees spatial locality.
    std::unique_ptr<double[]> mem_block_;

public:
    explicit Particles( std::size_t const num_particles )
    : N_{ num_particles }
    , mem_block_{ std::make_unique<double[]>( NUM_SUB_ARRAYS * num_particles ) }
    { }

    // Move-only: unique_ptr member prevents copying implicitly,
    // but explicit declarations make ownership semantics clear.
    Particles( Particles&& ) = default;
    Particles& operator=( Particles&& ) = default;
    Particles( Particles const& ) = delete;
    Particles& operator=( Particles const& ) = delete;

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