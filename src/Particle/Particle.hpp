#pragma once

#include "aligned_soa.hpp"

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
    AlignedSoA<double> mem_block_;

public:
    explicit Particles( std::size_t const num_particles )
    : N_{ num_particles }
    , mem_block_{ num_particles, NUM_SUB_ARRAYS }
    { }

    // Move-only: unique_ptr member prevents copying implicitly,
    // but explicit declarations make ownership semantics clear.
    Particles( Particles&& ) = default;
    Particles& operator=( Particles&& ) = default;
    Particles( Particles const& ) = delete;
    Particles& operator=( Particles const& ) = delete;

    [[nodiscard]] std::size_t num_particles() const { return N_; }

    // Mutable raw pointers:
    [[nodiscard]] double* pos_x() { return mem_block_[POS_X]; }
    [[nodiscard]] double* pos_y() { return mem_block_[POS_Y]; }
    [[nodiscard]] double* pos_z() { return mem_block_[POS_Z]; }
    [[nodiscard]] double* vel_x() { return mem_block_[VEL_X]; }
    [[nodiscard]] double* vel_y() { return mem_block_[VEL_Y]; }
    [[nodiscard]] double* vel_z() { return mem_block_[VEL_Z]; }
    [[nodiscard]] double* acc_x() { return mem_block_[ACC_X]; }
    [[nodiscard]] double* acc_y() { return mem_block_[ACC_Y]; }
    [[nodiscard]] double* acc_z() { return mem_block_[ACC_Z]; }
    [[nodiscard]] double* old_acc_x() { return mem_block_[OLD_ACC_X]; }
    [[nodiscard]] double* old_acc_y() { return mem_block_[OLD_ACC_Y]; }
    [[nodiscard]] double* old_acc_z() { return mem_block_[OLD_ACC_Z]; }
    [[nodiscard]] double* mass() { return mem_block_[MASS]; }

    // Const raw pointers:
    [[nodiscard]] double const* pos_x() const { return mem_block_[POS_X]; }
    [[nodiscard]] double const* pos_y() const { return mem_block_[POS_Y]; }
    [[nodiscard]] double const* pos_z() const { return mem_block_[POS_Z]; }
    [[nodiscard]] double const* vel_x() const { return mem_block_[VEL_X]; }
    [[nodiscard]] double const* vel_y() const { return mem_block_[VEL_Y]; }
    [[nodiscard]] double const* vel_z() const { return mem_block_[VEL_Z]; }
    [[nodiscard]] double const* acc_x() const { return mem_block_[ACC_X]; }
    [[nodiscard]] double const* acc_y() const { return mem_block_[ACC_Y]; }
    [[nodiscard]] double const* acc_z() const { return mem_block_[ACC_Z]; }
    [[nodiscard]] double const* old_acc_x() const { return mem_block_[OLD_ACC_X]; }
    [[nodiscard]] double const* old_acc_y() const { return mem_block_[OLD_ACC_Y]; }
    [[nodiscard]] double const* old_acc_z() const { return mem_block_[OLD_ACC_Z]; }
    [[nodiscard]] double const* mass() const { return mem_block_[MASS]; }
};