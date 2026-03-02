#pragma once

#include <memory>

class Particles {
private:
    static constexpr std::size_t num_fields_{ 13 }; // Number of pointers
    std::size_t N_;                                 // Number of particles
    std::unique_ptr<double[]> mem_block_;           // Size of memory block for each pointer

public:
    // Constructor:
    explicit Particles( std::size_t const num_particles )
    : N_{ num_particles } {
        mem_block_ = std::make_unique<double[]>( num_fields_*num_particles );
    }

    // Getters & Setters:
    [[nodiscard]] std::size_t num_particles() const { return N_; }

    // Mutable raw pointers:
    [[nodiscard]] double* pos_x() { return mem_block_.get(); }
    [[nodiscard]] double* pos_y() { return mem_block_.get() + N_; }
    [[nodiscard]] double* pos_z() { return mem_block_.get() + 2*N_; }
    [[nodiscard]] double* vel_x() { return mem_block_.get() + 3*N_; }
    [[nodiscard]] double* vel_y() { return mem_block_.get() + 4*N_; }
    [[nodiscard]] double* vel_z() { return mem_block_.get() + 5*N_; }
    [[nodiscard]] double* acc_x() { return mem_block_.get() + 6*N_; }
    [[nodiscard]] double* acc_y() { return mem_block_.get() + 7*N_; }
    [[nodiscard]] double* acc_z() { return mem_block_.get() + 8*N_; }
    [[nodiscard]] double* old_acc_x() { return mem_block_.get() + 9*N_; }
    [[nodiscard]] double* old_acc_y() { return mem_block_.get() + 10*N_; }
    [[nodiscard]] double* old_acc_z() { return mem_block_.get() + 11*N_; }
    [[nodiscard]] double* mass() { return mem_block_.get() + 12*N_; }

    // Const raw pointers:
    [[nodiscard]] const double* pos_x() const { return mem_block_.get(); }
    [[nodiscard]] const double* pos_y() const { return mem_block_.get() + N_; }
    [[nodiscard]] const double* pos_z() const { return mem_block_.get() + 2*N_; }
    [[nodiscard]] const double* vel_x() const { return mem_block_.get() + 3*N_; }
    [[nodiscard]] const double* vel_y() const { return mem_block_.get() + 4*N_; }
    [[nodiscard]] const double* vel_z() const { return mem_block_.get() + 5*N_; }
    [[nodiscard]] const double* acc_x() const { return mem_block_.get() + 6*N_; }
    [[nodiscard]] const double* acc_y() const { return mem_block_.get() + 7*N_; }
    [[nodiscard]] const double* acc_z() const { return mem_block_.get() + 8*N_; }
    [[nodiscard]] const double* old_acc_x() const { return mem_block_.get() + 9*N_; }
    [[nodiscard]] const double* old_acc_y() const { return mem_block_.get() + 10*N_; }
    [[nodiscard]] const double* old_acc_z() const { return mem_block_.get() + 11*N_; }
    [[nodiscard]] const double* mass() const { return mem_block_.get() + 12*N_; }
};