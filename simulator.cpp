#include <iostream>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <vector>

static constexpr double G{ 6.67e-11 };

class Vec_3D {
public:
    double x_, y_, z_;
    double norm() {
        return std::sqrt( x_*x_ + y_*y_ + z_*z_ );
    }
};

Vec_3D operator+( Vec_3D const &vec, Vec_3D const &other ) {
    return{ vec.x_ + other.x_, vec.y_ + other.y_, vec.z_ + other.z_ };
}

Vec_3D operator-( Vec_3D const &vec, Vec_3D const &other ) {
    return{ vec.x_ - other.x_, vec.y_ - other.y_, vec.z_ - other.z_ };
}

Vec_3D operator*( Vec_3D const &vec, double const &constant ) {
    return{ vec.x_ * constant, vec.y_ * constant, vec.z_ * constant };
}

Vec_3D &operator+=( Vec_3D &vec, Vec_3D const &other ) {
    vec.x_ += other.x_;
    vec.y_ += other.y_;
    vec.z_ += other.z_;
    return vec;
}

class Body {
private:
    Vec_3D pos_;
    Vec_3D vel_;
    Vec_3D acc_;
    double mass_;

public:
    Body ( Vec_3D pos, Vec_3D vel, double mass ) 
         : pos_( pos ), vel_( vel ), acc_{ 0, 0, 0 }, mass_( mass ) {}
    
    void calculate_new_acc( std::vector<Body> const &other_bodies ) {
        Vec_3D total_force{};

        for ( std::size_t idx{ 0 }; idx < other_bodies.size(); ++idx ) {
            Vec_3D R{ other_bodies[idx].get_pos() - pos_ };
            double dist{ R.norm() };

            if ( dist < 1.0e-10 ) continue;
            double force_mag{ ( G * mass_ * other_bodies[idx].get_mass() ) / ( dist * dist ) };

            total_force += R * ( force_mag / dist );
        }
        acc_ = total_force * ( 1.0 / mass_ );
    }

    void update( double dt ) {
        vel_ += acc_ * dt;
        pos_ += vel_ * dt;
    }

    Vec_3D get_pos() const { return pos_; }
    Vec_3D get_vel() const { return vel_; }
    Vec_3D get_acc() const { return acc_; }
    double get_mass() const { return mass_; }
};

int main() {
    constexpr double MASS{ 1.0e20 };
    constexpr double POS{ 1.0e4 };
    constexpr double dt{ 0.1 };
    int steps{ 0 };
    int current_step{ 0 };

    Body B_1{ { -1.0*POS, POS, POS }, { 0, 0, 0 }, MASS };
    Body B_2{ { POS, -1.0*POS, POS }, { 0, 0, 0 }, MASS };
    Body B_3{ { POS, POS, -1.0*POS }, { 0, 0, 0 }, MASS };

    std::vector<Body> bodies{ B_1, B_2, B_3 };

    std::cout << "<--- 3-Body Simulation --->" << std::endl;

    std::cout << "\nNumber of simulation steps: ";
    std::cin >> steps;

    std::cout << "\nStarting N-Body Simulation..." << std::endl;
    std::cout << std::fixed << std::setprecision( 2 );

    for( steps; steps > 0; --steps ) {
        // Calculates new acceleration.
        for ( std::size_t idx{ 0 }; idx < bodies.size(); ++idx ) {
            bodies[idx].calculate_new_acc( bodies );
        }

        // Updates position and velocity for all bodies.
        for ( std::size_t idx{ 0 }; idx < bodies.size(); ++idx ) {
            bodies[idx].update( dt );
        }

        std::cout << "\n<--- Step: " << ( current_step + 1 ) << " --->" << std::endl;

        // Displays the current position and velocity for all bodies.
        for ( std::size_t idx{ 0 }; idx < bodies.size(); ++idx ) {
            Vec_3D curr_body_pos{ bodies[idx].get_pos() };
            Vec_3D curr_body_vel{ bodies[idx].get_vel() };

            std::cout << "Body " << idx << ": Pos(" << curr_body_pos.x_ / 1000.0 << ", " 
                                                    << curr_body_pos.y_ / 1000.0 << ", " 
                                                    << curr_body_pos.z_ / 1000.0 << ") km, ";

            std::cout << "Vel(" << curr_body_vel.x_ / 3.6 << ", " 
                                                    << curr_body_vel.y_ / 3.6 << ", " 
                                                    << curr_body_vel.z_ / 3.6 << ") km/h, ";
            std::cout << "Speed: " << curr_body_vel.norm() / 3.6 << " km/h" << std::endl;
        }
        std::cout << std::endl;

        // Displays distance form all bodies to the other.
        for ( std::size_t idx{ 0 }; idx < bodies.size(); ++idx ) {
            for ( std::size_t second_idx{ 0 }; second_idx < bodies.size(); ++second_idx ) {
                if ( idx == second_idx ) continue;

                Vec_3D R{ bodies[idx].get_pos() - bodies[second_idx].get_pos() };
                std::cout << "Distance: B" << idx + 1 << "-B"
                                           << second_idx + 1 << ": "
                                           << R.norm() / 1000.0 << " km" << std::endl;
            }
            std::cout << std::endl;
        }
        ++current_step;
    }

    std::cout << "\n<--- End of Simulation --->" << std::endl;

    return 0;
}