#include "Output.hpp"

Binary_Output::Binary_Output( std::string const &path, 
                              std::vector<std::string> const &names, 
                              std::size_t const num )
: file_{ path, std::ios::binary | std::ios::out }
, num_bodies_{ num } {
    if ( !file_.is_open() ) {
        std::cout << "Failed to open file: " + path + ". Ensure JPL data has been fetched." << std::endl;
    }

    uint64_t const n{ static_cast<uint64_t>( num ) };
    file_.write( reinterpret_cast<char const*>( &n ), sizeof( n ) );

    for ( std::size_t i{}; i < num; ++i ) {
        char name_buf[32] = {};
        std::strncpy( name_buf, names[i].c_str(), 31 );
        file_.write( name_buf, 32 );
    }

    buffer_.reserve( 2 + num * 6 );
}

void Binary_Output::write(
    Particles const &particles,
    std::size_t const num,
    std::size_t const step,
    double const time ) {

    buffer_.clear();

    // Frame header - store step as reinterpreted double for uniform buffer
    uint64_t const s{ step };
    double step_bits;
    std::memcpy( &step_bits, &s, sizeof( double ) );
    buffer_.push_back( step_bits );
    buffer_.push_back( time );

    // Per-body state
    for ( std::size_t i = 0; i < num; ++i ) {
        buffer_.push_back( particles.pos_x()[i] );
        buffer_.push_back( particles.pos_y()[i] );
        buffer_.push_back( particles.pos_z()[i] );
        buffer_.push_back( particles.vel_x()[i] );
        buffer_.push_back( particles.vel_y()[i] );
        buffer_.push_back( particles.vel_z()[i] );
    }

    file_.write( reinterpret_cast<char const*>( buffer_.data() ),
                 buffer_.size() * sizeof( double ) );
}