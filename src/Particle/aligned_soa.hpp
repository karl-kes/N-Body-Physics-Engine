#pragma once

#include <algorithm>
#include <cstddef>
#include <memory>
#include <cstdlib>

#if defined(_WIN32)
    #include <malloc.h>
    inline void* AlignedAlloc( std::size_t alignment, std::size_t size ) {
        return _aligned_malloc( size, alignment );
    }
    inline void AlignedFree( void* ptr ) { _aligned_free( ptr ); }
#elif defined(__APPLE__)
    inline void* AlignedAlloc( std::size_t alignment, std::size_t size ) {
        // posix_memalign requires alignment >= sizeof(void*) and a power of two.
        alignment = std::max( alignment, sizeof( void* ) );
        void* ptr{};
        posix_memalign( &ptr, alignment, size );
        return ptr;
    }
    inline void AlignedFree( void* ptr ) { std::free( ptr ); }
#else
    inline void* AlignedAlloc( std::size_t alignment, std::size_t size ) {
        return std::aligned_alloc( alignment, size );
    }
    inline void AlignedFree( void* ptr ) { std::free( ptr ); }
#endif

struct AlignedDeleter {
    template <typename T> void operator()( T* ptr ) const { AlignedFree( ptr ); }
};

#if defined(__AVX512F__)
    constexpr std::size_t SIMD_BYTES{ 64 };
#elif defined(__AVX2__) || defined(__AVX__)
    constexpr std::size_t SIMD_BYTES{ 32 };
#elif defined(__SSE2__)
    constexpr std::size_t SIMD_BYTES{ 16 };
#else
    constexpr std::size_t SIMD_BYTES{ sizeof( double ) };
#endif

template <typename T> class AlignedSoA {
private:
    static constexpr std::size_t alignment_bytes{ SIMD_BYTES };
    static constexpr std::size_t elements_per_alignment{ SIMD_BYTES / sizeof( T ) };

    std::size_t num_elements_;
    std::size_t stride_length_;
    std::size_t num_arrays_;
    std::unique_ptr<T[], AlignedDeleter> memory_block_;

public:
    // Round N up to the nearest multiple of the SIMD-element width so every
    // sub-array start is naturally aligned for vectorized loads and stores.
    [[nodiscard]] static constexpr std::size_t round_up( std::size_t const unpadded ) {
        return ( unpadded + elements_per_alignment - 1 ) & ~( elements_per_alignment - 1 );
    }

    AlignedSoA()
    : num_elements_{}
    , stride_length_{}
    , num_arrays_{}
    , memory_block_{ nullptr }
    { }

    AlignedSoA( AlignedSoA&& ) noexcept = default;
    AlignedSoA& operator=( AlignedSoA&& ) noexcept = default;

    AlignedSoA( std::size_t const num_elements, std::size_t const num_arrays )
    : num_elements_{ num_elements }
    , stride_length_{ round_up( num_elements ) }
    , num_arrays_{ num_arrays } {
        std::size_t const total_elements{ num_arrays_ * stride_length_ };
        std::size_t const total_bytes{ total_elements * sizeof( T ) };

        T* ptr{ static_cast<T*>( AlignedAlloc( alignment_bytes, total_bytes ) ) };
        if ( !ptr ) throw std::bad_alloc();

        std::fill_n( ptr, total_elements, T{} );
        memory_block_.reset( ptr );
    }

    [[nodiscard]] std::size_t stride() const { return stride_length_; }
    [[nodiscard]] std::size_t num_elements() const { return num_elements_; }

    T* operator[]( std::size_t const array_index ) {
        return memory_block_.get() + array_index * stride();
    }
    T const* operator[]( std::size_t const array_index ) const {
        return memory_block_.get() + array_index * stride();
    }
};
