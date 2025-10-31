#pragma once

#include <cmath>
#include <chrono>

namespace vel {

// Mathematical constants
constexpr float PI = 3.14159265358979323846f;
constexpr float TWO_PI = 2.0f * PI;
constexpr float HALF_PI = PI * 0.5f;
constexpr float EPSILON = 1e-6f;
constexpr float DEG_TO_RAD = PI / 180.0f;
constexpr float RAD_TO_DEG = 180.0f / PI;

// Utility functions
namespace math {

// Angle conversion
constexpr float radians(float degrees) noexcept { return degrees * DEG_TO_RAD; }
constexpr float degrees(float radians) noexcept { return radians * RAD_TO_DEG; }

// Comparison with epsilon
inline bool isEqual(float a, float b, float epsilon = EPSILON) noexcept {
    return std::abs(a - b) <= epsilon;
}

inline bool isZero(float value, float epsilon = EPSILON) noexcept {
    return std::abs(value) <= epsilon;
}

// Clamping and interpolation
constexpr float clamp(float value, float min, float max) noexcept {
    return value < min ? min : (value > max ? max : value);
}

constexpr float lerp(float a, float b, float t) noexcept {
    return a + t * (b - a);
}

inline float smoothstep(float edge0, float edge1, float x) noexcept {
    // Handle edge case where edge0 == edge1 (zero-width transition)
    if (isEqual(edge0, edge1)) {
        return (x >= edge0) ? 1.0f : 0.0f;  // Step function at edge0
    }
    
    // Handle case where edge1 < edge0 (inverted range)
    if (edge1 < edge0) {
        // Swap edges and invert result for inverted range
        float t = clamp((x - edge1) / (edge0 - edge1), 0.0f, 1.0f);
        return 1.0f - (t * t * (3.0f - 2.0f * t));
    }
    
    // Normal case: edge0 < edge1
    float t = clamp((x - edge0) / (edge1 - edge0), 0.0f, 1.0f);
    return t * t * (3.0f - 2.0f * t);
}

// Sign function
constexpr int sign(float value) noexcept {
    return (value > 0.0f) ? 1 : ((value < 0.0f) ? -1 : 0);
}

// Min/Max
constexpr float min(float a, float b) noexcept { return a < b ? a : b; }
constexpr float max(float a, float b) noexcept { return a > b ? a : b; }

// Power of 2 checks
constexpr bool isPowerOfTwo(int value) noexcept {
    return value > 0 && (value & (value - 1)) == 0;
}

// Compile-time flag to control fast inverse square root usage
#ifndef USE_FAST_INV_SQRT
#define USE_FAST_INV_SQRT 1  // Default to enabled, can be overridden at compile time
#endif

// Fast inverse square root approximation with improved precision
inline float fastInverseSqrt(float x) noexcept {
    #if USE_FAST_INV_SQRT
    if (x <= 0.0f) {
        return 0.0f;  // Handle edge case
    }
    
    union { float f; int i; } conv;
    conv.f = x;
    conv.i = 0x5f3759df - (conv.i >> 1);
    
    // First Newton-Raphson iteration
    conv.f *= 1.5f - (x * 0.5f * conv.f * conv.f);
    
    // Second Newton-Raphson iteration for improved precision
    conv.f *= 1.5f - (x * 0.5f * conv.f * conv.f);
    
    return conv.f;
    #else
    // Fallback to standard library implementation
    return (x > 0.0f) ? (1.0f / std::sqrt(x)) : 0.0f;
    #endif
}

// Standard inverse square root for comparison
inline float standardInverseSqrt(float x) noexcept {
    return (x > 0.0f) ? (1.0f / std::sqrt(x)) : 0.0f;
}

// Benchmark function to compare performance (for testing purposes)
inline void benchmarkInverseSqrt(int iterations = 1000000) noexcept {
    // This function is primarily for development/testing
    // In production, the compiler flag should be used to select the best method
    volatile float result = 0.0f;  // Prevent optimization
    
    // Test values
    float test_values[] = {1.0f, 4.0f, 9.0f, 16.0f, 25.0f, 0.25f, 0.01f, 100.0f};
    int num_values = sizeof(test_values) / sizeof(test_values[0]);
    
    // Benchmark fast version
    auto start_fast = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        result += fastInverseSqrt(test_values[i % num_values]);
    }
    auto end_fast = std::chrono::high_resolution_clock::now();
    
    // Benchmark standard version
    auto start_std = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        result += standardInverseSqrt(test_values[i % num_values]);
    }
    auto end_std = std::chrono::high_resolution_clock::now();
    
    // Prevent optimization of result
    (void)result;
}

} // namespace math
} // namespace vel