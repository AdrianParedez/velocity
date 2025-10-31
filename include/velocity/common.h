#pragma once

#include <cmath>

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

constexpr float smoothstep(float edge0, float edge1, float x) noexcept {
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

// Fast inverse square root approximation
inline float fastInverseSqrt(float x) noexcept {
    union { float f; int i; } conv;
    conv.f = x;
    conv.i = 0x5f3759df - (conv.i >> 1);
    conv.f *= 1.5f - (x * 0.5f * conv.f * conv.f);
    return conv.f;
}

} // namespace math
} // namespace vel