#pragma once

#include "vec3.h"
#include <string>
#include <cmath>

namespace vel {

class Vec4 {
public:
    float x, y, z, w;

    // Constructors
    constexpr Vec4() noexcept : x(0.0f), y(0.0f), z(0.0f), w(0.0f) {}
    constexpr Vec4(float x, float y, float z, float w) noexcept : x(x), y(y), z(z), w(w) {}
    constexpr explicit Vec4(float scalar) noexcept : x(scalar), y(scalar), z(scalar), w(scalar) {}
    constexpr Vec4(const Vec3& xyz, float w) noexcept : x(xyz.x), y(xyz.y), z(xyz.z), w(w) {}
    constexpr Vec4(const Vec2& xy, float z, float w) noexcept : x(xy.x), y(xy.y), z(z), w(w) {}

    // Element access
    float& operator[](int index) noexcept;
    const float& operator[](int index) const noexcept;

    // Arithmetic operators
    constexpr Vec4 operator+(const Vec4& rhs) const noexcept { return Vec4(x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.w); }
    constexpr Vec4 operator-(const Vec4& rhs) const noexcept { return Vec4(x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w); }
    constexpr Vec4 operator*(float scalar) const noexcept { return Vec4(x * scalar, y * scalar, z * scalar, w * scalar); }
    Vec4 operator/(float scalar) const noexcept;
    constexpr Vec4 operator-() const noexcept { return Vec4(-x, -y, -z, -w); }

    // Compound assignment operators
    constexpr Vec4& operator+=(const Vec4& rhs) noexcept { x += rhs.x; y += rhs.y; z += rhs.z; w += rhs.w; return *this; }
    constexpr Vec4& operator-=(const Vec4& rhs) noexcept { x -= rhs.x; y -= rhs.y; z -= rhs.z; w -= rhs.w; return *this; }
    constexpr Vec4& operator*=(float scalar) noexcept { x *= scalar; y *= scalar; z *= scalar; w *= scalar; return *this; }
    Vec4& operator/=(float scalar) noexcept;

    // Comparison operators
    bool operator==(const Vec4& rhs) const noexcept;
    bool operator!=(const Vec4& rhs) const noexcept { return !(*this == rhs); }

    // Vector operations
    constexpr float dot(const Vec4& rhs) const noexcept { return x * rhs.x + y * rhs.y + z * rhs.z + w * rhs.w; }
    float length() const noexcept { return std::sqrt(x * x + y * y + z * z + w * w); }
    constexpr float lengthSquared() const noexcept { return x * x + y * y + z * z + w * w; }
    Vec4 normalized() const noexcept;
    Vec4& normalize() noexcept;
    constexpr Vec4 lerp(const Vec4& rhs, float t) const noexcept { return *this + (rhs - *this) * t; }

    // Swizzling
    constexpr Vec2 xy() const noexcept { return Vec2(x, y); }
    constexpr Vec2 zw() const noexcept { return Vec2(z, w); }
    constexpr Vec3 xyz() const noexcept { return Vec3(x, y, z); }

    // Utility
    std::string toString(int precision = 3) const;

    // Static constants
    static constexpr Vec4 zero() noexcept { return Vec4(0.0f, 0.0f, 0.0f, 0.0f); }
    static constexpr Vec4 one() noexcept { return Vec4(1.0f, 1.0f, 1.0f, 1.0f); }
    static constexpr Vec4 unitX() noexcept { return Vec4(1.0f, 0.0f, 0.0f, 0.0f); }
    static constexpr Vec4 unitY() noexcept { return Vec4(0.0f, 1.0f, 0.0f, 0.0f); }
    static constexpr Vec4 unitZ() noexcept { return Vec4(0.0f, 0.0f, 1.0f, 0.0f); }
    static constexpr Vec4 unitW() noexcept { return Vec4(0.0f, 0.0f, 0.0f, 1.0f); }
};

// Non-member operators
constexpr Vec4 operator*(float scalar, const Vec4& vec) noexcept { return vec * scalar; }

} // namespace vel