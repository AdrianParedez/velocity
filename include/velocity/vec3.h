#pragma once

#include "vec2.h"
#include <string>
#include <cmath>

namespace vel {

class Vec3 {
public:
    float x, y, z;

    // Constructors
    constexpr Vec3() noexcept : x(0.0f), y(0.0f), z(0.0f) {}
    constexpr Vec3(float x, float y, float z) noexcept : x(x), y(y), z(z) {}
    constexpr explicit Vec3(float scalar) noexcept : x(scalar), y(scalar), z(scalar) {}
    constexpr Vec3(const Vec2& xy, float z) noexcept : x(xy.x), y(xy.y), z(z) {}

    // Element access
    float& operator[](int index) noexcept;
    const float& operator[](int index) const noexcept;

    // Arithmetic operators
    constexpr Vec3 operator+(const Vec3& rhs) const noexcept { return Vec3(x + rhs.x, y + rhs.y, z + rhs.z); }
    constexpr Vec3 operator-(const Vec3& rhs) const noexcept { return Vec3(x - rhs.x, y - rhs.y, z - rhs.z); }
    constexpr Vec3 operator*(float scalar) const noexcept { return Vec3(x * scalar, y * scalar, z * scalar); }
    Vec3 operator/(float scalar) const noexcept;
    constexpr Vec3 operator-() const noexcept { return Vec3(-x, -y, -z); }

    // Compound assignment operators
    constexpr Vec3& operator+=(const Vec3& rhs) noexcept { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
    constexpr Vec3& operator-=(const Vec3& rhs) noexcept { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
    constexpr Vec3& operator*=(float scalar) noexcept { x *= scalar; y *= scalar; z *= scalar; return *this; }
    Vec3& operator/=(float scalar) noexcept;

    // Comparison operators
    bool operator==(const Vec3& rhs) const noexcept;
    bool operator!=(const Vec3& rhs) const noexcept { return !(*this == rhs); }

    // Vector operations
    constexpr float dot(const Vec3& rhs) const noexcept { return x * rhs.x + y * rhs.y + z * rhs.z; }
    constexpr Vec3 cross(const Vec3& rhs) const noexcept {
        return Vec3(y * rhs.z - z * rhs.y, z * rhs.x - x * rhs.z, x * rhs.y - y * rhs.x);
    }
    float length() const noexcept { return std::sqrt(x * x + y * y + z * z); }
    constexpr float lengthSquared() const noexcept { return x * x + y * y + z * z; }
    Vec3 normalized() const noexcept;
    Vec3& normalize() noexcept;
    float distance(const Vec3& rhs) const noexcept { return (*this - rhs).length(); }
    constexpr float distanceSquared(const Vec3& rhs) const noexcept { return (*this - rhs).lengthSquared(); }
    constexpr Vec3 lerp(const Vec3& rhs, float t) const noexcept { return *this + (rhs - *this) * t; }
    Vec3 reflect(const Vec3& normal) const noexcept;
    Vec3 project(const Vec3& onto) const noexcept;
    Vec3 reject(const Vec3& onto) const noexcept { return *this - project(onto); }

    // Swizzling
    constexpr Vec2 xy() const noexcept { return Vec2(x, y); }
    constexpr Vec2 xz() const noexcept { return Vec2(x, z); }
    constexpr Vec2 yz() const noexcept { return Vec2(y, z); }

    // Utility
    std::string toString(int precision = 3) const;

    // Static constants
    static constexpr Vec3 zero() noexcept { return Vec3(0.0f, 0.0f, 0.0f); }
    static constexpr Vec3 one() noexcept { return Vec3(1.0f, 1.0f, 1.0f); }
    static constexpr Vec3 unitX() noexcept { return Vec3(1.0f, 0.0f, 0.0f); }
    static constexpr Vec3 unitY() noexcept { return Vec3(0.0f, 1.0f, 0.0f); }
    static constexpr Vec3 unitZ() noexcept { return Vec3(0.0f, 0.0f, 1.0f); }
    static constexpr Vec3 forward() noexcept { return Vec3(0.0f, 0.0f, -1.0f); }
    static constexpr Vec3 back() noexcept { return Vec3(0.0f, 0.0f, 1.0f); }
    static constexpr Vec3 up() noexcept { return Vec3(0.0f, 1.0f, 0.0f); }
    static constexpr Vec3 down() noexcept { return Vec3(0.0f, -1.0f, 0.0f); }
    static constexpr Vec3 right() noexcept { return Vec3(1.0f, 0.0f, 0.0f); }
    static constexpr Vec3 left() noexcept { return Vec3(-1.0f, 0.0f, 0.0f); }
};

// Non-member operators
constexpr Vec3 operator*(float scalar, const Vec3& vec) noexcept { return vec * scalar; }

} // namespace vel