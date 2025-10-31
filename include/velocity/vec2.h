#pragma once

#include <string>
#include <cmath>

namespace vel {

class Vec2 {
public:
    float x, y;

    // Constructors
    constexpr Vec2() noexcept : x(0.0f), y(0.0f) {}
    constexpr Vec2(float x, float y) noexcept : x(x), y(y) {}
    constexpr explicit Vec2(float scalar) noexcept : x(scalar), y(scalar) {}

    // Element access
    float& operator[](int index) noexcept;
    const float& operator[](int index) const noexcept;

    // Arithmetic operators
    constexpr Vec2 operator+(const Vec2& rhs) const noexcept { return Vec2(x + rhs.x, y + rhs.y); }
    constexpr Vec2 operator-(const Vec2& rhs) const noexcept { return Vec2(x - rhs.x, y - rhs.y); }
    constexpr Vec2 operator*(float scalar) const noexcept { return Vec2(x * scalar, y * scalar); }
    Vec2 operator/(float scalar) const noexcept;
    constexpr Vec2 operator-() const noexcept { return Vec2(-x, -y); }

    // Compound assignment operators
    constexpr Vec2& operator+=(const Vec2& rhs) noexcept { x += rhs.x; y += rhs.y; return *this; }
    constexpr Vec2& operator-=(const Vec2& rhs) noexcept { x -= rhs.x; y -= rhs.y; return *this; }
    constexpr Vec2& operator*=(float scalar) noexcept { x *= scalar; y *= scalar; return *this; }
    Vec2& operator/=(float scalar) noexcept;

    // Comparison operators
    bool operator==(const Vec2& rhs) const noexcept;
    bool operator!=(const Vec2& rhs) const noexcept { return !(*this == rhs); }

    // Vector operations
    constexpr float dot(const Vec2& rhs) const noexcept { return x * rhs.x + y * rhs.y; }
    constexpr float cross(const Vec2& rhs) const noexcept { return x * rhs.y - y * rhs.x; }
    float length() const noexcept { return std::sqrt(x * x + y * y); }
    constexpr float lengthSquared() const noexcept { return x * x + y * y; }
    Vec2 normalized() const noexcept;
    Vec2& normalize() noexcept;
    float distance(const Vec2& rhs) const noexcept { return (*this - rhs).length(); }
    constexpr float distanceSquared(const Vec2& rhs) const noexcept { return (*this - rhs).lengthSquared(); }
    constexpr Vec2 lerp(const Vec2& rhs, float t) const noexcept { return *this + (rhs - *this) * t; }
    Vec2 reflect(const Vec2& normal) const noexcept;
    float angle() const noexcept { return std::atan2(y, x); }
    Vec2 rotate(float radians) const noexcept;
    Vec2 perpendicular() const noexcept { return Vec2(-y, x); }

    // Utility
    std::string toString(int precision = 3) const;

    // Static constants
    static constexpr Vec2 zero() noexcept { return Vec2(0.0f, 0.0f); }
    static constexpr Vec2 one() noexcept { return Vec2(1.0f, 1.0f); }
    static constexpr Vec2 unitX() noexcept { return Vec2(1.0f, 0.0f); }
    static constexpr Vec2 unitY() noexcept { return Vec2(0.0f, 1.0f); }
};

// Non-member operators
constexpr Vec2 operator*(float scalar, const Vec2& vec) noexcept { return vec * scalar; }

} // namespace vel