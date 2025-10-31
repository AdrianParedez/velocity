#include "velocity/vec2.h"
#include "velocity/common.h"
#include <sstream>
#include <iomanip>

namespace vel {

float& Vec2::operator[](int index) noexcept {
    #ifdef _DEBUG
    if (index < 0 || index >= 2) {
        // In debug mode, assert on invalid index
        // Could also log warning here
        static float dummy = 0.0f;
        return dummy;
    }
    #endif
    return (&x)[index];
}

const float& Vec2::operator[](int index) const noexcept {
    #ifdef _DEBUG
    if (index < 0 || index >= 2) {
        // In debug mode, return reference to static zero
        static const float zero = 0.0f;
        return zero;
    }
    #endif
    return (&x)[index];
}

bool Vec2::operator==(const Vec2& rhs) const noexcept {
    return math::isEqual(x, rhs.x) && math::isEqual(y, rhs.y);
}

Vec2 Vec2::normalized() const noexcept {
    float len = length();
    if (math::isZero(len)) {
        return Vec2::zero();
    }
    return *this / len;
}

Vec2& Vec2::normalize() noexcept {
    float len = length();
    if (!math::isZero(len)) {
        *this /= len;
    }
    return *this;
}

Vec2 Vec2::reflect(const Vec2& normal) const noexcept {
    return *this - normal * (2.0f * dot(normal));
}

Vec2 Vec2::operator/(float scalar) const noexcept {
    if (math::isZero(scalar)) {
        #ifdef _DEBUG
        // In debug mode, you could add logging here
        // std::cerr << "Warning: Division by zero in Vec2::operator/, returning zero vector" << std::endl;
        #endif
        return Vec2::zero();
    }
    return Vec2(x / scalar, y / scalar);
}

Vec2& Vec2::operator/=(float scalar) noexcept {
    if (math::isZero(scalar)) {
        #ifdef _DEBUG
        // In debug mode, you could add logging here
        // std::cerr << "Warning: Division by zero in Vec2::operator/=, setting to zero vector" << std::endl;
        #endif
        x = 0.0f;
        y = 0.0f;
        return *this;
    }
    x /= scalar;
    y /= scalar;
    return *this;
}

Vec2 Vec2::rotate(float radians) const noexcept {
    float cos_r = std::cos(radians);
    float sin_r = std::sin(radians);
    return Vec2(x * cos_r - y * sin_r, x * sin_r + y * cos_r);
}

std::string Vec2::toString(int precision) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision);
    oss << "Vec2(" << x << ", " << y << ")";
    return oss.str();
}

} // namespace vel