#include "velocity/vec4.h"
#include "velocity/common.h"
#include <sstream>
#include <iomanip>

namespace vel {

float& Vec4::operator[](int index) noexcept {
    #ifdef _DEBUG
    if (index < 0 || index >= 4) {
        // In debug mode, return reference to static dummy for invalid index
        static float dummy = 0.0f;
        return dummy;
    }
    #endif
    return (&x)[index];
}

const float& Vec4::operator[](int index) const noexcept {
    #ifdef _DEBUG
    if (index < 0 || index >= 4) {
        // In debug mode, return reference to static zero
        static const float zero = 0.0f;
        return zero;
    }
    #endif
    return (&x)[index];
}

bool Vec4::operator==(const Vec4& rhs) const noexcept {
    return math::isEqual(x, rhs.x) && 
           math::isEqual(y, rhs.y) && 
           math::isEqual(z, rhs.z) && 
           math::isEqual(w, rhs.w);
}

Vec4 Vec4::normalized() const noexcept {
    float len = length();
    if (math::isZero(len)) {
        return Vec4::zero();
    }
    return *this / len;
}

Vec4& Vec4::normalize() noexcept {
    float len = length();
    if (!math::isZero(len)) {
        *this /= len;
    }
    return *this;
}

Vec4 Vec4::operator/(float scalar) const noexcept {
    if (math::isZero(scalar)) {
        #ifdef _DEBUG
        // In debug mode, you could add logging here
        // std::cerr << "Warning: Division by zero in Vec4::operator/, returning zero vector" << std::endl;
        #endif
        return Vec4::zero();
    }
    return Vec4(x / scalar, y / scalar, z / scalar, w / scalar);
}

Vec4& Vec4::operator/=(float scalar) noexcept {
    if (math::isZero(scalar)) {
        #ifdef _DEBUG
        // In debug mode, you could add logging here
        // std::cerr << "Warning: Division by zero in Vec4::operator/=, setting to zero vector" << std::endl;
        #endif
        x = 0.0f;
        y = 0.0f;
        z = 0.0f;
        w = 0.0f;
        return *this;
    }
    x /= scalar;
    y /= scalar;
    z /= scalar;
    w /= scalar;
    return *this;
}

std::string Vec4::toString(int precision) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision);
    oss << "Vec4(" << x << ", " << y << ", " << z << ", " << w << ")";
    return oss.str();
}

} // namespace vel