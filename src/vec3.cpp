#include "velocity/vec3.h"
#include "velocity/common.h"
#include <sstream>
#include <iomanip>

namespace vel {

float& Vec3::operator[](int index) noexcept {
    #ifdef _DEBUG
    if (index < 0 || index >= 3) {
        // In debug mode, return reference to static dummy for invalid index
        static float dummy = 0.0f;
        return dummy;
    }
    #endif
    return (&x)[index];
}

const float& Vec3::operator[](int index) const noexcept {
    #ifdef _DEBUG
    if (index < 0 || index >= 3) {
        // In debug mode, return reference to static zero
        static const float zero = 0.0f;
        return zero;
    }
    #endif
    return (&x)[index];
}

bool Vec3::operator==(const Vec3& rhs) const noexcept {
    return math::isEqual(x, rhs.x) && 
           math::isEqual(y, rhs.y) && 
           math::isEqual(z, rhs.z);
}

Vec3 Vec3::normalized() const noexcept {
    float len = length();
    if (math::isZero(len)) {
        return Vec3::zero();
    }
    return *this / len;
}

Vec3& Vec3::normalize() noexcept {
    float len = length();
    if (!math::isZero(len)) {
        *this /= len;
    }
    return *this;
}

Vec3 Vec3::reflect(const Vec3& normal) const noexcept {
    return *this - normal * (2.0f * dot(normal));
}

Vec3 Vec3::project(const Vec3& onto) const noexcept {
    float dot_product = dot(onto);
    float onto_length_sq = onto.lengthSquared();
    if (math::isZero(onto_length_sq)) {
        return Vec3::zero();
    }
    return onto * (dot_product / onto_length_sq);
}

Vec3 Vec3::operator/(float scalar) const noexcept {
    if (math::isZero(scalar)) {
        #ifdef _DEBUG
        // In debug mode, you could add logging here
        // std::cerr << "Warning: Division by zero in Vec3::operator/, returning zero vector" << std::endl;
        #endif
        return Vec3::zero();
    }
    return Vec3(x / scalar, y / scalar, z / scalar);
}

Vec3& Vec3::operator/=(float scalar) noexcept {
    if (math::isZero(scalar)) {
        #ifdef _DEBUG
        // In debug mode, you could add logging here
        // std::cerr << "Warning: Division by zero in Vec3::operator/=, setting to zero vector" << std::endl;
        #endif
        x = 0.0f;
        y = 0.0f;
        z = 0.0f;
        return *this;
    }
    x /= scalar;
    y /= scalar;
    z /= scalar;
    return *this;
}

std::string Vec3::toString(int precision) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision);
    oss << "Vec3(" << x << ", " << y << ", " << z << ")";
    return oss.str();
}

} // namespace vel