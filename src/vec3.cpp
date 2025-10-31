#include "velocity/vec3.h"
#include "velocity/common.h"
#include <sstream>
#include <iomanip>

namespace vel {

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

std::string Vec3::toString(int precision) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision);
    oss << "Vec3(" << x << ", " << y << ", " << z << ")";
    return oss.str();
}

} // namespace vel