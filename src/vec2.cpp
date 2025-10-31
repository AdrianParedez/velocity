#include "velocity/vec2.h"
#include "velocity/common.h"
#include <sstream>
#include <iomanip>

namespace vel {

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