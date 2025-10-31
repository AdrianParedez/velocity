#include "velocity/vec4.h"
#include "velocity/common.h"
#include <sstream>
#include <iomanip>

namespace vel {

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

std::string Vec4::toString(int precision) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision);
    oss << "Vec4(" << x << ", " << y << ", " << z << ", " << w << ")";
    return oss.str();
}

} // namespace vel