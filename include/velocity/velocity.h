#pragma once

// Main header that includes all Velocity components
#include "common.h"
#include "vec2.h"
#include "vec3.h"
#include "vec4.h"
#include "mat4.h"

namespace vel {

// Version information
constexpr int VERSION_MAJOR = 1;
constexpr int VERSION_MINOR = 1;
constexpr int VERSION_PATCH = 0;

// Version string helper
inline const char* getVersionString() noexcept {
    return "1.1.0";
}

// Version as single integer for comparison (major * 10000 + minor * 100 + patch)
constexpr int VERSION_NUMBER = VERSION_MAJOR * 10000 + VERSION_MINOR * 100 + VERSION_PATCH;

} // namespace vel