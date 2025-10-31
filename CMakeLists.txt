cmake_minimum_required(VERSION 3.12)
project(velocity VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
endif()

# Include directories
include_directories(include)

# Source files
set(VELOCITY_SOURCES
    src/vec2.cpp
    src/vec3.cpp
    src/vec4.cpp
    src/mat4.cpp
)

# Header files
set(VELOCITY_HEADERS
    include/velocity/common.h
    include/velocity/vec2.h
    include/velocity/vec3.h
    include/velocity/vec4.h
    include/velocity/mat4.h
    include/velocity/velocity.h
)

# Create library
add_library(velocity STATIC ${VELOCITY_SOURCES} ${VELOCITY_HEADERS})

# Set target properties
set_target_properties(velocity PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    PUBLIC_HEADER "${VELOCITY_HEADERS}"
)

# Include directories for the library
target_include_directories(velocity PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Example executable (optional)
option(BUILD_EXAMPLES "Build example programs" ON)
if(BUILD_EXAMPLES)
    add_executable(velocity_example examples/example.cpp)
    target_link_libraries(velocity_example velocity)
endif()

# Tests (optional)
option(BUILD_TESTS "Build test programs" ON)
if(BUILD_TESTS)
    add_executable(velocity_tests tests/test_main.cpp)
    target_link_libraries(velocity_tests velocity)
endif()

# Install configuration
install(TARGETS velocity
    EXPORT velocityTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION include/velocity
)

install(EXPORT velocityTargets
    FILE velocityTargets.cmake
    NAMESPACE velocity::
    DESTINATION lib/cmake/velocity
)

# Create config file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    velocityConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/velocityConfigVersion.cmake"
    DESTINATION lib/cmake/velocity
)