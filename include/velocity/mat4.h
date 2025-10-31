#pragma once

#include "vec3.h"
#include "vec4.h"
#include <array>
#include <string>

namespace vel {

/**
 * 4x4 Matrix class with column-major storage (OpenGL/GLM compatible)
 * 
 * Storage layout:
 * - Data is stored in column-major order: [col0, col1, col2, col3]
 * - Each column contains 4 elements: [row0, row1, row2, row3]
 * - Memory layout: [c0r0, c0r1, c0r2, c0r3, c1r0, c1r1, c1r2, c1r3, ...]
 * 
 * Constructor usage:
 * - Mat4(c0r0, c0r1, c0r2, c0r3, c1r0, ...) - column-major parameters
 * - Mat4::fromRows(r0c0, r0c1, r0c2, r0c3, r1c0, ...) - row-major parameters
 * 
 * This layout is compatible with OpenGL and GLM libraries.
 */
class Mat4 {
private:
    // Column-major storage: m[col * 4 + row]
    std::array<float, 16> m;

public:
    // Constructors
    Mat4() noexcept;
    explicit Mat4(float diagonal) noexcept;
    
    // Column-major constructor: parameters are organized by columns
    // Mat4(col0_x, col0_y, col0_z, col0_w,    // First column
    //      col1_x, col1_y, col1_z, col1_w,    // Second column  
    //      col2_x, col2_y, col2_z, col2_w,    // Third column
    //      col3_x, col3_y, col3_z, col3_w)    // Fourth column
    Mat4(float c0r0, float c0r1, float c0r2, float c0r3,
         float c1r0, float c1r1, float c1r2, float c1r3,
         float c2r0, float c2r1, float c2r2, float c2r3,
         float c3r0, float c3r1, float c3r2, float c3r3) noexcept;
    
    explicit Mat4(const std::array<float, 16>& values) noexcept;

    // Element access (row, col)
    constexpr float& operator()(int row, int col) noexcept { return m[col * 4 + row]; }
    constexpr const float& operator()(int row, int col) const noexcept { return m[col * 4 + row]; }
    
    // Raw data access
    constexpr float* data() noexcept { return m.data(); }
    constexpr const float* data() const noexcept { return m.data(); }

    // Arithmetic operators
    Mat4 operator+(const Mat4& rhs) const noexcept;
    Mat4 operator-(const Mat4& rhs) const noexcept;
    Mat4 operator*(const Mat4& rhs) const noexcept;
    Mat4 operator*(float scalar) const noexcept;
    Vec4 operator*(const Vec4& vec) const noexcept;

    // Compound assignment operators
    Mat4& operator+=(const Mat4& rhs) noexcept;
    Mat4& operator-=(const Mat4& rhs) noexcept;
    Mat4& operator*=(const Mat4& rhs) noexcept;
    Mat4& operator*=(float scalar) noexcept;

    // Comparison operators
    bool operator==(const Mat4& rhs) const noexcept;
    bool operator!=(const Mat4& rhs) const noexcept { return !(*this == rhs); }

    // Matrix operations
    Mat4 transposed() const noexcept;
    Mat4& transpose() noexcept;
    float determinant() const noexcept;
    Mat4 inverse() const;
    Mat4& invert();

    // Transform operations
    Vec3 transformPoint(const Vec3& point) const noexcept;
    Vec3 transformDirection(const Vec3& direction) const noexcept;
    Vec3 transformNormal(const Vec3& normal) const;

    // Utility
    std::string toString(int precision = 3) const;

    // Static factory methods
    static Mat4 identity() noexcept;
    
    // Row-major constructor for convenience (parameters organized by rows)
    // This is often more intuitive when writing matrices by hand
    static Mat4 fromRows(float r0c0, float r0c1, float r0c2, float r0c3,
                         float r1c0, float r1c1, float r1c2, float r1c3,
                         float r2c0, float r2c1, float r2c2, float r2c3,
                         float r3c0, float r3c1, float r3c2, float r3c3) noexcept;
    
    static Mat4 translation(const Vec3& translation) noexcept;
    static Mat4 scale(const Vec3& scale) noexcept;
    static Mat4 scale(float uniform_scale) noexcept;
    static Mat4 rotationX(float radians) noexcept;
    static Mat4 rotationY(float radians) noexcept;
    static Mat4 rotationZ(float radians) noexcept;
    static Mat4 rotation(const Vec3& axis, float radians) noexcept;
    static Mat4 lookAt(const Vec3& eye, const Vec3& center, const Vec3& up) noexcept;
    static Mat4 perspective(float fovy, float aspect, float near, float far) noexcept;
    static Mat4 ortho(float left, float right, float bottom, float top, float near, float far) noexcept;
    static Mat4 frustum(float left, float right, float bottom, float top, float near, float far) noexcept;

    // Decomposition
    Vec3 getTranslation() const noexcept;
    Vec3 getScale() const noexcept;
    Mat4 getRotation() const noexcept;
};

// Non-member operators
Mat4 operator*(float scalar, const Mat4& mat) noexcept;

} // namespace vel