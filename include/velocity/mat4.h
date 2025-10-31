#pragma once

#include "vec3.h"
#include "vec4.h"
#include <array>
#include <string>

namespace vel {

class Mat4 {
private:
    // Column-major storage: m[col][row] or m[col * 4 + row]
    std::array<float, 16> m;

public:
    // Constructors
    Mat4() noexcept;
    explicit Mat4(float diagonal) noexcept;
    Mat4(float m00, float m01, float m02, float m03,
         float m10, float m11, float m12, float m13,
         float m20, float m21, float m22, float m23,
         float m30, float m31, float m32, float m33) noexcept;
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