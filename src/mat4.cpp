#include "velocity/mat4.h"
#include "velocity/common.h"
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace vel {

Mat4::Mat4() noexcept : m{} {}

Mat4::Mat4(float diagonal) noexcept : m{} {
    m[0] = diagonal;   // m[0][0]
    m[5] = diagonal;   // m[1][1]
    m[10] = diagonal;  // m[2][2]
    m[15] = diagonal;  // m[3][3]
}

Mat4::Mat4(float m00, float m01, float m02, float m03,
           float m10, float m11, float m12, float m13,
           float m20, float m21, float m22, float m23,
           float m30, float m31, float m32, float m33) noexcept {
    m[0] = m00; m[4] = m01; m[8] = m02;  m[12] = m03;
    m[1] = m10; m[5] = m11; m[9] = m12;  m[13] = m13;
    m[2] = m20; m[6] = m21; m[10] = m22; m[14] = m23;
    m[3] = m30; m[7] = m31; m[11] = m32; m[15] = m33;
}

Mat4::Mat4(const std::array<float, 16>& values) noexcept : m(values) {}

Mat4 Mat4::operator+(const Mat4& rhs) const noexcept {
    Mat4 result;
    for (int i = 0; i < 16; ++i) {
        result.m[i] = m[i] + rhs.m[i];
    }
    return result;
}

Mat4 Mat4::operator-(const Mat4& rhs) const noexcept {
    Mat4 result;
    for (int i = 0; i < 16; ++i) {
        result.m[i] = m[i] - rhs.m[i];
    }
    return result;
}

Mat4 Mat4::operator*(const Mat4& rhs) const noexcept {
    Mat4 result;
    for (int col = 0; col < 4; ++col) {
        for (int row = 0; row < 4; ++row) {
            float sum = 0.0f;
            for (int k = 0; k < 4; ++k) {
                sum += (*this)(row, k) * rhs(k, col);
            }
            result(row, col) = sum;
        }
    }
    return result;
}

Mat4 Mat4::operator*(float scalar) const noexcept {
    Mat4 result;
    for (int i = 0; i < 16; ++i) {
        result.m[i] = m[i] * scalar;
    }
    return result;
}

Vec4 Mat4::operator*(const Vec4& vec) const noexcept {
    return Vec4(
        (*this)(0, 0) * vec.x + (*this)(0, 1) * vec.y + (*this)(0, 2) * vec.z + (*this)(0, 3) * vec.w,
        (*this)(1, 0) * vec.x + (*this)(1, 1) * vec.y + (*this)(1, 2) * vec.z + (*this)(1, 3) * vec.w,
        (*this)(2, 0) * vec.x + (*this)(2, 1) * vec.y + (*this)(2, 2) * vec.z + (*this)(2, 3) * vec.w,
        (*this)(3, 0) * vec.x + (*this)(3, 1) * vec.y + (*this)(3, 2) * vec.z + (*this)(3, 3) * vec.w
    );
}

Mat4& Mat4::operator+=(const Mat4& rhs) noexcept {
    for (int i = 0; i < 16; ++i) {
        m[i] += rhs.m[i];
    }
    return *this;
}

Mat4& Mat4::operator-=(const Mat4& rhs) noexcept {
    for (int i = 0; i < 16; ++i) {
        m[i] -= rhs.m[i];
    }
    return *this;
}

Mat4& Mat4::operator*=(const Mat4& rhs) noexcept {
    *this = *this * rhs;
    return *this;
}

Mat4& Mat4::operator*=(float scalar) noexcept {
    for (int i = 0; i < 16; ++i) {
        m[i] *= scalar;
    }
    return *this;
}

bool Mat4::operator==(const Mat4& rhs) const noexcept {
    for (int i = 0; i < 16; ++i) {
        if (!math::isEqual(m[i], rhs.m[i])) {
            return false;
        }
    }
    return true;
}

Mat4 Mat4::transposed() const noexcept {
    Mat4 result;
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            result(col, row) = (*this)(row, col);
        }
    }
    return result;
}

Mat4& Mat4::transpose() noexcept {
    *this = transposed();
    return *this;
}

float Mat4::determinant() const noexcept {
    float det = 0.0f;
    
    // Calculate determinant using cofactor expansion along first row
    det += (*this)(0, 0) * (
        (*this)(1, 1) * ((*this)(2, 2) * (*this)(3, 3) - (*this)(2, 3) * (*this)(3, 2)) -
        (*this)(1, 2) * ((*this)(2, 1) * (*this)(3, 3) - (*this)(2, 3) * (*this)(3, 1)) +
        (*this)(1, 3) * ((*this)(2, 1) * (*this)(3, 2) - (*this)(2, 2) * (*this)(3, 1))
    );
    
    det -= (*this)(0, 1) * (
        (*this)(1, 0) * ((*this)(2, 2) * (*this)(3, 3) - (*this)(2, 3) * (*this)(3, 2)) -
        (*this)(1, 2) * ((*this)(2, 0) * (*this)(3, 3) - (*this)(2, 3) * (*this)(3, 0)) +
        (*this)(1, 3) * ((*this)(2, 0) * (*this)(3, 2) - (*this)(2, 2) * (*this)(3, 0))
    );
    
    det += (*this)(0, 2) * (
        (*this)(1, 0) * ((*this)(2, 1) * (*this)(3, 3) - (*this)(2, 3) * (*this)(3, 1)) -
        (*this)(1, 1) * ((*this)(2, 0) * (*this)(3, 3) - (*this)(2, 3) * (*this)(3, 0)) +
        (*this)(1, 3) * ((*this)(2, 0) * (*this)(3, 1) - (*this)(2, 1) * (*this)(3, 0))
    );
    
    det -= (*this)(0, 3) * (
        (*this)(1, 0) * ((*this)(2, 1) * (*this)(3, 2) - (*this)(2, 2) * (*this)(3, 1)) -
        (*this)(1, 1) * ((*this)(2, 0) * (*this)(3, 2) - (*this)(2, 2) * (*this)(3, 0)) +
        (*this)(1, 2) * ((*this)(2, 0) * (*this)(3, 1) - (*this)(2, 1) * (*this)(3, 0))
    );
    
    return det;
}

Mat4 Mat4::inverse() const {
    float det = determinant();
    if (math::isZero(det)) {
        throw std::runtime_error("Matrix is not invertible (determinant is zero)");
    }
    
    Mat4 result;
    float inv_det = 1.0f / det;
    
    // Calculate adjugate matrix and multiply by inverse determinant
    // This is a simplified version - full implementation would use cofactor matrix
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            // Calculate minor determinant for position (row, col)
            float minor[9];
            int minor_idx = 0;
            
            for (int i = 0; i < 4; ++i) {
                if (i == row) continue;
                for (int j = 0; j < 4; ++j) {
                    if (j == col) continue;
                    minor[minor_idx++] = (*this)(i, j);
                }
            }
            
            float minor_det = minor[0] * (minor[4] * minor[8] - minor[5] * minor[7]) -
                             minor[1] * (minor[3] * minor[8] - minor[5] * minor[6]) +
                             minor[2] * (minor[3] * minor[7] - minor[4] * minor[6]);
            
            result(col, row) = ((row + col) % 2 == 0 ? 1.0f : -1.0f) * minor_det * inv_det;
        }
    }
    
    return result;
}

Mat4& Mat4::invert() {
    *this = inverse();
    return *this;
}

Vec3 Mat4::transformPoint(const Vec3& point) const noexcept {
    Vec4 result = *this * Vec4(point, 1.0f);
    if (math::isZero(result.w)) {
        return Vec3(result.x, result.y, result.z);
    }
    return Vec3(result.x / result.w, result.y / result.w, result.z / result.w);
}

Vec3 Mat4::transformDirection(const Vec3& direction) const noexcept {
    Vec4 result = *this * Vec4(direction, 0.0f);
    return Vec3(result.x, result.y, result.z);
}

Vec3 Mat4::transformNormal(const Vec3& normal) const {
    // Transform normal using inverse transpose
    Mat4 inv_transpose = inverse().transposed();
    Vec4 result = inv_transpose * Vec4(normal, 0.0f);
    return Vec3(result.x, result.y, result.z).normalized();
}

std::string Mat4::toString(int precision) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision);
    oss << "Mat4(\n";
    for (int row = 0; row < 4; ++row) {
        oss << "  [";
        for (int col = 0; col < 4; ++col) {
            oss << (*this)(row, col);
            if (col < 3) oss << ", ";
        }
        oss << "]";
        if (row < 3) oss << ",";
        oss << "\n";
    }
    oss << ")";
    return oss.str();
}

Mat4 Mat4::identity() noexcept {
    return Mat4(1.0f);
}

Mat4 Mat4::translation(const Vec3& translation) noexcept {
    Mat4 result = identity();
    result(0, 3) = translation.x;
    result(1, 3) = translation.y;
    result(2, 3) = translation.z;
    return result;
}

Mat4 Mat4::scale(const Vec3& scale) noexcept {
    Mat4 result;
    result(0, 0) = scale.x;
    result(1, 1) = scale.y;
    result(2, 2) = scale.z;
    result(3, 3) = 1.0f;
    return result;
}

Mat4 Mat4::scale(float uniform_scale) noexcept {
    return scale(Vec3(uniform_scale));
}

Mat4 Mat4::rotationX(float radians) noexcept {
    Mat4 result = identity();
    float cos_r = std::cos(radians);
    float sin_r = std::sin(radians);
    
    result(1, 1) = cos_r;
    result(1, 2) = -sin_r;
    result(2, 1) = sin_r;
    result(2, 2) = cos_r;
    
    return result;
}

Mat4 Mat4::rotationY(float radians) noexcept {
    Mat4 result = identity();
    float cos_r = std::cos(radians);
    float sin_r = std::sin(radians);
    
    result(0, 0) = cos_r;
    result(0, 2) = sin_r;
    result(2, 0) = -sin_r;
    result(2, 2) = cos_r;
    
    return result;
}

Mat4 Mat4::rotationZ(float radians) noexcept {
    Mat4 result = identity();
    float cos_r = std::cos(radians);
    float sin_r = std::sin(radians);
    
    result(0, 0) = cos_r;
    result(0, 1) = -sin_r;
    result(1, 0) = sin_r;
    result(1, 1) = cos_r;
    
    return result;
}

Mat4 Mat4::rotation(const Vec3& axis, float radians) noexcept {
    Vec3 normalized_axis = axis.normalized();
    float cos_r = std::cos(radians);
    float sin_r = std::sin(radians);
    float one_minus_cos = 1.0f - cos_r;
    
    float x = normalized_axis.x;
    float y = normalized_axis.y;
    float z = normalized_axis.z;
    
    Mat4 result;
    result(0, 0) = cos_r + x * x * one_minus_cos;
    result(0, 1) = x * y * one_minus_cos - z * sin_r;
    result(0, 2) = x * z * one_minus_cos + y * sin_r;
    result(0, 3) = 0.0f;
    
    result(1, 0) = y * x * one_minus_cos + z * sin_r;
    result(1, 1) = cos_r + y * y * one_minus_cos;
    result(1, 2) = y * z * one_minus_cos - x * sin_r;
    result(1, 3) = 0.0f;
    
    result(2, 0) = z * x * one_minus_cos - y * sin_r;
    result(2, 1) = z * y * one_minus_cos + x * sin_r;
    result(2, 2) = cos_r + z * z * one_minus_cos;
    result(2, 3) = 0.0f;
    
    result(3, 0) = 0.0f;
    result(3, 1) = 0.0f;
    result(3, 2) = 0.0f;
    result(3, 3) = 1.0f;
    
    return result;
}

Mat4 Mat4::lookAt(const Vec3& eye, const Vec3& center, const Vec3& up) noexcept {
    Vec3 f = (center - eye).normalized();
    Vec3 s = f.cross(up).normalized();
    Vec3 u = s.cross(f);
    
    Mat4 result = identity();
    result(0, 0) = s.x;
    result(1, 0) = s.y;
    result(2, 0) = s.z;
    result(0, 1) = u.x;
    result(1, 1) = u.y;
    result(2, 1) = u.z;
    result(0, 2) = -f.x;
    result(1, 2) = -f.y;
    result(2, 2) = -f.z;
    result(3, 0) = -s.dot(eye);
    result(3, 1) = -u.dot(eye);
    result(3, 2) = f.dot(eye);
    
    return result;
}

Mat4 Mat4::perspective(float fovy, float aspect, float near, float far) noexcept {
    float tan_half_fovy = std::tan(fovy * 0.5f);
    
    Mat4 result;
    result(0, 0) = 1.0f / (aspect * tan_half_fovy);
    result(1, 1) = 1.0f / tan_half_fovy;
    result(2, 2) = -(far + near) / (far - near);
    result(2, 3) = -1.0f;
    result(3, 2) = -(2.0f * far * near) / (far - near);
    
    return result;
}

Mat4 Mat4::ortho(float left, float right, float bottom, float top, float near, float far) noexcept {
    Mat4 result = identity();
    result(0, 0) = 2.0f / (right - left);
    result(1, 1) = 2.0f / (top - bottom);
    result(2, 2) = -2.0f / (far - near);
    result(3, 0) = -(right + left) / (right - left);
    result(3, 1) = -(top + bottom) / (top - bottom);
    result(3, 2) = -(far + near) / (far - near);
    
    return result;
}

Mat4 Mat4::frustum(float left, float right, float bottom, float top, float near, float far) noexcept {
    Mat4 result;
    result(0, 0) = (2.0f * near) / (right - left);
    result(1, 1) = (2.0f * near) / (top - bottom);
    result(2, 0) = (right + left) / (right - left);
    result(2, 1) = (top + bottom) / (top - bottom);
    result(2, 2) = -(far + near) / (far - near);
    result(2, 3) = -1.0f;
    result(3, 2) = -(2.0f * far * near) / (far - near);
    
    return result;
}

Vec3 Mat4::getTranslation() const noexcept {
    return Vec3((*this)(0, 3), (*this)(1, 3), (*this)(2, 3));
}

Vec3 Mat4::getScale() const noexcept {
    Vec3 scale_x((*this)(0, 0), (*this)(1, 0), (*this)(2, 0));
    Vec3 scale_y((*this)(0, 1), (*this)(1, 1), (*this)(2, 1));
    Vec3 scale_z((*this)(0, 2), (*this)(1, 2), (*this)(2, 2));
    
    return Vec3(scale_x.length(), scale_y.length(), scale_z.length());
}

Mat4 Mat4::getRotation() const noexcept {
    Vec3 scale = getScale();
    Mat4 result = *this;
    
    // Remove scale from rotation matrix
    result(0, 0) /= scale.x; result(0, 1) /= scale.y; result(0, 2) /= scale.z;
    result(1, 0) /= scale.x; result(1, 1) /= scale.y; result(1, 2) /= scale.z;
    result(2, 0) /= scale.x; result(2, 1) /= scale.y; result(2, 2) /= scale.z;
    
    // Remove translation
    result(0, 3) = 0.0f;
    result(1, 3) = 0.0f;
    result(2, 3) = 0.0f;
    result(3, 3) = 1.0f;
    
    return result;
}

Mat4 operator*(float scalar, const Mat4& mat) noexcept {
    return mat * scalar;
}

} // namespace vel