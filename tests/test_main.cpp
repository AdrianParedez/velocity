#include <velocity/velocity.h>
#include <iostream>
#include <cassert>
#include <cmath>

using namespace vel;

void testVec2() {
    std::cout << "Testing Vec2..." << std::endl;
    
    // Constructor tests
    Vec2 v1;
    assert(v1.x == 0.0f && v1.y == 0.0f);
    
    Vec2 v2(3.0f, 4.0f);
    assert(v2.x == 3.0f && v2.y == 4.0f);
    
    Vec2 v3(5.0f);
    assert(v3.x == 5.0f && v3.y == 5.0f);
    
    // Arithmetic operations
    Vec2 sum = v2 + Vec2(1.0f, 2.0f);
    assert(sum.x == 4.0f && sum.y == 6.0f);
    
    Vec2 diff = v2 - Vec2(1.0f, 2.0f);
    assert(diff.x == 2.0f && diff.y == 2.0f);
    
    Vec2 scaled = v2 * 2.0f;
    assert(scaled.x == 6.0f && scaled.y == 8.0f);
    
    // Vector operations
    assert(std::abs(v2.length() - 5.0f) < EPSILON);
    assert(v2.lengthSquared() == 25.0f);
    
    Vec2 normalized = v2.normalized();
    assert(std::abs(normalized.length() - 1.0f) < EPSILON);
    
    // Dot product
    assert(v2.dot(Vec2(1.0f, 0.0f)) == 3.0f);
    
    // Cross product (2D)
    assert(v2.cross(Vec2(1.0f, 0.0f)) == -4.0f);
    
    // Static constants
    assert(Vec2::zero() == Vec2(0.0f, 0.0f));
    assert(Vec2::one() == Vec2(1.0f, 1.0f));
    assert(Vec2::unitX() == Vec2(1.0f, 0.0f));
    assert(Vec2::unitY() == Vec2(0.0f, 1.0f));
    
    std::cout << "Vec2 tests passed!" << std::endl;
}

void testVec3() {
    std::cout << "Testing Vec3..." << std::endl;
    
    // Constructor tests
    Vec3 v1;
    assert(v1.x == 0.0f && v1.y == 0.0f && v1.z == 0.0f);
    
    Vec3 v2(1.0f, 2.0f, 3.0f);
    assert(v2.x == 1.0f && v2.y == 2.0f && v2.z == 3.0f);
    
    Vec3 v3(5.0f);
    assert(v3.x == 5.0f && v3.y == 5.0f && v3.z == 5.0f);
    
    Vec3 v4(Vec2(1.0f, 2.0f), 3.0f);
    assert(v4.x == 1.0f && v4.y == 2.0f && v4.z == 3.0f);
    
    // Arithmetic operations
    Vec3 sum = v2 + Vec3(1.0f, 1.0f, 1.0f);
    assert(sum.x == 2.0f && sum.y == 3.0f && sum.z == 4.0f);
    
    // Vector operations
    assert(std::abs(v2.length() - std::sqrt(14.0f)) < EPSILON);
    assert(v2.lengthSquared() == 14.0f);
    
    // Cross product
    Vec3 cross = Vec3::unitX().cross(Vec3::unitY());
    assert(cross == Vec3::unitZ());
    
    // Dot product
    assert(v2.dot(Vec3(1.0f, 0.0f, 0.0f)) == 1.0f);
    
    // Swizzling
    Vec2 xy = v2.xy();
    assert(xy.x == 1.0f && xy.y == 2.0f);
    
    // Static constants
    assert(Vec3::zero() == Vec3(0.0f, 0.0f, 0.0f));
    assert(Vec3::up() == Vec3(0.0f, 1.0f, 0.0f));
    assert(Vec3::forward() == Vec3(0.0f, 0.0f, -1.0f));
    
    std::cout << "Vec3 tests passed!" << std::endl;
}

void testVec4() {
    std::cout << "Testing Vec4..." << std::endl;
    
    // Constructor tests
    Vec4 v1;
    assert(v1.x == 0.0f && v1.y == 0.0f && v1.z == 0.0f && v1.w == 0.0f);
    
    Vec4 v2(1.0f, 2.0f, 3.0f, 4.0f);
    assert(v2.x == 1.0f && v2.y == 2.0f && v2.z == 3.0f && v2.w == 4.0f);
    
    Vec4 v3(Vec3(1.0f, 2.0f, 3.0f), 4.0f);
    assert(v3.x == 1.0f && v3.y == 2.0f && v3.z == 3.0f && v3.w == 4.0f);
    
    // Vector operations
    assert(v2.lengthSquared() == 30.0f);
    assert(v2.dot(Vec4(1.0f, 0.0f, 0.0f, 0.0f)) == 1.0f);
    
    // Swizzling
    Vec3 xyz = v2.xyz();
    assert(xyz.x == 1.0f && xyz.y == 2.0f && xyz.z == 3.0f);
    
    std::cout << "Vec4 tests passed!" << std::endl;
}

void testMat4() {
    std::cout << "Testing Mat4..." << std::endl;
    
    // Identity matrix
    Mat4 identity = Mat4::identity();
    assert(identity(0, 0) == 1.0f && identity(1, 1) == 1.0f && 
           identity(2, 2) == 1.0f && identity(3, 3) == 1.0f);
    assert(identity(0, 1) == 0.0f && identity(1, 0) == 0.0f);
    
    // Translation matrix
    Mat4 translation = Mat4::translation(Vec3(1.0f, 2.0f, 3.0f));
    Vec3 point(0.0f, 0.0f, 0.0f);
    Vec3 translated = translation.transformPoint(point);
    assert(std::abs(translated.x - 1.0f) < EPSILON &&
           std::abs(translated.y - 2.0f) < EPSILON &&
           std::abs(translated.z - 3.0f) < EPSILON);
    
    // Scale matrix
    Mat4 scale = Mat4::scale(Vec3(2.0f, 3.0f, 4.0f));
    Vec3 scaled_point = scale.transformPoint(Vec3(1.0f, 1.0f, 1.0f));
    assert(std::abs(scaled_point.x - 2.0f) < EPSILON &&
           std::abs(scaled_point.y - 3.0f) < EPSILON &&
           std::abs(scaled_point.z - 4.0f) < EPSILON);
    
    // Matrix multiplication
    Mat4 combined = translation * scale;
    Vec3 combined_result = combined.transformPoint(Vec3(1.0f, 1.0f, 1.0f));
    assert(std::abs(combined_result.x - 3.0f) < EPSILON &&
           std::abs(combined_result.y - 5.0f) < EPSILON &&
           std::abs(combined_result.z - 7.0f) < EPSILON);
    
    // Determinant
    assert(std::abs(identity.determinant() - 1.0f) < EPSILON);
    
    // Transpose
    Mat4 test_matrix(1, 2, 3, 4,
                     5, 6, 7, 8,
                     9, 10, 11, 12,
                     13, 14, 15, 16);
    Mat4 transposed = test_matrix.transposed();
    assert(transposed(0, 1) == 5.0f && transposed(1, 0) == 2.0f);
    
    std::cout << "Mat4 tests passed!" << std::endl;
}

void testMathUtils() {
    std::cout << "Testing Math Utils..." << std::endl;
    
    // Angle conversion
    assert(std::abs(math::radians(180.0f) - PI) < EPSILON);
    assert(std::abs(math::degrees(PI) - 180.0f) < EPSILON);
    
    // Comparison
    assert(math::isEqual(1.0f, 1.0f + EPSILON * 0.5f));
    assert(!math::isEqual(1.0f, 1.0f + EPSILON * 2.0f));
    assert(math::isZero(EPSILON * 0.5f));
    
    // Clamping
    assert(math::clamp(1.5f, 0.0f, 1.0f) == 1.0f);
    assert(math::clamp(-0.5f, 0.0f, 1.0f) == 0.0f);
    assert(math::clamp(0.5f, 0.0f, 1.0f) == 0.5f);
    
    // Interpolation
    assert(math::lerp(0.0f, 10.0f, 0.5f) == 5.0f);
    assert(math::lerp(0.0f, 10.0f, 0.0f) == 0.0f);
    assert(math::lerp(0.0f, 10.0f, 1.0f) == 10.0f);
    
    // Sign function
    assert(math::sign(5.0f) == 1);
    assert(math::sign(-3.0f) == -1);
    assert(math::sign(0.0f) == 0);
    
    // Power of two
    assert(math::isPowerOfTwo(1));
    assert(math::isPowerOfTwo(2));
    assert(math::isPowerOfTwo(16));
    assert(!math::isPowerOfTwo(3));
    assert(!math::isPowerOfTwo(15));
    
    std::cout << "Math Utils tests passed!" << std::endl;
}

int main() {
    std::cout << "Running Velocity Math Library Tests" << std::endl;
    std::cout << "===========================" << std::endl;
    
    try {
        testVec2();
        testVec3();
        testVec4();
        testMat4();
        testMathUtils();
        
        std::cout << std::endl << "All tests passed successfully!" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}