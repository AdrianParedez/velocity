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
    
    // Division by zero protection tests
    Vec2 div_test(3.0f, 4.0f);
    Vec2 div_result = div_test / 0.0f;
    assert(div_result == Vec2::zero());
    
    Vec2 div_assign_test(5.0f, 6.0f);
    div_assign_test /= 0.0f;
    assert(div_assign_test == Vec2::zero());
    
    // Bounds checking tests (valid indices)
    Vec2 bounds_test(10.0f, 20.0f);
    assert(bounds_test[0] == 10.0f);
    assert(bounds_test[1] == 20.0f);
    
    // Invalid indices should return 0.0f in debug mode, but we can't test assertions here
    // In release mode, behavior is undefined but won't crash with our implementation
    
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
    
    // Division by zero protection tests
    Vec3 div_test(1.0f, 2.0f, 3.0f);
    Vec3 div_result = div_test / 0.0f;
    assert(div_result == Vec3::zero());
    
    Vec3 div_assign_test(4.0f, 5.0f, 6.0f);
    div_assign_test /= 0.0f;
    assert(div_assign_test == Vec3::zero());
    
    // Bounds checking tests (valid indices)
    Vec3 bounds_test(10.0f, 20.0f, 30.0f);
    assert(bounds_test[0] == 10.0f);
    assert(bounds_test[1] == 20.0f);
    assert(bounds_test[2] == 30.0f);
    
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
    
    // Division by zero protection tests
    Vec4 div_test(1.0f, 2.0f, 3.0f, 4.0f);
    Vec4 div_result = div_test / 0.0f;
    assert(div_result == Vec4::zero());
    
    Vec4 div_assign_test(5.0f, 6.0f, 7.0f, 8.0f);
    div_assign_test /= 0.0f;
    assert(div_assign_test == Vec4::zero());
    
    // Bounds checking tests (valid indices)
    Vec4 bounds_test(10.0f, 20.0f, 30.0f, 40.0f);
    assert(bounds_test[0] == 10.0f);
    assert(bounds_test[1] == 20.0f);
    assert(bounds_test[2] == 30.0f);
    assert(bounds_test[3] == 40.0f);
    
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
    
    // Transpose - using fromRows for intuitive row-major input
    Mat4 test_matrix = Mat4::fromRows(1, 2, 3, 4,
                                      5, 6, 7, 8,
                                      9, 10, 11, 12,
                                      13, 14, 15, 16);
    Mat4 transposed = test_matrix.transposed();
    assert(transposed(0, 1) == 5.0f && transposed(1, 0) == 2.0f);
    
    // Test matrix inverse with known matrices
    // Test 1: Identity matrix inverse should be identity
    Mat4 identity_inv = identity.inverse();
    assert((identity_inv == identity));
    
    // Test 2: Translation matrix inverse
    Mat4 trans = Mat4::translation(Vec3(5.0f, 10.0f, 15.0f));
    Mat4 trans_inv = trans.inverse();
    Mat4 trans_product = trans * trans_inv;
    // Product should be identity (within epsilon)
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            float expected = (i == j) ? 1.0f : 0.0f;
            assert(std::abs(trans_product(i, j) - expected) < EPSILON * 10.0f);
        }
    }
    
    // Test 3: Scale matrix inverse
    Mat4 scale_mat = Mat4::scale(Vec3(2.0f, 3.0f, 4.0f));
    Mat4 scale_inv = scale_mat.inverse();
    Mat4 scale_product = scale_mat * scale_inv;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            float expected = (i == j) ? 1.0f : 0.0f;
            assert(std::abs(scale_product(i, j) - expected) < EPSILON * 10.0f);
        }
    }
    
    // Test 4: Rotation matrix inverse (should equal transpose for rotation matrices)
    Mat4 rot = Mat4::rotationZ(math::radians(45.0f));
    Mat4 rot_inv = rot.inverse();
    Mat4 rot_transpose = rot.transposed();
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            assert(std::abs(rot_inv(i, j) - rot_transpose(i, j)) < EPSILON * 10.0f);
        }
    }
    
    // Test 5: General invertible matrix - using fromRows for clarity
    Mat4 general = Mat4::fromRows(2, 0, 0, 1,
                                  0, 3, 0, 2,
                                  0, 0, 4, 3,
                                  0, 0, 0, 1);
    Mat4 general_inv = general.inverse();
    Mat4 general_product = general * general_inv;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            float expected = (i == j) ? 1.0f : 0.0f;
            assert(std::abs(general_product(i, j) - expected) < EPSILON * 10.0f);
        }
    }
    
    // Test 6: Constructor parameter order verification
    // Test column-major constructor
    Mat4 col_major(1, 5, 9, 13,   // Column 0: (1,5,9,13)
                   2, 6, 10, 14,  // Column 1: (2,6,10,14)
                   3, 7, 11, 15,  // Column 2: (3,7,11,15)
                   4, 8, 12, 16); // Column 3: (4,8,12,16)
    
    // Verify the matrix looks like this when accessed by (row, col):
    // [1  2  3  4 ]
    // [5  6  7  8 ]
    // [9  10 11 12]
    // [13 14 15 16]
    assert(col_major(0, 0) == 1.0f && col_major(0, 1) == 2.0f && col_major(0, 2) == 3.0f && col_major(0, 3) == 4.0f);
    assert(col_major(1, 0) == 5.0f && col_major(1, 1) == 6.0f && col_major(1, 2) == 7.0f && col_major(1, 3) == 8.0f);
    assert(col_major(2, 0) == 9.0f && col_major(2, 1) == 10.0f && col_major(2, 2) == 11.0f && col_major(2, 3) == 12.0f);
    assert(col_major(3, 0) == 13.0f && col_major(3, 1) == 14.0f && col_major(3, 2) == 15.0f && col_major(3, 3) == 16.0f);
    
    // Test fromRows constructor (should create the same matrix)
    Mat4 row_major = Mat4::fromRows(1, 2, 3, 4,
                                    5, 6, 7, 8,
                                    9, 10, 11, 12,
                                    13, 14, 15, 16);
    
    // Both matrices should be identical
    assert(col_major == row_major);
    
    // Test 7: Verify storage order matches OpenGL/GLM convention
    // Create a simple translation matrix using column-major constructor
    Mat4 translation_col(1, 0, 0, 0,  // Column 0: (1,0,0,0)
                         0, 1, 0, 0,  // Column 1: (0,1,0,0)
                         0, 0, 1, 0,  // Column 2: (0,0,1,0)
                         5, 3, 2, 1); // Column 3: (5,3,2,1) - translation
    
    // This should be equivalent to Mat4::translation(Vec3(5,3,2))
    Mat4 translation_factory = Mat4::translation(Vec3(5.0f, 3.0f, 2.0f));
    assert(translation_col == translation_factory);
    
    // Test point transformation
    Vec3 test_point(1.0f, 1.0f, 1.0f);
    Vec3 transformed1 = translation_col.transformPoint(test_point);
    Vec3 transformed2 = translation_factory.transformPoint(test_point);
    assert(std::abs(transformed1.x - 6.0f) < EPSILON && 
           std::abs(transformed1.y - 4.0f) < EPSILON && 
           std::abs(transformed1.z - 3.0f) < EPSILON);
    assert(transformed1.x == transformed2.x && transformed1.y == transformed2.y && transformed1.z == transformed2.z);
    
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
    
    // Fast inverse square root tests
    float test_values[] = {1.0f, 4.0f, 9.0f, 16.0f, 25.0f, 0.25f};
    for (float x : test_values) {
        float fast_result = math::fastInverseSqrt(x);
        float standard_result = math::standardInverseSqrt(x);
        
        // With 2 Newton-Raphson iterations, error should be very small
        float error = std::abs(fast_result - standard_result);
        float relative_error = (standard_result != 0.0f) ? error / standard_result : 0.0f;
        
        // Allow for small numerical differences
        assert(relative_error < 0.001f);  // Less than 0.1% error
    }
    
    // Test edge cases
    assert(math::fastInverseSqrt(0.0f) == 0.0f);
    assert(math::fastInverseSqrt(-1.0f) == 0.0f);
    
    // Smoothstep function tests
    // Normal case: edge0 < edge1
    assert(math::smoothstep(0.0f, 1.0f, -0.5f) == 0.0f);  // Below range
    assert(math::smoothstep(0.0f, 1.0f, 0.0f) == 0.0f);   // At edge0
    assert(math::smoothstep(0.0f, 1.0f, 0.5f) == 0.5f);   // Middle (should be 0.5 for smoothstep)
    assert(math::smoothstep(0.0f, 1.0f, 1.0f) == 1.0f);   // At edge1
    assert(math::smoothstep(0.0f, 1.0f, 1.5f) == 1.0f);   // Above range
    
    // Edge case: edge0 == edge1 (zero-width transition)
    assert(math::smoothstep(0.5f, 0.5f, 0.0f) == 0.0f);   // Below step
    assert(math::smoothstep(0.5f, 0.5f, 0.5f) == 1.0f);   // At step
    assert(math::smoothstep(0.5f, 0.5f, 1.0f) == 1.0f);   // Above step
    
    // Edge case: edge1 < edge0 (inverted range)
    assert(math::smoothstep(1.0f, 0.0f, -0.5f) == 1.0f);  // Below inverted range
    assert(math::smoothstep(1.0f, 0.0f, 0.0f) == 1.0f);   // At edge1 (inverted)
    assert(math::smoothstep(1.0f, 0.0f, 0.5f) == 0.5f);   // Middle (inverted)
    assert(math::smoothstep(1.0f, 0.0f, 1.0f) == 0.0f);   // At edge0 (inverted)
    assert(math::smoothstep(1.0f, 0.0f, 1.5f) == 0.0f);   // Above inverted range
    
    // Test that no division by zero occurs
    float result1 = math::smoothstep(2.0f, 2.0f, 1.0f);   // Should not crash
    float result2 = math::smoothstep(3.0f, 3.0f, 5.0f);   // Should not crash
    assert(result1 == 0.0f);
    assert(result2 == 1.0f);
    
    std::cout << "Math Utils tests passed!" << std::endl;
}

int main() {
    std::cout << "Running Velocity Math Library Tests" << std::endl;
    std::cout << "Version: " << getVersionString() << std::endl;
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