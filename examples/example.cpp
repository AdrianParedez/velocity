#include <velocity/velocity.h>
#include <iostream>

using namespace vel;

void demonstrateVec2() {
    std::cout << "=== Vec2 Demo ===" << std::endl;
    
    Vec2 a(3.0f, 4.0f);
    Vec2 b(1.0f, 2.0f);
    
    std::cout << "a = " << a.toString() << std::endl;
    std::cout << "b = " << b.toString() << std::endl;
    std::cout << "a + b = " << (a + b).toString() << std::endl;
    std::cout << "a - b = " << (a - b).toString() << std::endl;
    std::cout << "a * 2 = " << (a * 2.0f).toString() << std::endl;
    std::cout << "a.dot(b) = " << a.dot(b) << std::endl;
    std::cout << "a.cross(b) = " << a.cross(b) << std::endl;
    std::cout << "a.length() = " << a.length() << std::endl;
    std::cout << "a.normalized() = " << a.normalized().toString() << std::endl;
    std::cout << "a.distance(b) = " << a.distance(b) << std::endl;
    std::cout << "a.angle() = " << math::degrees(a.angle()) << " degrees" << std::endl;
    std::cout << std::endl;
}

void demonstrateVec3() {
    std::cout << "=== Vec3 Demo ===" << std::endl;
    
    Vec3 a(1.0f, 2.0f, 3.0f);
    Vec3 b(4.0f, 5.0f, 6.0f);
    
    std::cout << "a = " << a.toString() << std::endl;
    std::cout << "b = " << b.toString() << std::endl;
    std::cout << "a + b = " << (a + b).toString() << std::endl;
    std::cout << "a - b = " << (a - b).toString() << std::endl;
    std::cout << "a * 2 = " << (a * 2.0f).toString() << std::endl;
    std::cout << "a.dot(b) = " << a.dot(b) << std::endl;
    std::cout << "a.cross(b) = " << a.cross(b).toString() << std::endl;
    std::cout << "a.length() = " << a.length() << std::endl;
    std::cout << "a.normalized() = " << a.normalized().toString() << std::endl;
    std::cout << "a.distance(b) = " << a.distance(b) << std::endl;
    
    // Demonstrate static constants
    std::cout << "Vec3::up() = " << Vec3::up().toString() << std::endl;
    std::cout << "Vec3::forward() = " << Vec3::forward().toString() << std::endl;
    std::cout << std::endl;
}

void demonstrateMat4() {
    std::cout << "=== Mat4 Demo ===" << std::endl;
    
    // Identity matrix
    Mat4 identity = Mat4::identity();
    std::cout << "Identity matrix:" << std::endl;
    std::cout << identity.toString(1) << std::endl;
    
    // Translation matrix
    Mat4 translation = Mat4::translation(Vec3(1.0f, 2.0f, 3.0f));
    std::cout << "Translation matrix:" << std::endl;
    std::cout << translation.toString(1) << std::endl;
    
    // Scale matrix
    Mat4 scale = Mat4::scale(Vec3(2.0f, 3.0f, 4.0f));
    std::cout << "Scale matrix:" << std::endl;
    std::cout << scale.toString(1) << std::endl;
    
    // Rotation matrix (45 degrees around Y axis)
    Mat4 rotation = Mat4::rotationY(math::radians(45.0f));
    std::cout << "Rotation matrix (45° around Y):" << std::endl;
    std::cout << rotation.toString(3) << std::endl;
    
    // Combined transformation
    Mat4 transform = translation * rotation * scale;
    std::cout << "Combined transformation (T * R * S):" << std::endl;
    std::cout << transform.toString(3) << std::endl;
    
    // Transform a point
    Vec3 point(1.0f, 0.0f, 0.0f);
    Vec3 transformed = transform.transformPoint(point);
    std::cout << "Point " << point.toString() << " transformed to " << transformed.toString() << std::endl;
    std::cout << std::endl;
}

void demonstrateProjectionMatrices() {
    std::cout << "=== Projection Matrices ===" << std::endl;
    
    // Perspective projection
    Mat4 perspective = Mat4::perspective(math::radians(60.0f), 16.0f/9.0f, 0.1f, 100.0f);
    std::cout << "Perspective projection (60° FOV, 16:9 aspect):" << std::endl;
    std::cout << perspective.toString(3) << std::endl;
    
    // Orthographic projection
    Mat4 ortho = Mat4::ortho(-10.0f, 10.0f, -10.0f, 10.0f, 0.1f, 100.0f);
    std::cout << "Orthographic projection:" << std::endl;
    std::cout << ortho.toString(3) << std::endl;
    
    // Look-at matrix
    Mat4 lookAt = Mat4::lookAt(Vec3(0.0f, 0.0f, 5.0f), Vec3::zero(), Vec3::up());
    std::cout << "Look-at matrix (eye at (0,0,5), looking at origin):" << std::endl;
    std::cout << lookAt.toString(3) << std::endl;
    std::cout << std::endl;
}

void demonstrateMathUtils() {
    std::cout << "=== Math Utils Demo ===" << std::endl;
    
    std::cout << "PI = " << PI << std::endl;
    std::cout << "45 degrees in radians = " << math::radians(45.0f) << std::endl;
    std::cout << "PI/4 radians in degrees = " << math::degrees(PI/4.0f) << std::endl;
    std::cout << "clamp(1.5, 0, 1) = " << math::clamp(1.5f, 0.0f, 1.0f) << std::endl;
    std::cout << "lerp(0, 10, 0.5) = " << math::lerp(0.0f, 10.0f, 0.5f) << std::endl;
    std::cout << "smoothstep(0, 1, 0.5) = " << math::smoothstep(0.0f, 1.0f, 0.5f) << std::endl;
    std::cout << "sign(-3.14) = " << math::sign(-3.14f) << std::endl;
    std::cout << "isPowerOfTwo(16) = " << (math::isPowerOfTwo(16) ? "true" : "false") << std::endl;
    std::cout << std::endl;
}

int main() {
    std::cout << "Velocity Math Library Demo" << std::endl;
    std::cout << "Version " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << std::endl;
    std::cout << "==================" << std::endl << std::endl;
    
    demonstrateVec2();
    demonstrateVec3();
    demonstrateMat4();
    demonstrateProjectionMatrices();
    demonstrateMathUtils();
    
    return 0;
}