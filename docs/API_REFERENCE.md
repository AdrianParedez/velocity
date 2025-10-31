# Velocity Math Library - API Reference

## Table of Contents

1. [Namespace](#namespace)
2. [Vec2 - 2D Vector](#vec2---2d-vector)
3. [Vec3 - 3D Vector](#vec3---3d-vector)
4. [Vec4 - 4D Vector](#vec4---4d-vector)
5. [Mat4 - 4x4 Matrix](#mat4---4x4-matrix)
6. [Math Utilities](#math-utilities)
7. [Constants](#constants)
8. [Examples](#examples)

## Namespace

All Velocity math library components are contained within the `vel` namespace:

```cpp
#include <velocity/velocity.h>
using namespace vel;  // Optional: brings all types into scope
```

## Vec2 - 2D Vector

### Declaration
```cpp
class Vec2 {
public:
    float x, y;
};
```

### Constructors
```cpp
Vec2()                          // Zero vector (0, 0)
Vec2(float x, float y)          // Component constructor
Vec2(float scalar)              // Uniform constructor (scalar, scalar)
```

### Arithmetic Operations
```cpp
Vec2 operator+(const Vec2& rhs) const    // Vector addition
Vec2 operator-(const Vec2& rhs) const    // Vector subtraction
Vec2 operator*(float scalar) const       // Scalar multiplication
Vec2 operator/(float scalar) const       // Scalar division
Vec2 operator-() const                   // Negation

Vec2& operator+=(const Vec2& rhs)        // Compound addition
Vec2& operator-=(const Vec2& rhs)        // Compound subtraction
Vec2& operator*=(float scalar)           // Compound scalar multiplication
Vec2& operator/=(float scalar)           // Compound scalar division
```

### Comparison Operations
```cpp
bool operator==(const Vec2& rhs) const   // Equality (with epsilon tolerance)
bool operator!=(const Vec2& rhs) const   // Inequality
```

### Vector Operations
```cpp
float dot(const Vec2& rhs) const         // Dot product
float cross(const Vec2& rhs) const       // 2D cross product (scalar result)
float length() const                     // Vector magnitude
float lengthSquared() const              // Squared magnitude (faster)
Vec2 normalized() const                  // Unit vector
Vec2& normalize()                        // Normalize in-place
float distance(const Vec2& rhs) const    // Distance to another vector
float distanceSquared(const Vec2& rhs) const // Squared distance (faster)
Vec2 lerp(const Vec2& rhs, float t) const    // Linear interpolation
Vec2 reflect(const Vec2& normal) const   // Reflection across normal
float angle() const                      // Angle in radians
Vec2 rotate(float radians) const         // Rotation by angle
Vec2 perpendicular() const               // Perpendicular vector (-y, x)
```

### Element Access
```cpp
float& operator[](int index)             // Mutable element access
const float& operator[](int index) const // Const element access
```

### Static Constants
```cpp
static Vec2 zero()                       // (0, 0)
static Vec2 one()                        // (1, 1)
static Vec2 unitX()                      // (1, 0)
static Vec2 unitY()                      // (0, 1)
```

### Utility
```cpp
std::string toString(int precision = 3) const  // String representation
```

## Vec3 - 3D Vector

### Declaration
```cpp
class Vec3 {
public:
    float x, y, z;
};
```

### Constructors
```cpp
Vec3()                                   // Zero vector (0, 0, 0)
Vec3(float x, float y, float z)          // Component constructor
Vec3(float scalar)                       // Uniform constructor
Vec3(const Vec2& xy, float z)            // From Vec2 + z component
```

### Arithmetic Operations
```cpp
Vec3 operator+(const Vec3& rhs) const    // Vector addition
Vec3 operator-(const Vec3& rhs) const    // Vector subtraction
Vec3 operator*(float scalar) const       // Scalar multiplication
Vec3 operator/(float scalar) const       // Scalar division
Vec3 operator-() const                   // Negation

Vec3& operator+=(const Vec3& rhs)        // Compound operations
Vec3& operator-=(const Vec3& rhs)
Vec3& operator*=(float scalar)
Vec3& operator/=(float scalar)
```

### Comparison Operations
```cpp
bool operator==(const Vec3& rhs) const   // Equality (with epsilon tolerance)
bool operator!=(const Vec3& rhs) const   // Inequality
```

### Vector Operations
```cpp
float dot(const Vec3& rhs) const         // Dot product
Vec3 cross(const Vec3& rhs) const        // Cross product
float length() const                     // Vector magnitude
float lengthSquared() const              // Squared magnitude
Vec3 normalized() const                  // Unit vector
Vec3& normalize()                        // Normalize in-place
float distance(const Vec3& rhs) const    // Distance to another vector
float distanceSquared(const Vec3& rhs) const // Squared distance
Vec3 lerp(const Vec3& rhs, float t) const    // Linear interpolation
Vec3 reflect(const Vec3& normal) const   // Reflection across normal
Vec3 project(const Vec3& onto) const     // Vector projection
Vec3 reject(const Vec3& onto) const      // Vector rejection
```

### Swizzling
```cpp
Vec2 xy() const                          // Extract (x, y)
Vec2 xz() const                          // Extract (x, z)
Vec2 yz() const                          // Extract (y, z)
```

### Element Access
```cpp
float& operator[](int index)             // Mutable element access
const float& operator[](int index) const // Const element access
```

### Static Constants
```cpp
static Vec3 zero()                       // (0, 0, 0)
static Vec3 one()                        // (1, 1, 1)
static Vec3 unitX()                      // (1, 0, 0)
static Vec3 unitY()                      // (0, 1, 0)
static Vec3 unitZ()                      // (0, 0, 1)
static Vec3 forward()                    // (0, 0, -1) - Forward direction
static Vec3 back()                       // (0, 0, 1)
static Vec3 up()                         // (0, 1, 0)
static Vec3 down()                       // (0, -1, 0)
static Vec3 right()                      // (1, 0, 0)
static Vec3 left()                       // (-1, 0, 0)
```

### Utility
```cpp
std::string toString(int precision = 3) const  // String representation
```

## Vec4 - 4D Vector

### Declaration
```cpp
class Vec4 {
public:
    float x, y, z, w;
};
```

### Constructors
```cpp
Vec4()                                   // Zero vector (0, 0, 0, 0)
Vec4(float x, float y, float z, float w) // Component constructor
Vec4(float scalar)                       // Uniform constructor
Vec4(const Vec3& xyz, float w)           // From Vec3 + w component
Vec4(const Vec2& xy, float z, float w)   // From Vec2 + z, w components
```

### Arithmetic Operations
```cpp
Vec4 operator+(const Vec4& rhs) const    // Vector addition
Vec4 operator-(const Vec4& rhs) const    // Vector subtraction
Vec4 operator*(float scalar) const       // Scalar multiplication
Vec4 operator/(float scalar) const       // Scalar division
Vec4 operator-() const                   // Negation

Vec4& operator+=(const Vec4& rhs)        // Compound operations
Vec4& operator-=(const Vec4& rhs)
Vec4& operator*=(float scalar)
Vec4& operator/=(float scalar)
```

### Comparison Operations
```cpp
bool operator==(const Vec4& rhs) const   // Equality (with epsilon tolerance)
bool operator!=(const Vec4& rhs) const   // Inequality
```

### Vector Operations
```cpp
float dot(const Vec4& rhs) const         // Dot product
float length() const                     // Vector magnitude
float lengthSquared() const              // Squared magnitude
Vec4 normalized() const                  // Unit vector
Vec4& normalize()                        // Normalize in-place
Vec4 lerp(const Vec4& rhs, float t) const    // Linear interpolation
```

### Swizzling
```cpp
Vec2 xy() const                          // Extract (x, y)
Vec2 zw() const                          // Extract (z, w)
Vec3 xyz() const                         // Extract (x, y, z)
```

### Element Access
```cpp
float& operator[](int index)             // Mutable element access
const float& operator[](int index) const // Const element access
```

### Static Constants
```cpp
static Vec4 zero()                       // (0, 0, 0, 0)
static Vec4 one()                        // (1, 1, 1, 1)
static Vec4 unitX()                      // (1, 0, 0, 0)
static Vec4 unitY()                      // (0, 1, 0, 0)
static Vec4 unitZ()                      // (0, 0, 1, 0)
static Vec4 unitW()                      // (0, 0, 0, 1)
```

### Utility
```cpp
std::string toString(int precision = 3) const  // String representation
```

## Mat4 - 4x4 Matrix

### Declaration
```cpp
class Mat4 {
    // Column-major storage: m[col * 4 + row]
};
```

### Constructors
```cpp
Mat4()                                   // Zero matrix
Mat4(float diagonal)                     // Diagonal matrix
Mat4(float m00, float m01, ..., float m33) // Component constructor
Mat4(const std::array<float, 16>& values)  // Array constructor
```

### Element Access
```cpp
float& operator()(int row, int col)      // Element access (row, col)
const float& operator()(int row, int col) const
float* data()                            // Raw data pointer
const float* data() const
```

### Arithmetic Operations
```cpp
Mat4 operator+(const Mat4& rhs) const    // Matrix addition
Mat4 operator-(const Mat4& rhs) const    // Matrix subtraction
Mat4 operator*(const Mat4& rhs) const    // Matrix multiplication
Mat4 operator*(float scalar) const       // Scalar multiplication
Vec4 operator*(const Vec4& vec) const    // Matrix-vector multiplication

Mat4& operator+=(const Mat4& rhs)        // Compound operations
Mat4& operator-=(const Mat4& rhs)
Mat4& operator*=(const Mat4& rhs)
Mat4& operator*=(float scalar)
```

### Comparison Operations
```cpp
bool operator==(const Mat4& rhs) const   // Equality (with epsilon tolerance)
bool operator!=(const Mat4& rhs) const   // Inequality
```

### Matrix Operations
```cpp
Mat4 transposed() const                  // Transpose matrix
Mat4& transpose()                        // Transpose in-place
float determinant() const                // Matrix determinant
Mat4 inverse() const                     // Inverse matrix
Mat4& invert()                           // Invert in-place
```

### Transform Operations
```cpp
Vec3 transformPoint(const Vec3& point) const      // Transform point (w=1)
Vec3 transformDirection(const Vec3& direction) const // Transform direction (w=0)
Vec3 transformNormal(const Vec3& normal) const    // Transform normal (inverse transpose)
```

### Static Factory Methods
```cpp
static Mat4 identity()                   // Identity matrix
static Mat4 translation(const Vec3& translation) // Translation matrix
static Mat4 scale(const Vec3& scale)     // Scale matrix
static Mat4 scale(float uniform_scale)   // Uniform scale matrix
static Mat4 rotationX(float radians)     // Rotation around X-axis
static Mat4 rotationY(float radians)     // Rotation around Y-axis
static Mat4 rotationZ(float radians)     // Rotation around Z-axis
static Mat4 rotation(const Vec3& axis, float radians) // Rotation around arbitrary axis
static Mat4 lookAt(const Vec3& eye, const Vec3& center, const Vec3& up) // View matrix
static Mat4 perspective(float fovy, float aspect, float near, float far) // Perspective projection
static Mat4 ortho(float left, float right, float bottom, float top, float near, float far) // Orthographic projection
static Mat4 frustum(float left, float right, float bottom, float top, float near, float far) // Frustum projection
```

### Decomposition
```cpp
Vec3 getTranslation() const              // Extract translation
Vec3 getScale() const                    // Extract scale
Mat4 getRotation() const                 // Extract rotation matrix
```

### Utility
```cpp
std::string toString(int precision = 3) const  // String representation
```

## Math Utilities

### Namespace
```cpp
namespace vel::math {
    // Utility functions
}
```

### Angle Conversion
```cpp
float radians(float degrees)             // Convert degrees to radians
float degrees(float radians)             // Convert radians to degrees
```

### Comparison Functions
```cpp
bool isEqual(float a, float b, float epsilon = EPSILON) // Epsilon-based equality
bool isZero(float value, float epsilon = EPSILON)       // Check if value is zero
```

### Clamping and Interpolation
```cpp
float clamp(float value, float min, float max)          // Clamp value to range
float lerp(float a, float b, float t)                   // Linear interpolation
float smoothstep(float edge0, float edge1, float x)     // Smooth interpolation
```

### Utility Functions
```cpp
int sign(float value)                    // Sign function (-1, 0, 1)
float min(float a, float b)              // Minimum value
float max(float a, float b)              // Maximum value
bool isPowerOfTwo(int value)             // Check if value is power of 2
float fastInverseSqrt(float x)           // Fast inverse square root approximation
```

## Constants

### Mathematical Constants
```cpp
namespace vel {
    constexpr float PI = 3.14159265358979323846f;        // Pi
    constexpr float TWO_PI = 2.0f * PI;                  // 2 * Pi
    constexpr float HALF_PI = PI * 0.5f;                 // Pi / 2
    constexpr float EPSILON = 1e-6f;                     // Default epsilon
    constexpr float DEG_TO_RAD = PI / 180.0f;            // Degree to radian conversion
    constexpr float RAD_TO_DEG = 180.0f / PI;            // Radian to degree conversion
}
```

## Examples

### Basic Vector Operations
```cpp
#include <velocity/velocity.h>
using namespace vel;

// Create vectors
Vec3 a(1.0f, 2.0f, 3.0f);
Vec3 b(4.0f, 5.0f, 6.0f);

// Basic operations
Vec3 sum = a + b;                        // (5, 7, 9)
Vec3 scaled = a * 2.0f;                  // (2, 4, 6)
float dot_product = a.dot(b);            // 32
Vec3 cross_product = a.cross(b);         // (-3, 6, -3)

// Vector properties
float length = a.length();               // 3.74166
Vec3 unit = a.normalized();              // (0.267, 0.535, 0.802)
```

### Matrix Transformations
```cpp
// Create transformation matrices
Mat4 translation = Mat4::translation(Vec3(1, 2, 3));
Mat4 rotation = Mat4::rotationY(math::radians(45.0f));
Mat4 scale = Mat4::scale(Vec3(2, 2, 2));

// Combine transformations (order matters: T * R * S)
Mat4 transform = translation * rotation * scale;

// Transform points
Vec3 point(1, 0, 0);
Vec3 transformed = transform.transformPoint(point);
```

### Projection Matrices
```cpp
// Perspective projection
Mat4 perspective = Mat4::perspective(
    math::radians(60.0f),  // Field of view
    16.0f / 9.0f,          // Aspect ratio
    0.1f,                  // Near plane
    100.0f                 // Far plane
);

// Orthographic projection
Mat4 ortho = Mat4::ortho(-10, 10, -10, 10, 0.1f, 100.0f);

// View matrix
Mat4 view = Mat4::lookAt(
    Vec3(0, 0, 5),         // Eye position
    Vec3(0, 0, 0),         // Look at target
    Vec3(0, 1, 0)          // Up vector
);
```

### Epsilon-Based Comparisons
```cpp
// Floating-point comparison
float a = 0.1f + 0.2f;
float b = 0.3f;

bool equal = (a == b);                   // May be false due to precision
bool epsilon_equal = math::isEqual(a, b); // True with epsilon tolerance

// Vector comparison
Vec3 v1(1.0f, 2.0f, 3.0f);
Vec3 v2(1.0000001f, 2.0f, 3.0f);
bool vectors_equal = (v1 == v2);         // True with default epsilon
```

### Performance Tips
```cpp
// Use squared operations when possible (faster)
float dist_sq = a.distanceSquared(b);    // Faster than distance()
float len_sq = a.lengthSquared();        // Faster than length()

// In-place operations avoid temporaries
Vec3 velocity(1, 2, 3);
velocity.normalize();                    // Faster than velocity = velocity.normalized()
velocity *= deltaTime;                   // Faster than velocity = velocity * deltaTime
```