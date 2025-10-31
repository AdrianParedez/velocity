# Documentation Overview

This directory contains comprehensive documentation for the Velocity Math Library, including API references, performance benchmarks, and usage examples.

### Available Documentation

- **[API Reference](API_REFERENCE.md)** - Complete API documentation for all classes and functions
- **[Benchmark Results](BENCHMARK_RESULTS.md)** - Performance analysis and benchmark data

## Library Features

- **High Performance**: Optimized vector and matrix operations with up to 704M+ operations per second
- **Modern C++**: C++17 with constexpr operations and zero-cost abstractions
- **Graphics Ready**: Column-major matrices, projection utilities, and transformation functions
- **Cross-Platform**: Windows, Linux, and macOS support with multiple build systems
- **Comprehensive**: Vec2, Vec3, Vec4 vectors and Mat4 matrices with full mathematical operations
- **Production Ready**: Robust error handling, bounds checking, and numerical stability
- **Flexible**: Compile-time configuration options for performance vs safety trade-offs

## Quick Reference

### Vector Operations
```cpp
#include <velocity/velocity.h>
using namespace vel;

// Vector creation and basic operations
Vec3 position(1.0f, 2.0f, 3.0f);
Vec3 velocity(0.5f, 0.0f, -1.0f);
Vec3 newPosition = position + velocity * deltaTime;

// Vector utilities
float length = position.length();
Vec3 normalized = position.normalized();
float dotProduct = position.dot(velocity);
Vec3 crossProduct = position.cross(velocity);
```

### Matrix Operations
```cpp
// Matrix transformations
Mat4 transform = Mat4::translation(position) * 
                 Mat4::rotationY(math::radians(45.0f)) * 
                 Mat4::scale(Vec3(2.0f));

// Transform points and vectors
Vec3 transformedPoint = transform.transformPoint(Vec3(1, 0, 0));
Vec3 transformedVector = transform.transformVector(Vec3(0, 1, 0));

// Projection matrices
Mat4 projection = Mat4::perspective(math::radians(60.0f), 16.0f/9.0f, 0.1f, 100.0f);
Mat4 ortho = Mat4::orthographic(-10.0f, 10.0f, -10.0f, 10.0f, 0.1f, 100.0f);
```

## Performance Highlights

Latest benchmark results on Intel i5-11320H @ 3.20GHz with GCC 15.2 (-O3):

### Top Performers
- **Vec3 Length**: 704M ops/sec (10.7 GB/s) - Fastest vector operation
- **Vec3 Scalar Multiplication**: 456M ops/sec (10.4 GB/s)
- **Mat4 Determinant**: 44M ops/sec (2.9 GB/s) - Fastest matrix operation
- **Mat4 Transpose**: 44M ops/sec (5.4 GB/s)

### Real-World Performance
- **Vector Addition**: 329M ops/sec - Critical for physics simulations
- **Matrix Multiplication**: 30M ops/sec - Core graphics transformation
- **Point Transformation**: 45M ops/sec - 3D rendering pipeline
- **Cross Product**: 316M ops/sec - Lighting and collision detection

### Version 1.1.0 Improvements
- ✅ **Enhanced Safety**: Zero division and bounds checking protection
- ✅ **Better Accuracy**: Improved matrix inverse with proper cofactor calculation
- ✅ **Optimized Math**: Dual Newton-Raphson iteration for fast inverse sqrt
- ✅ **Robust Edge Cases**: Comprehensive error handling and validation

For complete performance analysis and detailed benchmarks, see [BENCHMARK_RESULTS.md](BENCHMARK_RESULTS.md).

## Building and Integration

### Using Make
```bash
make all          # Build library, examples, tests, and benchmarks
make run-tests    # Build and run tests
make run-benchmarks # Build and run performance benchmarks
```

### Using CMake
```bash
mkdir build && cd build
cmake ..
make
```

### Integration
```cpp
// Include the main header
#include <velocity/velocity.h>

// Use the velocity namespace
using namespace vel;

// All vector and matrix types are now available
Vec2 screenPos(100.0f, 200.0f);
Vec3 worldPos(1.0f, 2.0f, 3.0f);
Vec4 color(1.0f, 0.5f, 0.2f, 1.0f);
Mat4 transform = Mat4::identity();
```