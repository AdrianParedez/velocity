# Documentation Overview

This directory contains comprehensive documentation for the Velocity Math Library, including API references, performance benchmarks, and usage examples.

### Available Documentation

- **[API Reference](API_REFERENCE.md)** - Complete API documentation for all classes and functions
- **[Benchmark Results](BENCHMARK_RESULTS.md)** - Performance analysis and benchmark data

## Library Features

- **High Performance**: Optimized vector and matrix operations with up to 400M+ operations per second
- **Modern C++**: C++17 with constexpr operations and zero-cost abstractions
- **Graphics Ready**: Column-major matrices, projection utilities, and transformation functions
- **Cross-Platform**: Windows, Linux, and macOS support with multiple build systems
- **Comprehensive**: Vec2, Vec3, Vec4 vectors and Mat4 matrices with full mathematical operations

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

Benchmark results on Intel i5-11320H @ 3.20GHz:

- **Vector Operations**: Up to 408M operations/second
- **Matrix Operations**: Up to 48M operations/second
- **Memory Throughput**: Up to 9.2GB/second

For detailed performance analysis, see [BENCHMARK_RESULTS.md](BENCHMARK_RESULTS.md).

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