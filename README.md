# Velocity Math Library

A high-performance C++17 vector and matrix math library designed for real-time graphics, game development, and scientific computing.

## Features

- **High Performance**: Optimized vector and matrix operations with up to 400M+ operations per second
- **Modern C++**: C++17 with constexpr operations and zero-cost abstractions
- **Graphics Ready**: Column-major matrices, projection utilities, and transformation functions
- **Cross-Platform**: Windows, Linux, and macOS support with multiple build systems
- **Comprehensive**: Vec2, Vec3, Vec4 vectors and Mat4 matrices with full mathematical operations

## Quick Start

```cpp
#include <velocity/velocity.h>
using namespace vel;

// Vector operations
Vec3 position(1.0f, 2.0f, 3.0f);
Vec3 velocity(0.5f, 0.0f, -1.0f);
Vec3 newPosition = position + velocity * deltaTime;

// Matrix transformations
Mat4 transform = Mat4::translation(position) * 
                 Mat4::rotationY(math::radians(45.0f)) * 
                 Mat4::scale(Vec3(2.0f));

// Transform points
Vec3 transformedPoint = transform.transformPoint(Vec3(1, 0, 0));
```

## Performance

Benchmark results on Intel i5-11320H @ 3.20GHz with GCC 15.2 (-O3):

### Vector Operations (1M iterations)
- **Vec3 Length**: 704M ops/sec (10.7 GB/s)
- **Vec3 Scalar Multiplication**: 456M ops/sec (10.4 GB/s)
- **Vec3 Addition**: 329M ops/sec (11.3 GB/s)
- **Vec3 Cross Product**: 316M ops/sec (10.8 GB/s)

### Matrix Operations (100K iterations)
- **Mat4 Determinant**: 44M ops/sec (2.9 GB/s)
- **Mat4 Transpose**: 44M ops/sec (5.4 GB/s)
- **Mat4 Addition**: 38M ops/sec (6.9 GB/s)
- **Mat4 Multiplication**: 30M ops/sec (5.5 GB/s)

### Key Improvements in v1.1.0
- ✅ **Zero Division Protection**: Safe vector division operations
- ✅ **Bounds Checking**: Debug-mode array access validation
- ✅ **Improved Matrix Inverse**: Proper cofactor-based calculation
- ✅ **Enhanced Precision**: Dual Newton-Raphson fast inverse sqrt
- ✅ **Better Error Handling**: Robust edge case management

See [docs/BENCHMARK_RESULTS.md](docs/BENCHMARK_RESULTS.md) for complete performance analysis.

## Building

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

## Documentation

Complete documentation is available in the [docs/](docs/) directory:

- [Benchmark Results](docs/BENCHMARK_RESULTS.md) - Latest performance data
- [Documentation Index](docs/README.md) - Complete documentation overview

## Requirements

- **C++17** compatible compiler
- **CMake 3.12+** (for CMake builds)
- **Make** (for Makefile builds)

## License

MIT License - see [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome! Please ensure all tests pass and include benchmarks for performance-related changes.

## Author

Adrian Paredez - [GitHub](https://github.com/AdrianParedez)