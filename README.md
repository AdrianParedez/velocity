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

Benchmark results on Intel i5-11320H @ 3.20GHz:

- **Vector Operations**: Up to 408M operations/second
- **Matrix Operations**: Up to 48M operations/second
- **Memory Throughput**: Up to 9.2GB/second

See [docs/BENCHMARK_RESULTS.md](docs/BENCHMARK_RESULTS.md) for detailed performance analysis.

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