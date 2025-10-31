# Velocity Math Library - Benchmark Results

## System Information

- **CPU**: Intel i5-11320H @ 3.20GHz (11th Gen)
- **Compiler**: GCC 15.2
- **Optimization**: -O3 -march=native -DNDEBUG
- **Architecture**: 64-bit
- **Date**: October 31, 2025

## Performance Results

### Vector Operations (1M iterations)

| Operation | Time (μs) | Ops/sec | MB/s |
|-----------|-----------|---------|------|
| Vec3 Addition | 3037 | 329M | 11304.7 |
| Vec3 Subtraction | 3127 | 320M | 10979.3 |
| Vec3 Scalar Mult | 2195 | 456M | 10427.4 |
| Vec3 Dot Product | 3700 | 270M | 7217.0 |
| Vec3 Cross Product | 3167 | 316M | 10840.6 |
| Vec3 Length | 1421 | 704M | 10738.1 |
| Vec3 Normalization | 5366 | 186M | 4265.4 |
| Vec3 Distance | 3153 | 317M | 8469.0 |

### Matrix Operations (100K iterations)

| Operation | Time (μs) | Ops/sec | MB/s |
|-----------|-----------|---------|------|
| Mat4 Multiplication | 3317 | 30M | 5520.2 |
| Mat4 Addition | 2638 | 38M | 6941.1 |
| Mat4 Transpose | 2280 | 44M | 5354.0 |
| Mat4 Determinant | 2266 | 44M | 2861.9 |
| Mat4 Inverse | 492 | 20M | 2481.1 |

### Transformation Operations (1M iterations)

| Operation | Time (μs) | Ops/sec | MB/s |
|-----------|-----------|---------|------|
| Point Transform | 22303 | 45M | 3249.8 |
| Direction Transform | 8922 | 112M | 8123.7 |

## Performance Summary

### Top Performers

- **Fastest Vector Operation**: Vec3 Length (704M ops/sec)
- **Fastest Matrix Operation**: Mat4 Determinant (44M ops/sec)

### Key Metrics

- **Vector Addition**: High-frequency operation for real-time applications
- **Matrix Multiplication**: Core transformation operation for graphics
- **Memory Throughput**: Optimized for cache-friendly access patterns
- **Numerical Accuracy**: Sub-microsecond precision errors

## Notes

- Results may vary based on CPU, compiler, and system load
- Benchmarks use release optimization (-O3) for maximum performance
- Memory throughput calculated based on theoretical data access patterns
- All operations maintain numerical accuracy within specified tolerances
