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
| Vec3 Addition | 5071 | 197M | 6770.3 |
| Vec3 Subtraction | 5538 | 181M | 6199.4 |
| Vec3 Scalar Mult | 3998 | 250M | 5724.9 |
| Vec3 Dot Product | 5001 | 200M | 5339.5 |
| Vec3 Cross Product | 5095 | 196M | 6738.4 |
| Vec3 Length | 2452 | 408M | 6223.0 |
| Vec3 Normalization | 7049 | 142M | 3247.0 |
| Vec3 Distance | 2888 | 346M | 9246.1 |

### Matrix Operations (100K iterations)

| Operation | Time (μs) | Ops/sec | MB/s |
|-----------|-----------|---------|------|
| Mat4 Multiplication | 3929 | 25M | 4660.4 |
| Mat4 Addition | 2392 | 42M | 7654.9 |
| Mat4 Transpose | 2893 | 35M | 4219.5 |
| Mat4 Determinant | 2086 | 48M | 3108.8 |
| Mat4 Inverse | 3538 | 3M | 345.0 |

### Transformation Operations (1M iterations)

| Operation | Time (μs) | Ops/sec | MB/s |
|-----------|-----------|---------|------|
| Point Transform | 27791 | 36M | 2608.0 |
| Direction Transform | 21699 | 46M | 3340.2 |

## Performance Summary

### Top Performers

- **Fastest Vector Operation**: Vec3 Length (408M ops/sec)
- **Fastest Matrix Operation**: Mat4 Determinant (48M ops/sec)

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
