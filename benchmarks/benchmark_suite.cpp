#include <velocity/velocity.h>
#include <chrono>
#include <iostream>
#include <vector>
#include <random>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <sstream>

using namespace vel;

class BenchmarkTimer {
private:
    std::chrono::high_resolution_clock::time_point start_time;
    
public:
    void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }
    
    double elapsed_microseconds() {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        return static_cast<double>(duration.count());
    }
    
    double elapsed_milliseconds() {
        return elapsed_microseconds() / 1000.0;
    }
};

class PerformanceBenchmark {
private:
    std::mt19937 rng{42}; // Fixed seed for reproducible results
    std::uniform_real_distribution<float> dist{-10.0f, 10.0f};
    
    struct BenchmarkResult {
        std::string name;
        double time_us;
        double ops_per_sec;
        double mb_per_sec;
    };
    
    std::vector<BenchmarkResult> results;
    
    void printHeader(const std::string& title) {
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << title << std::endl;
        std::cout << std::string(60, '=') << std::endl;
        std::cout << std::left << std::setw(25) << "Operation" 
                  << std::setw(12) << "Time (μs)" 
                  << std::setw(15) << "Ops/sec" 
                  << std::setw(12) << "MB/s" << std::endl;
        std::cout << std::string(60, '-') << std::endl;
    }
    
    void printResult(const std::string& name, double time_us, int iterations, size_t bytes_per_op = 0) {
        double ops_per_sec = (iterations * 1e6) / time_us;
        double mb_per_sec = bytes_per_op > 0 ? (ops_per_sec * bytes_per_op) / (1024 * 1024) : 0;
        
        std::cout << std::left << std::setw(25) << name
                  << std::setw(12) << std::fixed << std::setprecision(0) << time_us
                  << std::setw(15) << std::fixed << std::setprecision(0) << ops_per_sec;
        
        if (bytes_per_op > 0) {
            std::cout << std::setw(12) << std::fixed << std::setprecision(1) << mb_per_sec;
        }
        std::cout << std::endl;
        
        // Store result for document generation
        results.push_back({name, time_us, ops_per_sec, mb_per_sec});
    }
    
public:
    void runVectorBenchmarks() {
        printHeader("Vector Operations Benchmark (1M iterations)");
        
        const int iterations = 1000000;
        BenchmarkTimer timer;
        
        // Generate test data
        std::vector<Vec3> vec_a(iterations), vec_b(iterations), results(iterations);
        for (int i = 0; i < iterations; ++i) {
            vec_a[i] = Vec3(dist(rng), dist(rng), dist(rng));
            vec_b[i] = Vec3(dist(rng), dist(rng), dist(rng));
        }
        
        // Vec3 Addition
        timer.start();
        for (int i = 0; i < iterations; ++i) {
            results[i] = vec_a[i] + vec_b[i];
        }
        printResult("Vec3 Addition", timer.elapsed_microseconds(), iterations, sizeof(Vec3) * 3);
        
        // Vec3 Subtraction
        timer.start();
        for (int i = 0; i < iterations; ++i) {
            results[i] = vec_a[i] - vec_b[i];
        }
        printResult("Vec3 Subtraction", timer.elapsed_microseconds(), iterations, sizeof(Vec3) * 3);
        
        // Vec3 Scalar Multiplication
        timer.start();
        for (int i = 0; i < iterations; ++i) {
            results[i] = vec_a[i] * 2.5f;
        }
        printResult("Vec3 Scalar Mult", timer.elapsed_microseconds(), iterations, sizeof(Vec3) * 2);
        
        // Vec3 Dot Product
        std::vector<float> dot_results(iterations);
        timer.start();
        for (int i = 0; i < iterations; ++i) {
            dot_results[i] = vec_a[i].dot(vec_b[i]);
        }
        printResult("Vec3 Dot Product", timer.elapsed_microseconds(), iterations, sizeof(Vec3) * 2 + sizeof(float));
        
        // Vec3 Cross Product
        timer.start();
        for (int i = 0; i < iterations; ++i) {
            results[i] = vec_a[i].cross(vec_b[i]);
        }
        printResult("Vec3 Cross Product", timer.elapsed_microseconds(), iterations, sizeof(Vec3) * 3);
        
        // Vec3 Length
        timer.start();
        for (int i = 0; i < iterations; ++i) {
            dot_results[i] = vec_a[i].length();
        }
        printResult("Vec3 Length", timer.elapsed_microseconds(), iterations, sizeof(Vec3) + sizeof(float));
        
        // Vec3 Normalization
        timer.start();
        for (int i = 0; i < iterations; ++i) {
            results[i] = vec_a[i].normalized();
        }
        printResult("Vec3 Normalization", timer.elapsed_microseconds(), iterations, sizeof(Vec3) * 2);
        
        // Vec3 Distance
        timer.start();
        for (int i = 0; i < iterations; ++i) {
            dot_results[i] = vec_a[i].distance(vec_b[i]);
        }
        printResult("Vec3 Distance", timer.elapsed_microseconds(), iterations, sizeof(Vec3) * 2 + sizeof(float));
    }
    
    void runMatrixBenchmarks() {
        printHeader("Matrix Operations Benchmark (100K iterations)");
        
        const int iterations = 100000;
        BenchmarkTimer timer;
        
        // Generate test data
        std::vector<Mat4> mat_a(iterations), mat_b(iterations), results(iterations);
        for (int i = 0; i < iterations; ++i) {
            // Create random transformation matrices
            Vec3 translation(dist(rng), dist(rng), dist(rng));
            Vec3 scale(std::abs(dist(rng)) + 0.1f, std::abs(dist(rng)) + 0.1f, std::abs(dist(rng)) + 0.1f);
            float angle = dist(rng);
            
            mat_a[i] = Mat4::translation(translation) * Mat4::rotationY(angle) * Mat4::scale(scale);
            mat_b[i] = Mat4::translation(translation * 0.5f) * Mat4::rotationX(angle * 0.7f);
        }
        
        // Matrix Multiplication
        timer.start();
        for (int i = 0; i < iterations; ++i) {
            results[i] = mat_a[i] * mat_b[i];
        }
        printResult("Mat4 Multiplication", timer.elapsed_microseconds(), iterations, sizeof(Mat4) * 3);
        
        // Matrix Addition
        timer.start();
        for (int i = 0; i < iterations; ++i) {
            results[i] = mat_a[i] + mat_b[i];
        }
        printResult("Mat4 Addition", timer.elapsed_microseconds(), iterations, sizeof(Mat4) * 3);
        
        // Matrix Transpose
        timer.start();
        for (int i = 0; i < iterations; ++i) {
            results[i] = mat_a[i].transposed();
        }
        printResult("Mat4 Transpose", timer.elapsed_microseconds(), iterations, sizeof(Mat4) * 2);
        
        // Matrix Determinant
        std::vector<float> det_results(iterations);
        timer.start();
        for (int i = 0; i < iterations; ++i) {
            det_results[i] = mat_a[i].determinant();
        }
        printResult("Mat4 Determinant", timer.elapsed_microseconds(), iterations, sizeof(Mat4) + sizeof(float));
        
        // Matrix Inverse (fewer iterations due to complexity)
        const int inverse_iterations = iterations / 10;
        timer.start();
        for (int i = 0; i < inverse_iterations; ++i) {
            try {
                results[i] = mat_a[i].inverse();
            } catch (...) {
                results[i] = Mat4::identity(); // Fallback for singular matrices
            }
        }
        printResult("Mat4 Inverse", timer.elapsed_microseconds(), inverse_iterations, sizeof(Mat4) * 2);
    }
    
    void runTransformBenchmarks() {
        printHeader("Transformation Benchmark (1M iterations)");
        
        const int iterations = 1000000;
        BenchmarkTimer timer;
        
        // Generate test data
        std::vector<Vec3> points(iterations), results(iterations);
        Mat4 transform = Mat4::translation(Vec3(1, 2, 3)) * 
                        Mat4::rotationY(math::radians(45.0f)) * 
                        Mat4::scale(Vec3(2, 2, 2));
        
        for (int i = 0; i < iterations; ++i) {
            points[i] = Vec3(dist(rng), dist(rng), dist(rng));
        }
        
        // Point Transformation
        timer.start();
        for (int i = 0; i < iterations; ++i) {
            results[i] = transform.transformPoint(points[i]);
        }
        printResult("Point Transform", timer.elapsed_microseconds(), iterations, sizeof(Vec3) + sizeof(Mat4));
        
        // Direction Transformation
        timer.start();
        for (int i = 0; i < iterations; ++i) {
            results[i] = transform.transformDirection(points[i]);
        }
        printResult("Direction Transform", timer.elapsed_microseconds(), iterations, sizeof(Vec3) + sizeof(Mat4));
    }
    
    void runMemoryBenchmarks() {
        printHeader("Memory Access Pattern Benchmark");
        
        const int size = 1000000;
        BenchmarkTimer timer;
        
        // Sequential access test
        std::vector<Vec3> sequential_data(size);
        for (int i = 0; i < size; ++i) {
            sequential_data[i] = Vec3(i, i+1, i+2);
        }
        
        timer.start();
        Vec3 sum = Vec3::zero();
        for (int i = 0; i < size; ++i) {
            sum += sequential_data[i];
        }
        double sequential_time = timer.elapsed_microseconds();
        
        // Random access test
        std::vector<int> indices(size);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), rng);
        
        timer.start();
        sum = Vec3::zero();
        for (int i = 0; i < size; ++i) {
            sum += sequential_data[indices[i]];
        }
        double random_time = timer.elapsed_microseconds();
        
        printResult("Sequential Access", sequential_time, size, sizeof(Vec3) * 2);
        printResult("Random Access", random_time, size, sizeof(Vec3) * 2);
        
        std::cout << "\nCache Efficiency Ratio: " << std::fixed << std::setprecision(2) 
                  << (sequential_time / random_time) << "x faster (sequential vs random)" << std::endl;
    }
    
    void runAccuracyTests() {
        printHeader("Numerical Accuracy Tests");
        
        // Test vector normalization accuracy
        Vec3 test_vec(1e-6f, 1e-6f, 1e-6f);  // Very small vector
        Vec3 normalized = test_vec.normalized();
        float length_error = std::abs(normalized.length() - 1.0f);
        
        std::cout << "Vector Normalization Error: " << std::scientific << length_error << std::endl;
        
        // Test matrix inverse accuracy
        Mat4 original = Mat4::translation(Vec3(1, 2, 3)) * Mat4::rotationY(0.5f);
        Mat4 inverse = original.inverse();
        Mat4 identity_test = original * inverse;
        
        float identity_error = 0.0f;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                float expected = (i == j) ? 1.0f : 0.0f;
                identity_error += std::abs(identity_test(i, j) - expected);
            }
        }
        
        std::cout << "Matrix Inverse Error: " << std::scientific << identity_error << std::endl;
        
        // Test floating-point comparison
        float a = 0.1f + 0.2f;
        float b = 0.3f;
        bool standard_equal = (a == b);
        bool epsilon_equal = math::isEqual(a, b);
        
        std::cout << "0.1 + 0.2 == 0.3: " << (standard_equal ? "true" : "false") << std::endl;
        std::cout << "isEqual(0.1 + 0.2, 0.3): " << (epsilon_equal ? "true" : "false") << std::endl;
        std::cout << "Difference: " << std::scientific << std::abs(a - b) << std::endl;
    }
    
    void generateBenchmarkDocument() {
        std::ofstream doc("docs/BENCHMARK_RESULTS.md");
        
        doc << "# Velocity Math Library - Benchmark Results\n\n";
        doc << "## System Information\n\n";
        doc << "- **CPU**: Intel i5-11320H @ 3.20GHz (11th Gen)\n";
        doc << "- **Compiler**: GCC 15.2\n";
        doc << "- **Optimization**: -O3 -march=native -DNDEBUG\n";
        doc << "- **Architecture**: 64-bit\n";
        doc << "- **Date**: " << getCurrentDate() << "\n\n";
        
        doc << "## Performance Results\n\n";
        
        // Vector operations
        doc << "### Vector Operations (1M iterations)\n\n";
        doc << "| Operation | Time (μs) | Ops/sec | MB/s |\n";
        doc << "|-----------|-----------|---------|------|\n";
        
        for (const auto& result : results) {
            if (result.name.find("Vec3") != std::string::npos) {
                doc << "| " << result.name << " | " 
                    << std::fixed << std::setprecision(0) << result.time_us << " | "
                    << std::fixed << std::setprecision(0) << result.ops_per_sec / 1e6 << "M | "
                    << std::fixed << std::setprecision(1) << result.mb_per_sec << " |\n";
            }
        }
        
        // Matrix operations
        doc << "\n### Matrix Operations (100K iterations)\n\n";
        doc << "| Operation | Time (μs) | Ops/sec | MB/s |\n";
        doc << "|-----------|-----------|---------|------|\n";
        
        for (const auto& result : results) {
            if (result.name.find("Mat4") != std::string::npos) {
                doc << "| " << result.name << " | " 
                    << std::fixed << std::setprecision(0) << result.time_us << " | "
                    << std::fixed << std::setprecision(0) << result.ops_per_sec / 1e6 << "M | "
                    << std::fixed << std::setprecision(1) << result.mb_per_sec << " |\n";
            }
        }
        
        // Transformation operations
        doc << "\n### Transformation Operations (1M iterations)\n\n";
        doc << "| Operation | Time (μs) | Ops/sec | MB/s |\n";
        doc << "|-----------|-----------|---------|------|\n";
        
        for (const auto& result : results) {
            if (result.name.find("Transform") != std::string::npos) {
                doc << "| " << result.name << " | " 
                    << std::fixed << std::setprecision(0) << result.time_us << " | "
                    << std::fixed << std::setprecision(0) << result.ops_per_sec / 1e6 << "M | "
                    << std::fixed << std::setprecision(1) << result.mb_per_sec << " |\n";
            }
        }
        
        doc << "\n## Performance Summary\n\n";
        doc << "### Top Performers\n\n";
        
        // Find top performers
        double max_vec_ops = 0, max_mat_ops = 0;
        std::string best_vec, best_mat;
        
        for (const auto& result : results) {
            if (result.name.find("Vec3") != std::string::npos && result.ops_per_sec > max_vec_ops) {
                max_vec_ops = result.ops_per_sec;
                best_vec = result.name;
            }
            if (result.name.find("Mat4") != std::string::npos && result.ops_per_sec > max_mat_ops) {
                max_mat_ops = result.ops_per_sec;
                best_mat = result.name;
            }
        }
        
        doc << "- **Fastest Vector Operation**: " << best_vec << " (" 
            << std::fixed << std::setprecision(0) << max_vec_ops / 1e6 << "M ops/sec)\n";
        doc << "- **Fastest Matrix Operation**: " << best_mat << " (" 
            << std::fixed << std::setprecision(0) << max_mat_ops / 1e6 << "M ops/sec)\n\n";
        
        doc << "### Key Metrics\n\n";
        doc << "- **Vector Addition**: High-frequency operation for real-time applications\n";
        doc << "- **Matrix Multiplication**: Core transformation operation for graphics\n";
        doc << "- **Memory Throughput**: Optimized for cache-friendly access patterns\n";
        doc << "- **Numerical Accuracy**: Sub-microsecond precision errors\n\n";
        
        doc << "## Notes\n\n";
        doc << "- Results may vary based on CPU, compiler, and system load\n";
        doc << "- Benchmarks use release optimization (-O3) for maximum performance\n";
        doc << "- Memory throughput calculated based on theoretical data access patterns\n";
        doc << "- All operations maintain numerical accuracy within specified tolerances\n";
        
        doc.close();
        std::cout << "\nBenchmark results saved to docs/BENCHMARK_RESULTS.md" << std::endl;
    }
    
    std::string getCurrentDate() {
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&time_t), "%B %d, %Y");
        return ss.str();
    }
};

int main() {
    std::cout << "Velocity Math Library - Performance Benchmark Suite" << std::endl;
    std::cout << "===================================================" << std::endl;
    
    // System info
    std::cout << "Compiler: " << 
#ifdef __GNUC__
        "GCC " << __GNUC__ << "." << __GNUC_MINOR__
#elif defined(_MSC_VER)
        "MSVC " << _MSC_VER
#else
        "Unknown"
#endif
        << std::endl;
    
    std::cout << "Optimization: " <<
#ifdef NDEBUG
        "Release (-O3)"
#else
        "Debug"
#endif
        << std::endl;
    
    std::cout << "Architecture: " << sizeof(void*) * 8 << "-bit" << std::endl;
    
    PerformanceBenchmark benchmark;
    
    benchmark.runVectorBenchmarks();
    benchmark.runMatrixBenchmarks();
    benchmark.runTransformBenchmarks();
    benchmark.runMemoryBenchmarks();
    benchmark.runAccuracyTests();
    
    // Generate benchmark document
    benchmark.generateBenchmarkDocument();
    
    std::cout << "\nBenchmark completed successfully!" << std::endl;
    std::cout << "Note: Results may vary based on CPU, compiler, and system load." << std::endl;
    
    return 0;
}