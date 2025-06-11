/*
The MIT License (MIT)

Copyright (c) 2024 Elowyn Fearne

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef ACOREX_ENABLE_BENCHMARK_MAIN
#define ACOREX_ENABLE_BENCHMARK_MAIN 0
#endif

#if ACOREX_ENABLE_BENCHMARK_MAIN
// Standalone benchmark runner
// Compile with: g++ -std=c++17 -O3 -march=native RunBenchmarks.cpp Explorer/SIMDUtils.cpp Explorer/BenchmarkSIMD.cpp -o run_benchmarks

#include <iostream>

// Forward declaration
void runSIMDBenchmarks();

int main(int argc, char* argv[]) {
    std::cout << "Starting SIMD Benchmarks for ACorEx Audio Processing..." << std::endl;
    
    try {
        runSIMDBenchmarks();
    } catch (const std::exception& e) {
        std::cerr << "Error running benchmarks: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
#endif // ACOREX_ENABLE_BENCHMARK_MAIN