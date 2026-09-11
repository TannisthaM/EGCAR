# Standalone native checks

`native_core.hpp` is exactly the core section of `src/egcar_native.cpp`, preceded
by the system Armadillo include. `reference_020.hpp` preserves the 0.2.0 core in
a separate namespace. Tests here do not exercise Rcpp, R installation or CV.

With system Armadillo, Eigen headers and a C++14 compiler, set EIGEN_INCLUDE to
the directory containing Eigen/Core (commonly /usr/include/eigen3):

```bash
EIGEN_INCLUDE=/usr/include/eigen3
g++ -std=c++14 -O3 -I"$EIGEN_INCLUDE" test_native.cpp -larmadillo -o test_native
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./test_native

g++ -std=c++14 -O3 -I"$EIGEN_INCLUDE" benchmark_021.cpp -larmadillo -o benchmark_021
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./benchmark_021 > timings.csv 2> correctness.json
```

The first program compares 120 cases with the independent dense implementation.
The second compares 864 cases with 0.2.0, followed by an interleaved median-timing
experiment on identical data and solver controls. `test_mapped_context_021.cpp`
additionally checks eight native non-owning-view / input-preservation cases
(compile in the same way as `benchmark_021.cpp`), without exercising Rcpp. `--validate-only` skips its
benchmark. Timings exclude preparation, R/C++ conversion, loading extraction,
CV and parallel startup. For sanitizer checks use `-O1 -g
-fsanitize=address,undefined` instead of `-O3`.

Installing the R package uses RcppArmadillo and RcppEigen's headers instead;
users do not need separate Armadillo or Eigen installations. No executable or
third-party header files are distributed in this repository.
