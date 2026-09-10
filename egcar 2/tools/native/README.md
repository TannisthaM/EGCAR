# Standalone native regression checks

The code between the EGCAR_CORE_HPP guards is identical to the corresponding
part of src/egcar_native.cpp. These tests do not exercise the Rcpp interface.
Run from this directory with system Armadillo headers and a C++ compiler:

```bash
g++ -std=c++14 -O1 -g -fsanitize=address,undefined test_native.cpp -larmadillo -o test_native
OPENBLAS_NUM_THREADS=1 ./test_native
```

System Armadillo is needed only for these standalone tests. Installing the R
package instead uses RcppArmadillo. No executable is distributed.
