# Validation of egcar 0.2.13

Date: 24 September 2026.

The package was compiled and installed using R 4.3.3 on Ubuntu 24.04.3,
x86-64, with g++ 13.3.0. The full installed-package test suite completed:

- 815 expectations passed; zero failures and zero test warnings.
- One opt-in third-party comparison suite was skipped.
- `R CMD check --no-manual`: zero errors, zero warnings, four NOTEs.

The NOTEs concern unavailable optional suggested packages (RhpcBLASctl,
RGCCA, PMA, ggplot2, remotes, roxygen2), the explicit C++14 specification,
installed binary size, and retaining the requested LaTeX source in inst/doc.

The 48 new expectations cover budget controls and old configuration objects,
unlimited-iteration numerical stopping, successful/error/timeout subprocess
outcomes, actual short-budget termination, spawned R descendant cleanup,
supervised vs in-process numerical equality, serial vs parallel SGCA results,
the bundled reference initializer, and public timeout output.

A separate small experiment integration run (n=18, p=(3,3,3), r=1, two
folds, a 0.01-second SGCA budget) confirmed that the experiment result row
contains status=time_limit, converged=FALSE, unavailable subspace error and
the requested timeout message. Its observed SGCA method elapsed time was
0.046 seconds, including startup/cleanup overhead. Other unavailable external
packages were skipped in that integration test. It was not a performance
benchmark or a successful all-methods example.

No six-hour wait, six-hour convergence study, full all-methods example, or
SLURM/CV sweep was performed. The previous 0.2.12 performance figures came
from the dedicated fixed-work benchmark, with three repetitions per case.
They are not new 0.2.13 speed measurements and do not measure SGCA supervision.

The LaTeX write-up passed pdflatex in draft mode without LaTeX errors or
layout-overflow warnings; no PDF was generated. The final delivery adds this
report and the check evidence after validation. The executable source,
namespace, help files, built DESCRIPTION (apart from the build timestamp),
and tests were byte-checked against the package snapshot that passed
R CMD check.

Evidence is included under inst/validation/performance_0213_*.
