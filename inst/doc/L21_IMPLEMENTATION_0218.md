# Global row L21 with exact spectral compression

Write the directed off-diagonal coefficient blocks as C_kl. The objective is

\[
\frac14\sum_{k\ne l}\operatorname{tr}(C_{kl}^T S_{kk} C_{kl}S_{ll})
-\frac12\sum_{k\ne l}\langle S_{kl},C_{kl}\rangle
+\lambda\sum_i\|C_{i,:}\|_2.
\]

The half weight per directed pair matches the earlier smooth loss when C is symmetric. Row and column norms can differ during directed fitting. After the sparse ADMM copy Z is complete, return (Z+Z^T)/2. Convergence certifies the directed solve, not optimality of the final projection for the earlier constrained objective.

For a right covariance basis Q_l, let R_kl=S_kl Q_l and N_kl=S_kl-R_kl Q_l^T. Every row of N_kl is orthogonal to Q_l in exact arithmetic. Represent the compact DRS/ADMM state H=Z+U as

\[
H_{kl}=B_{kl}Q_l^T+\operatorname{diag}(h_k)N_{kl}.
\]

A single h_ki works across every outgoing block l of row i: cold states start at zero, and the complementary update adds the same 1/(2*mu) times that row of N_kl in each block. Row proximal maps and dual rescaling apply the same row multiplier to every block. This invariant remains true for subsequent warm starts from the representation. Genuine nonzero N_kl is retained, including in the tests.

Consequently, each global squared row norm is the sum over outgoing blocks of ||B_kl[i,]||^2 plus h_ki^2 times the sum of ||N_kl[i,]||^2. The latter base norms are cached as vectors. No complementary coefficient matrix is carried through iterations.

Let a be the previous row shrinkage multiplier. Then Z=aH and U=(1-a)H. Form T=(2a-1)H. In right coordinates, project only on the left:

\[
P_t=Q_k^T T_B,\qquad
C_t=(S_{t,kl}+2\mu P_t)/(d_kd_l^T+2\mu).
\]

For a full left basis, C_B=Q_k C_t. For a thin left basis,

\[
C_B=T_B+(R_{kl}-Q_kS_{t,kl})/(2\mu)+Q_k(C_t-P_t).
\]

The nullspace coefficient is C_h=(2a-1)h+1/(2*mu). These expressions retain both left and right complementary terms. All comparisons with the independent explicit R solver include full-space residuals and objective histories.

Use alpha=1.6 to form W=alpha*C+(1-alpha*a)*H and shrink every global row by

\[
a_{new}=\max(0,1-\lambda/(\mu\|W_{i,:}\|_2)).
\]

The unrelaxed primal residual C-Z_new and dual residual mu*(Z_new-Z_old) are recovered from accumulated row sums of W^2, (W-H_old)^2 and W*(W-H_old). Stop thresholds use the ORIGINAL full-space dimension 2q. Relative residuals control adaptive mu; this does not replace the stopping criteria or change lambda. Dual rescaling preserves Z exactly: multiply H and h by a+(1-a)/factor and divide a by that quantity.

After convergence or the iteration cap, apply the existing numerical zero thresholds to the directed proximal rows, reconstruct each unordered symmetric edge as the average of the two directed edges, and compute loadings through the existing partial EVD. Symmetry is implicit in the one-block-per-unordered-pair representation. No full P-by-P reconstruction or full EVD is introduced.

## Functions changed or added

- `fit_l21_admm`: chooses the compact native path, constructs diagnostics and exposes symmetric blocks. Its timer is inside the existing complete-fit timer.
- `solve_rotated_l21` (C++): spectral/complementary ADMM, residuals, adaptive mu, final block symmetry, native cold allocation and state release.
- `fit_l21_full_admm` / `egcar_l21_direct_R`: independent full-space directed diagnostic path and compatibility path for directed full-coordinate states.
- `egcar_l21_expand_rotated`: reconstructs full H only for cross-backend diagnostics/warm starts; it is not used inside the fast loop.
- `.egcar_validate_init` / `egcar_warm_state`: validate and retain Hk/Hl in right spectral coordinates plus a and hnull, tagged `state_basis="right_spectral"` and `formulation="global_rows_posthoc"`.
- `fit_l21_admm_legacy`: internal, preserved symmetry-constrained baseline, not selected by new public L21 calls.

The covariance generator, L11/benchmark algorithms, loading algorithm and error metric implementations are unchanged. The rate coefficient is inherited rather than recalibrated. See VALIDATION.md for empirical runtime and total-memory limitations.
