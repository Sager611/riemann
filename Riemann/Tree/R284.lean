/-!
# Round 284 — node N1 (Weil positivity): our certified half-width 1.55 still the record; its constants re-checked

* Node: N1 Weil / Connes positivity.
* State of the art: no 2024-2026 result beats our certified half-width 1.55 (2608.24827 certifies an f-window L = 0.8,
  i.e. autocorrelation support 1.6; 2607.02828 and 2606.09096 give finite/operator results); best analytic window remains
  Yoshida / Connes-Consani (2L <= log 2).
* Our record re-checked: explicit-formula consistency PASS (primes and zeros sides agree to 1.8e-31 and 8.1e-29); PROOF.md
  section 1(B) constants recomputed independently (C = a(640) - 9/2 = 0.12359..., margins 3.9e-4 and 6.2e-4; finite tail
  2.82e-4 < 3.59e-4); full Arb re-run not possible for scouts (no python-flint). Extending beyond 1.55 needs the out-of-core
  rewrite noted in REPORT4 (memory 4-5 GB at 1.55).
* Heuristic computation: smallest eigenvalues of a 20-mode Legendre discretisation stay positive but tiny up to L = 2.0
  (5.9e-12 at 1.5, 1.7e-15 at 2.0); cos^k bump families give no useful certificates.
* Exploration ideas: SDP / rational PSD blocks with Bernstein subdivision for Weil kernels (from coding-theory SDP, N1);
  2610.41115 claims function-field Weil-form Hilbert spaces whose self-adjoint extensions have spectra converging to zeta-zero
  ordinates (to audit under N13); finite idele-class quotients as Gram tests (N1).
-/

namespace Fayed.R284
end Fayed.R284
