# R3.1 deep dive — Weil positivity beyond the log 2 window

**Status 2026-10-22:** historical. Window sizes here are FULL widths (1.75 = half-width 0.875). Superseded by Riemann/Deep/WeilRecord/ (certified half-width 3/2, DISCOVERIES entry 3).

**Verdict: DEAD as a route to RH. Our rigorous window 1.75 is NOT a record** (it beats Zhu's 1.6, but unrefereed 2026 preprints already claim 17/8 and log 15). No DISCOVERY.

Convention: Q(f) = Σ_ρ |f̂(γ_ρ)|²-side functional; f supported in [−L/2, L/2] (window L); only prime powers n < e^L enter.
Q(f) = (1/2π)∫|f̂(r)|² S(r) dr + 2Re(A·conj B), S(r) = Re ψ(1/4+ir/2) − log π − 2Σ_{n<e^L} Λ(n) n^{-1/2} cos(r log n), A = ∫f e^{x/2}, B = ∫f e^{−x/2}.

## Literature (exact windows)
- log 2: Yoshida 1992; Connes–Consani, Selecta 2021 (arXiv 2006.13771). Bombieri, Rend. Lincei 2000 (criterion).
- 1.6: Zhu, arXiv 2608.24827 (v2, 2 Sep 2026), certified, λ_min ∈ [8.9e-18, 2.27e-17]; argues cost is doubly exponential in L.
- 2 and 17/8: V. Liu, alphaXiv 2609 "Certified Weil positivity beyond the unit window" (preprint): 2^-151 at half-width 1, 2^-49162 at 17/16.
- log 15 ≈ 2.708: R. Mori, "paper VIII", Zenodo 10.5281/zenodo.23281175 (2026-10-10), Arb certificate, 2^-2904. Largest claimed window.
- log 3: Angelone (Zenodo 22062590). Desogus, arXiv 2609.20367: claims all windows and so RH; unverified, treat as not credible.
- Numerics for small eigenvalues and zeta zeros: Connes–Consani, "Spectral triples and ζ-cycles" (arXiv 2106.01715). Suzuki, arXiv 2606.09096: framework only.

## λ_min(L) = inf Q(f)/‖f‖² (Legendre Galerkin, Arb 1000 bits; upper bounds; even sector, odd ≈10³× larger)
| L | λ_min | L | λ_min | L | λ_min |
|---|---|---|---|---|---|
| 0.5 | 3.34e-2 | log3 | 5.54e-8 | 1.9 | 2.52e-26 |
| 0.6 | 7.57e-3 | 1.2 | 1.61e-9 | log7 | 6.75e-28 |
| log2−.01 | 1.53e-3 | 1.3 | 3.97e-11 | 2.0 | 5.86e-30 |
| log2 | 1.33e-3 | log4 | 7.64e-13 | log8 | 4.08e-33 |
| log2+.01 | 1.13e-3 | 1.5 | 3.65e-15 | log9 | 2.59e-38 |
| 0.8 | 1.81e-4 | 1.6 | 1.65e-17 | log11 | 9.55e-49 |
| 0.9 | 1.62e-5 | 1.7 | 4.25e-20 | 2.5 | ≤2.9e-53* |
| 1.0 | 9.35e-7 | 1.8 | 4.13e-23 | 3.0 | ≤5.4e-70* |
\* not converged (N=140; true value is smaller). Monotone decreasing, so the smallest margin is always at the largest L. No visible kink where a prime enters. Lowest state is even (as Zhu found).
Checks: explicit formula vs zeta zeros to 1e-28; Galerkin matrix vs Σ_γ|f̂(γ)|² to 15 digits; λ(1.6) sits inside Zhu's enclosure.

## Rigorous certificate (cert/cert.py, Arb ball arithmetic)
Q ≥ M = εI + C: ε = Φ(R) − P, using Φ increasing and S ≥ Φ − P (P = 2ΣΛ(n)n^{-1/2}); C = band-limited part (|r| ≤ R) plus rank-2 pole term. Legendre basis of size N.
Bounded pieces: Legendre tail ‖E‖ ≤ 2‖P⊥C‖_HS (bound |j_n(x)| ≤ x^n/(2n+1)!!); composite Gauss quadrature error (Bernstein ellipse inside |Im r| < 1/2, explicit bound on ψ); node displacement; then Cholesky of the midpoint matrix + residual and radius bounds.
Proven (both parity sectors): **Q(f) ≥ 9.99e-22 ‖f‖² for every complex f supported in [−7/8, 7/8]** (window 1.75, R=800, N=1050, ~8 min, 330 MB). Envelope estimate 1.38e-21 vs true about 1.3e-21, so the bound is nearly sharp.
Also certified: window 1.5 (≥2.0e-15), window 1.0 (≥5.0e-7). Consistency check: M_R ≤ Q against the independent x-space matrix.
Window 1.875 (R=900, N=1230): approximate λ_min(M) = 1.44e-25 > 0, but the certificate failed because the assembled matrix radii came out infinite (a numerical assembly bug, not a mathematical failure; open).

## Exact obstruction to going further with this method
- R must satisfy Φ(R) > P(L), i.e. R ≈ 2π·e^{P(L)}; basis size N ≈ 1.5·(L/2)·R; cost grows like N²·R.
- P jumps at each prime power: 2.94 (n≤4), 4.38 (n≤5, up to log 7), 5.85 (n≤8), 8.52 (n≤12), 9.94 (n≤15).
- Window log 15 needs R ≈ 2·10⁵ and N ≈ 3·10⁵: out of reach (our window 1.75 needed N = 10³). Beating Mori needs a non-pointwise tail method (prolate bases, Schur complement on the band complement).
- Precision is not the bottleneck: λ ≈ 1e-27 needs only 192 bits, because entry radii are about 1e-55.
- Why no window settles RH: λ_min(L) falls doubly exponentially (log λ ≈ −2π²N(T*)/log N(T*), T* = 2π e^L; Zhu), and each window costs exp(c e^{L/2}). No uniformity in L, so this face of RH is finite: it is circular in the sense of R3.7.

## Lean
Riemann/Deep/Weil/Defs.lean (builds, no sorry): weilQ, WeilPos, monotonicity, only n ≤ e^{2a} enter, no prime below log 2. Cited results as Props: WeilCriterion, Yoshida/Connes–Consani, Zhu, Liu, Mori, CertWindow (7/4). Implications Mori → Liu → Cert → Zhu → Yoshida are proved.

## Next steps
1. Keep R3.1 closed as an RH route (NOGO: "finite windows, doubly exponential cost, no uniformity").
2. Optional, small value: Lean-check a certificate (verify the Cholesky residual in Lean) or re-check Mori's log 15 data independently.

## Progress
- M1 literature (Zhu 1.6, Liu 17/8, Mori log 15). M2 conventions checked against zeros. M3 Lean builds. M4 λ_min table. M5 rigorous certificate at 1.75 (cert/). M6 report final.
