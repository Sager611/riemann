/-!
# Round 100 — false sufficient conditions for RH and certifiable counterexamples

| conjecture (⇒ RH) | status | Lean certification |
|-------------------|--------|--------------------|
| Turán 1948: ζ_N(s) = Σ_{n≤N} n^{-s} ≠ 0 for Re s > 1 | false: Montgomery 1983 (zeros right of 1 + (4/π − 1)·loglog N / log N for large N); smallest N = 19 (Platt–Trudgian 2016; a late scout confirms none right of 1 for N ≤ 18 and N = 20, 21, 28) | TARGET: one zero of ζ_19 in Re s > 1 (finite Dirichlet polynomial + box count), needs its location |
| Pólya 1919: L(x) = Σ_{n≤x} λ(n) ≤ 0 | false: Haselgrove 1958; smallest x = 906,150,257 (Tanaka 1980) | SKIP: kernel sieve to 9·10⁸ is out of reach without native_decide |
| Turán: T(x) = Σ_{n≤x} λ(n)/n ≥ 0 | false: Haselgrove 1958; smallest ≈ 7.2·10¹³ (Borwein–Ferguson–Mossinghoff 2008, recalled) | SKIP |
| Mertens: |M(x)| < √x | false: Odlyzko–te Riele 1985; no explicit x; first counterexample < exp(3.21·10⁶⁴) (Pintz 1987); none below 10¹⁶ (Hurst 2018) | SKIP |
| Jensen polynomials of ξ (Griffin–Ono–Rolen–Zagier 2019) | hyperbolic for each degree d and large shifts; no false finite statement | SKIP |

Two scouts (smallest N with location; Turán T(x)) timed out.
-/
