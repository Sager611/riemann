FORMALIZED: Fasi–Higham's conjecture for m = 3 and every k — every Taylor coefficient of 1/q_{k,3} at 0 is positive (FasiLean.fasi_higham_m3, in Mathlib's PowerSeries ℚ), plus the quantitative form c_n ≥ (R₀ − 2|R₊|) z₀^{−n} with 2|R₊| < R₀ (FasiLean.c_quant, real form); from NO hypotheses (the real root d is constructed in Lean from sSup on ℝ), sorry-free, standard axioms only.

# Lean: Fasi–Higham positivity, m = 3, every k

Discovery: [discoveries/fasi-higham-m3.md](../../../discoveries/fasi-higham-m3.md); proof: [NOTE.md](NOTE.md) §2. Files: [Lean/](Lean/) — Defs, Core, KeyIneq, Bound, Root, PS, Main, MainPS, Quant, Quant2, Check (≈ 760 lines). Namespace FasiLean. Lean 4 v4.35.0-rc2, Mathlib v4.35.0-rc2.

## Main statements (verbatim)

    -- Defs.lean: Fasi–Higham's q_{k,m} coefficient and the recurrence for c_n(k)
    def qcoef (k m j : ℕ) : ℚ :=
      if j ≤ m then (-1) ^ j * (m.choose j : ℚ) * ((k + m - j).factorial : ℚ) / ((k + m).factorial : ℚ)
      else 0
    def cs (k : ℕ) : ℕ → ℚ          -- cs k 0 = c_{−2}, cs k 1 = c_{−1}, cs k 2 = c_0
      | 0 => 0
      | 1 => 0
      | 2 => 1
      | n + 3 => 3 / (k + 3) * cs k (n + 2) - 3 / ((k + 2) * (k + 3)) * cs k (n + 1)
          + 1 / ((k + 1) * (k + 2) * (k + 3)) * cs k n
    def c (k n : ℕ) : ℚ := cs k (n + 2)
    theorem conv (k n : ℕ) :
        ∑ j ∈ range (n + 1), qcoef k 3 j * c k (n - j) = if n = 0 then 1 else 0
    -- PS.lean: c is the coefficient sequence of the power-series inverse of q_{k,3}
    noncomputable def qPS (k m : ℕ) : PowerSeries ℚ := PowerSeries.mk (qcoef k m)
    theorem coeff_inv (k n : ℕ) : PowerSeries.coeff n (qPS k 3)⁻¹ = c k n
    theorem conv_unique (k : ℕ) (b : ℕ → ℚ)
        (hb : ∀ n, ∑ j ∈ range (n + 1), qcoef k 3 j * b (n - j) = if n = 0 then 1 else 0) : b = c k
    -- Main.lean, MainPS.lean: MAIN THEOREM
    theorem c_pos (k n : ℕ) : 0 < c k n
    theorem fasi_higham_m3 (k n : ℕ) : 0 < PowerSeries.coeff n (qPS k 3)⁻¹
    -- Root.lean: the root of NOTE §2 (proved, no IVT)
    theorem exists_root (N : ℝ) (hN : 1 ≤ N) :
        ∃ d : ℝ, 1 / 2 < d ∧ d < 2 / 3 ∧ d ^ 3 + 3 * N * d = 2 * N
    -- Quant.lean, Quant2.lean: quantitative form
    noncomputable def qR (k : ℕ) (z : ℝ) : ℝ := ∑ j ∈ range 4, (qcoef k 3 j : ℝ) * z ^ j
    noncomputable def qR' (k : ℕ) (z : ℝ) : ℝ := ∑ j ∈ range 4, (j : ℝ) * (qcoef k 3 j : ℝ) * z ^ (j - 1)
    theorem c_quant (k : ℕ) (N : ℝ) (hN : N = k + 1) :
        ∃ d z₀ u v R₀ K : ℝ,
          (1 / 2 < d ∧ d < 2 / 3 ∧ d ^ 3 + 3 * N * d = 2 * N) ∧
          z₀ = N + d ∧ qR k z₀ = 0 ∧
          (∀ z, N * (N + 1) * (N + 2) * qR k z = (z₀ - z) * (z ^ 2 - 2 * u * z + v)) ∧
          u ^ 2 < v ∧ z₀ ^ 2 < v ∧
          R₀ = -1 / (z₀ * qR' k z₀) ∧ 0 < R₀ ∧
          K = R₀ ^ 2 * (8 * N * (N + d) ^ 2 * (1 - d) ^ 2
            / (d ^ 3 * ((N + d) ^ 2 + 3 * N * (1 - d)) * (3 * N + d ^ 2 + d + 1))) ∧
          0 ≤ K ∧ K < R₀ ^ 2 ∧
          (∀ n, ((c k n : ℝ) - R₀ * (1 / z₀) ^ n) ^ 2 ≤ K * ((1 / z₀) ^ n) ^ 2) ∧
          (∀ n r, 0 ≤ r → K ≤ r ^ 2 → (R₀ - r) * (1 / z₀) ^ n ≤ (c k n : ℝ))

Reading of c_quant: K is NOTE §2's (2|R₊|)² (its closed form of (2|R₊|/R₀)² times R₀²), so with r = √K = 2|R₊| < R₀ the last clause is c_n ≥ (R₀ − 2|R₊|) z₀^{−n} > 0; the clause before it is the stronger two-sided bound |c_n − R₀z₀^{−n}| ≤ 2|R₊| z₀^{−n}. The factorization with u² < v says z₀ is the only real zero and the complex pair z_± has |z_±|² = v > z₀² (D1).

## Hypotheses

None. The brief allowed one named hypothesis for the root (∀ N ≥ 1, ∃ d ∈ (1/2, 2/3), d³ + 3Nd = 2N; intended proof Mathlib's intermediate_value_Icc); it is proved instead (Root.lean, exists_root: d = sSup {x ≤ 2/3 : x³ + 3Nx ≤ 2N}, with explicit polynomial continuity estimates below/above and F_mono), so no uncompiled instance file is needed.

## Proof route (NOTE §2 in real variables; no ℂ, no square roots, no IVT)

- N = k + 1, a = z₀ = N + d, u = N − d/2, v = N² − dN + d² + 3N, R = R₀ = v/(3(d² + N)), ρ = 1/a. Exact identities: N(N+1)(N+2) = av (uses d³ + 3Nd = 2N), 3N(N+1) = 2ua + v, 3N = a + 2u, a² − 2ua + v = 3(d² + N); so N(N+1)(N+2)q(z) = (a − z)(z² − 2uz + v), v − u² = (3/4)d² + 3N > 0, v − a² = 3N(1 − d) > 0 (D1).
- Core.lean, decomp: with the shifted sequence C (C 0 = c_{−2} = 0, C 1 = c_{−1} = 0, C 2 = c_0 = 1), C_m = R a² ρ^m + f_m for all m, where f_0 = −Ra², f_1 = −Ra, v f_{m+2} = 2u f_{m+1} − f_m (f_m = c_{m−2} − Rρ^{m−2}, the remainder extended to n = −2, −1); one linear_combination per step.
- Gf_succ, Gf_pow, sq_le_Gf: G_m = v f_{m+1}² − 2u f_{m+1} f_m + f_m² satisfies v G_{m+1} = G_m, so G_m v^m = G_0 = R²a²·3(d² + N); and (v − u²) f_m² ≤ v G_m (completing the square: v G_m − (v − u²)f_m² = (v f_{m+1} − u f_m)²).
- Bound.lean, core: (c_n − Rρ^n)² = f_{n+2}² ≤ K ρ^{2n}, K = R²·3a²(d² + N)/(v(v − u²)), using v^n ≥ a^{2n}; K < R² is KeyIneq.lean's a²(a² − 2ua + v) < v(v − u²).
- KeyIneq.lean, key_ineq: v(v − u²) − 3a²(d² + N) = (9/4)((4 − 4d − d²)N² + (d² − 3d³)N − d⁴) > (9/4)((8/9)N² − (4/9)N − 16/81) > 0 for every N ≥ 1, d ∈ (1/2, 2/3).
- Main.lean, c_pos: pos_of_bound ((x − Rw)² ≤ Kw², K < R² ⇒ x > 0) with w = ρ^n, then exact_mod_cast to ℚ.
- Quant.lean, K_note: 3a²(d² + N)/(v(v − u²)) = NOTE's 8N(N+d)²(1−d)²/(d³[(N+d)² + 3N(1−d)][3N+d²+d+1]) on the curve: one linear_combination with cofactor 3(N+d)²v(Nd² − 2Nd + 4N + d⁴ + d³ + d²) (found with sympy). Quant2.lean: q(z₀) = 0, the factorization, R₀ = −1/(z₀q'(z₀)) from N(N+1)(N+2)q'(z₀) = −3(d² + N), and lower_of_bound.
- Simplification found while formalizing (worth a line in NOTE §2): extending the remainder to n = −2, −1 gives G_0 in closed form, (2|R₊|/R₀)² = 4(N+d)²(d²+N)/(v(d² + 4N)), and 2|R₊| < R₀ becomes (4 − 4d − d²)N² + (d² − 3d³)N − d⁴ > 0, which holds for ALL N ≥ 1, d ∈ (1/2, 2/3) without the cubic relation; NOTE's (2 − 3d)³-identity and g(d) < 0 are not needed. The relation enters only through N(N+1)(N+2) = z₀v (the factorization).

## Build

Each module separately, one at a time, lean directly (no lake process, no lake build), oleans to .scratch/fasilean/olean (build output, deleted at the end); LEAN_PATH = lake's (as in Riemann/Deep/DPCUncond/Lean/audit/leanpath.txt). Per module m in the order Defs Core KeyIneq Bound Root PS Main MainPS Quant Quant2 (d = Riemann/Deep/FasiHigham/Lean), then Check without -o:

    LEAN_PATH=.scratch/fasilean/olean:$(cat Riemann/Deep/DPCUncond/Lean/audit/leanpath.txt) bash .scratch/memcap.sh 290 1800 \
      /usr/bin/time -l ~/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean -j1 -o .scratch/fasilean/olean/$d/$m.olean $d/$m.lean

Clean rebuild of all 11 modules (Riemann/Deep/FasiHigham/Lean/audit/verify.sh, 12 Nov 2026 17:15–17:22 UTC), all rc = 0; peak memory footprint (/usr/bin/time -l) and memcap's sampled top-mem peak:

| module | imports (Mathlib) | time | footprint | memcap |
|---|---|---|---|---|
| Defs | Algebra.Order.Field.Rat, Data.Nat.Choose.Basic, BigOperators.Group.Finset.Basic, tactics | 49 s | 252 MiB | 228 MB |
| Core | Algebra.Order.Archimedean.Real.Basic, tactics | 25 s | 260 MiB | 248 MB |
| KeyIneq | same | 23 s | 262 MiB | 245 MB |
| Bound | Core, KeyIneq | 29 s | 268 MiB | 254 MB |
| Root | Algebra.Order.Archimedean.Real.Basic, tactics | 34 s | 279 MiB | 268 MB |
| PS | Defs, RingTheory.MvPowerSeries.Inverse, RingTheory.PowerSeries.Basic | 46 s | 282 MiB | 271 MB |
| Main | Defs, Bound, Root | 39 s | 252 MiB | 240 MB |
| MainPS | Main, PS | 26 s | 280 MiB | 215 MB |
| Quant | Main | 30 s | 274 MiB | 244 MB |
| Quant2 | Quant | 59 s | 273 MiB | 244 MB |
| Check | MainPS, Quant2 | 83 s | 273 MiB | 263 MB |

No Mathlib import beyond these (never import Mathlib; no topology, rpow, log, ℂ or IVT); nlinarith only once (Core.sq_le_Gf, one hint). Logs: Riemann/Deep/FasiHigham/Lean/audit/peaks.log (every compile of the session), logs/log_<module>.txt, review_prompt.txt, review_out.md, STATE.md.

## Static scan

    rg -n -w 'sorry|admit|native_decide|skipKernelTC' Riemann/Deep/FasiHigham/Lean   → no match (rc = 1)
    rg -n '^\s*axiom\s' Riemann/Deep/FasiHigham/Lean                                → no match (rc = 1)

## #print axioms (Check.lean)

    'FasiLean.fasi_higham_m3' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.c_pos' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.coeff_inv' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.conv' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.conv_unique' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.qcoef_zero' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.qcoef_one' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.qcoef_two' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.qcoef_three' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.qcoef_big' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.exists_root' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.decomp' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.key_ineq' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.core' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.main_real' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.qR_eq' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.qR'_eq' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.K_note' depends on axioms: [propext, Classical.choice, Quot.sound]
    'FasiLean.c_quant' depends on axioms: [propext, Classical.choice, Quot.sound]

## Fidelity (each Lean object against the source and NOTE §2)

Source: Scardino, arXiv 2610.39051, e-print TeX conjecture_higham_arxiv.tex (fetched 12 Nov 2026, deleted after reading): l.50 (eq:qdef) q_{k,m}(z) = Σ_{j=0}^{m} (−1)^j C(m,j) (k+m−j)!/(k+m)! z^j; l.46 q_{k,m}(0) = 1, degree ≤ m; l.67–69 1/q_{k,m}(z) = Σ_{n≥0} c_n(k,m) z^n, the Taylor expansion at the origin; l.74–76 Conjecture (Fasi–Higham Conj. 3.2): k ≥ 0, m ≥ 1 integers, m odd ⇒ all coefficients of the Taylor expansion of q_{k,m}^{−1}(z) at 0 are positive.

| Lean | source | match |
|---|---|---|
| qcoef k m j | TeX l.50 | exact: (−1)^j C(m,j)(k+m−j)!/(k+m)! for j ≤ m (no ℕ-truncation, j ≤ m ≤ k+m), 0 for j > m (degree m). qcoef_zero…three, qcoef_big: 1, −3/(k+3), 3/((k+2)(k+3)), −1/((k+1)(k+2)(k+3)), 0; i.e. the brief's (k+1)(k+2)(k+3)q = (k+1)(k+2)(k+3) − 3(k+1)(k+2)z + 3(k+1)z² − z³ (Quant.qR_eq) |
| qPS k m, (qPS k 3)⁻¹ | TeX l.46, l.67–69 | q_{k,m} as PowerSeries.mk of its coefficients; Mathlib's inverse in ℚ⟦X⟧ (constant coefficient 1, qPS_const); for a polynomial with q(0) ≠ 0 these are the Taylor coefficients of the analytic 1/q at 0 |
| cs, c, conv, conv_unique, coeff_inv | brief's recurrence, c_{<0} = 0 | conv: Σ_{j≤n} q_j c_{n−j} = [n = 0]; conv_unique: the only such sequence; coeff_inv: c k n = coeff n (qPS k 3)⁻¹ |
| fasi_higham_m3, c_pos | TeX l.74–76 | exact for m = 3: every k, n ∈ ℕ (including k = 0, n = 0); no claim for other odd m |
| exists_root | NOTE §2: d the root of d³ + 3Nd − 2N in (1/2, 2/3) | exact, every real N ≥ 1 |
| z₀, u, v in c_quant | NOTE §2: z₀ = N + d, z_± = N − d/2 ∓ i(√3/2)s, s² = d² + 4N | u = Re z_± = N − d/2, v = (N − d/2)² + (3/4)(d² + 4N) = N² − dN + d² + 3N = \|z_±\|²; the proved factorization N(N+1)(N+2)q(z) = (z₀ − z)(z² − 2uz + v) with u² < v makes z_± the roots of the quadratic; z₀² < v is D1 (\|z_±\|² − z₀² = 3N(1 − d)) |
| R₀ | NOTE: R_j = −1/(z_j q'(z_j)), c_n = Σ_j R_j z_j^{−n} | stated literally with qR' (termwise derivative; its j = 0 term is 0·…); R₀ = v/(3(d² + N)) = Q(0)d/(6N(N+d)(1−d)) |
| K | NOTE: (2\|R₊\|/R₀)² = 8N(N+d)²(1−d)²/(d³[(N+d)² + 3N(1−d)][3N + d² + d + 1]) | K = R₀² × NOTE's expression, i.e. K = (2\|R₊\|)²; K < R₀² ⟺ 2\|R₊\| < R₀ |
| last clause of c_quant | NOTE §1.1: c_n ≥ (R₀ − 2\|R₊\|) z₀^{−n}, R₀ − 2\|R₊\| > 0 | with r = √K = 2\|R₊\| (0 ≤ r < R₀) it is exactly this bound; sqrt-free equivalent form; the clause before it is the stronger \|c_n − R₀z₀^{−n}\| ≤ 2\|R₊\| z₀^{−n} |
| normalisation | brief: z^{−n} vs z^{−n−1} | re-derived: P/((z₀ − z)S(z)) has the term P/(S(z₀)(z₀ − z)) = Σ_n P/(S(z₀)z₀) (z/z₀)^n, so the coefficient is R₀ z₀^{−n} with R₀ = P/(z₀ S(z₀)) = −1/(z₀ q'(z₀)) (q'(z₀) = −S(z₀)/P, S(z₀) = 3(d² + N)); decomp's base case C 2 = 1 = R₀ + f_2 checks it |

Reviewer (statement fidelity; reviewer subagent, 12 Nov 2026; prompt Riemann/Deep/FasiHigham/Lean/audit/review_prompt.txt with the Lean statements, the TeX definition and NOTE §2 inline; output review_out.md): VERDICT FAITHFUL — qcoef incl. the j > m cut-off and ℕ-subtraction; Taylor identification incl. n = 0 and the z^n normalisation; scope exactly m = 3, all k, n; qR' at j = 0 vanishes by its factor 0; R₀ = −1/(z₀q'(z₀)) is the coefficient of 1/(1 − z/z₀); c_quant's last clause with r = √K is exactly c_n ≥ (R₀ − 2|R₊|)z₀^{−n} with positive margin; no hidden hypothesis or vacuity.

## FORMALIZATIONS.md entry

- Fasi–Higham positivity for m = 3 and every k (our side discovery [fasi-higham-m3](discoveries/fasi-higham-m3.md)) machine-checked in Lean 4 + Mathlib: FasiLean.fasi_higham_m3 : 0 < PowerSeries.coeff n (qPS k 3)⁻¹ for all k, n (q_{k,3} from the general q_{k,m} formula, Mathlib's power-series inverse), plus the quantitative form c_quant (c_n ≥ (R₀ − 2|R₊|)z₀^{−n}, 2|R₊| < R₀, real form); no hypotheses (the root is built from sSup on ℝ), sorry-free, standard axioms, no native_decide. [LEAN.md](Riemann/Deep/FasiHigham/LEAN.md) Reviewer: FAITHFUL.

## .lake

Left as it is: the repo-root .lake is the real Mathlib cache shared with the other Lean agents; not deleted, no lake exe cache get, no lake build, placeholder not restored. The orchestrator restores the placeholder.
