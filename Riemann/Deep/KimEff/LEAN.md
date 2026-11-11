PARTIAL: proved in Lean (sorry-free, standard axioms, no native_decide): DERIVATION §2.1 and §3 (zeros of Θ in the closed convex parabola, from ζ's half-plane via Θ(z²) = 8ξ(1/2+z) and the functional equation; the sector bound; the Obreschkoff limit step with Hurwitz for real polynomials) and all of §10 ((ii) ⇒ d ≤ 6.45 n²/log² n, the four (c, C) pairs, n_* = ⌈Z(12)⌉, the two-thirds wedge inside for d > 1.04710e8), γ real, and the a = 1/2 variant (3.59) with no zero-location input; hypotheses: (i) itself (§4–§9 + certify.py), Obreschkoff–Kim–Lee Thm 3, Hadamard/Weierstrass partial products, Gauss–Lucas for genus 0, γ(m) > 0, Bebel's half-plane.

# KimEff in Lean 4 + Mathlib

Files: [Lean/](Lean/) (1,430 lines): Core (definitions, Hurwitz, geometry), Step (named hypotheses, §3), Arith (§10: H2, H3, c = 1), Theorem (abstract main theorem), Xi (Riemann's ξ), Power (c = 0.51, 0.55, 0.6), Half (a = 1/2), NStar (n_* = ⌈Z(12)⌉), Real (γ real), TwoThirds (old record inside), Main (imports all), Axioms. Lean v4.35.0-rc2, Mathlib 0653561 (lake-manifest.json).

## Main theorem

```lean
theorem KimEff.kimEff_main (hH : HalfPlane a0) (hO : ObreschkoffKimLee)
    (hHad : ∀ n, HadamardApprox xiGamma ThetaDeriv n)
    (hGL : ∀ n, GaussLucasGenus0 ThetaDeriv n) (hpos : ∀ m, 0 < xiGamma m)
    (hi : ZeroFreeRegion a0 ThetaDeriv) :
    (∀ n : ℕ, nStar ≤ n → ∀ u : ℝ, 0 < u → Zfun u = n → ∀ d : ℕ,
        (d : ℝ) ≤ Dfun a0 (Rfun n u) → Hyperbolic (jensen xiGamma d n)) ∧
    (∀ n : ℕ, nStar ≤ n → ∀ d : ℕ, (d : ℝ) ≤ 6.45 * (n : ℝ) ^ 2 / Real.log n ^ 2 →
        Hyperbolic (jensen xiGamma d n)) ∧
    (∀ n d : ℕ, nStar ≤ n → 7.6890 * (d : ℝ) ^ (0.51 : ℝ) ≤ n → Hyperbolic (jensen xiGamma d n)) ∧
    (∀ n d : ℕ, nStar ≤ n → 0.87160 * (d : ℝ) ^ (0.55 : ℝ) ≤ n → Hyperbolic (jensen xiGamma d n)) ∧
    (∀ n d : ℕ, nStar ≤ n → 0.069829 * (d : ℝ) ^ (0.6 : ℝ) ≤ n → Hyperbolic (jensen xiGamma d n)) ∧
    (∀ n d : ℕ, nStar ≤ n → 1.1853e-10 * (d : ℝ) ≤ n → Hyperbolic (jensen xiGamma d n))
```

Also proved: xi_jensen_hyperbolic_half (a = 1/2: no hH; (ii) and d ≤ 3.59 n²/log² n); jensen_hyperbolic_ii/iii for any real γ, any F n (in the role of Θ^{(n)}) and any a > 0 (bundle Inputs); Zfun_12_bounds (998616377770 < Z(12) ≤ 998616377771) and twelve_le_u_iff (u > 0, Z(u) = n: u ≥ 12 ⇔ n ≥ n_*); Zfun_strictMonoOn (u(n) unique), exists_u; xiGammaC_eq, xiF_eq_tsum (8ξ(1/2+z) = Σ γ(m) z^{2m}/m! with γ real); twoThirds_sub, jensen_hyperbolic_twoThirds (d > 1.04710e8 and n² log(n+2) ≥ 24 d³ ⇒ n > n_* and d ≤ 6.45 n²/log² n, so the previous record lies inside).

## Definitions

| Lean | DERIVATION |
|---|---|
| Hyperbolic p | all complex zeros of the real polynomial p are real (p = 0 is not hyperbolic) |
| jensen γ d n, jensenOf P d | J^{d,n} = Σ_{j≤d} C(d,j) γ(n+j) X^j; J(P; d) = Σ_{j≤d} C(d,j) P^{(j)}(0) X^j |
| riemannXi s | (s(s−1)Λ₀(s) + 1)/2 = s(s−1)π^{−s/2}Γ(s/2)ζ(s)/2 (riemannXi_eq_completed; Λ₀ = Mathlib's completedRiemannZeta₀) |
| xiGamma m | m!/(2m)! · (d/dz)^{2m}[8ξ(1/2+z)] at 0 (Re taken; imaginary part proved 0) |
| ThetaDeriv n w | Σ_k γ(n+k) w^k/k! = Θ^{(n)}(w); ThetaDeriv_zero_sq: Θ(z²) = 8ξ(1/2+z) |
| parabola a, wRegion a R, sector δ | P = {(Im w)² ≤ 4a²(a² − Re w)}, W_R = P ∩ {Re w ≥ −R}, S(δ) = {\|Im z\| ≤ δ\|z\|} |
| sectorDelta a R, Dfun a R | δ(R) = 2a√(R+a²)/R, D(R) = R²/(4a²(R+a²)) (Dfun_mul_sq: D δ² = 1) |
| Zfun, Rfun n u, nStar, a0 | Z(u) = u(πe^{2u} − 9/4), (9/10)(n − 1/2)²/(u+1)², 998616377771, 0.37338196661320752 |

## Hypotheses

| name | statement | source | why not proved |
|---|---|---|---|
| HalfPlane a0 | ζ(s) ≠ 0 for Re s > 1/2 + a0 = 0.87338196661320752 (zerosSq_xi then proves \|Re ρ − 1/2\| ≤ a0 via ξ(1−s) = ξ(s)) | (I2): Bebel, QuasiRiemannTracker, Lean-checked (B_2 = 0.87338196661320751993 < 1/2 + a0) | external Lean development, not ported; the repo's port (OpenAI 7/8) gives only a = 3/8 > a0. With a = 1/2 Mathlib proves it (halfPlane_half) |
| ObreschkoffKimLee | P real, P(0) ≠ 0, all zeros in S(δ), d δ² ≤ 1 ⇒ J(P; d) hyperbolic | (I3): Obreschkoff 1961 as Kim–Lee arXiv 2105.05386 Thm 3 and corollary (TeX l.168–189, R423); P(0) ≠ 0 only weakens it | polynomial theorem of its own (Obreschkoff's sector theorem); not in Mathlib |
| HadamardApprox γ F n | ∃ real P_k, zeros of P_k ⊆ zeros of Θ^{(n)}, P_k^{(j)}(0) → γ(n+j) ∀ j | (2.2) Hadamard (order 1/2, genus 0) + §3 partial products and Weierstrass; (I5) | Hadamard factorization is not in Mathlib |
| GaussLucasGenus0 F n | zeros of Θ^{(n)} ⊆ closure(conv(zeros of Θ)) | (2.3) | needs the genus-0 product and its logarithmic derivative; Mathlib has polynomial Gauss–Lucas only |
| gammaPos | γ(m) > 0 | (2.2), Riemann's integral (I1) with φ > 0 (5.1) | Riemann's integral (I1) is not in Mathlib |
| ZeroFreeRegion a0 F | (i): n ≥ n_*, u > 0, Z(u) = n, w ∈ W_{R(n,u)} ⇒ Θ^{(n)}(w) ≠ 0 (u = u(n) ≥ 12: Zfun_strictMonoOn, twelve_le_u_iff) | DERIVATION §4–§9 + certify.py (Arb, u_* = 12) | needs Riemann's integral (I1) for the mixture identity, the Bessel facts (I4: DLMF 10.21.15, 10.16.9, 10.21.3), and Laplace bounds with about 15 certified majorants monotone in u: weeks of Lean. The §8 numbers at u = 12 (M1, Q2, η, Φ, E+, c6, FINAL) are not re-checked in Lean; they only feed (i) |

Each hypothesis is a true statement about the actual ξ (Θ^{(n)} is the n-th derivative of the entire Θ; γ_C = γ is proved), so the main theorem is not vacuous. HadamardApprox and GaussLucasGenus0 are consequences of DERIVATION (2.2)–(2.3), stated in the form §3 uses.

## What the Lean proofs cover

- Hurwitz, polynomial case (hyperbolic_of_tendsto): coefficientwise limit of hyperbolic real polynomials of degree ≤ d with nonzero d-th coefficient is hyperbolic (|p(z)| ≥ |lc| |Im z|^d).
- §2.1: z² ∈ P for \|Re z\| ≤ a; P convex and closed; closure(conv Z(Θ)) ⊆ P.
- §3: w ∈ P, Re w < −R ⇒ \|Im w\| ≤ δ(R)\|w\|; D(R)δ(R)² = 1; the limit argument (P_k(0) ≠ 0 eventually, Obreschkoff, Hurwitz).
- §10: H2 (u(n) + 1 ≤ (log n)/2 for n ≥ n_*, by IVT: (log n/2 − 1)(π n e^{−2} − 9/4) ≥ n), D(R) ≥ (R − a²)/(4a²), H3 (with 3.6/(4a0²) ≥ 6.4555, log² n ≤ 2n), power forms via c = p/q and log^{2p} n ≤ 6.45^p C^q n^{2p−q} (x ≤ e^{x−1}, e > 2.7182818283, log n_* ≤ 27.62964 from Taylor lower bounds; exact rationals up to ~1,000 digits).
- n_*: e to 18 digits (Real.exp_bound, 20 terms), π to 20 digits (Mathlib), Z(12) ∈ (998616377770, 998616377771].
- Two-thirds wedge: n ≤ n_* would give n² log(n+2) ≤ n_*² · 27.62964 ≤ 24 · (1.04710e8)³; then log(n+2) ≤ 2n, log⁶ n ≤ 8n³.
- ξ: Θ(z²) = 8ξ(1/2+z) (Taylor series of the entire even 8ξ(1/2+z)); zeros of ξ have \|Re s − 1/2\| ≤ a from the half-plane and ξ(1−s) = ξ(s); γ real (ζ(s̄) = conj ζ(s), Γ real on ℝ, derivatives of a function real on an interval are real).

## Build and check

```
rm .lake                      # placeholder file; restored afterwards
lake exe cache get            # Mathlib oleans
lake build Riemann.Deep.KimEff.Lean.Main
lake env lean Riemann/Deep/KimEff/Lean/Axioms.lean
rg -n "sorry|admit|native_decide|axiom |implemented_by|extern" Riemann/Deep/KimEff/Lean/   # no hits
```

```
✔ Built Riemann.Deep.KimEff.Lean.Main — Build completed successfully (3300 jobs); 0 warnings in every KimEff module
build times (one module at a time): Core 149 s, Step 29, Arith 33, Theorem 132, Xi 120, Power 182, Half 254,
NStar 127, Real 172, TwoThirds 116, Main 227 (≈ 26 min sequential, after the 4 min Mathlib cache)
lake env lean Riemann/Deep/KimEff/Lean/Axioms.lean:
'KimEff.kimEff_main' depends on axioms: [propext, Classical.choice, Quot.sound]
'KimEff.xi_jensen_hyperbolic_half' depends on axioms: [propext, Classical.choice, Quot.sound]
'KimEff.xi_jensen_hyperbolic' depends on axioms: [propext, Classical.choice, Quot.sound]
'KimEff.jensen_hyperbolic_ii' depends on axioms: [propext, Classical.choice, Quot.sound]
'KimEff.jensen_hyperbolic_iii' depends on axioms: [propext, Classical.choice, Quot.sound]
'KimEff.hyperbolic_of_zeroFree_W' depends on axioms: [propext, Classical.choice, Quot.sound]
'KimEff.hyperbolic_of_tendsto' depends on axioms: [propext, Classical.choice, Quot.sound]
'KimEff.mem_sector_of_parabola' depends on axioms: [propext, Classical.choice, Quot.sound]
'KimEff.zerosSq_xi' depends on axioms: [propext, Classical.choice, Quot.sound]
'KimEff.halfPlane_half' depends on axioms: [propext, Classical.choice, Quot.sound]
'KimEff.Zfun_12_bounds' depends on axioms: [propext, Classical.choice, Quot.sound]
'KimEff.twelve_le_u_iff' depends on axioms: [propext, Classical.choice, Quot.sound]
'KimEff.xiGammaC_eq' depends on axioms: [propext, Classical.choice, Quot.sound]
'KimEff.xiF_eq_tsum' depends on axioms: [propext, Classical.choice, Quot.sound]
'KimEff.ThetaDeriv_eq' depends on axioms: [propext, Classical.choice, Quot.sound]
'KimEff.twoThirds_sub' depends on axioms: [propext, Classical.choice, Quot.sound]
'KimEff.jensen_hyperbolic_twoThirds' depends on axioms: [propext, Classical.choice, Quot.sound]
```

Not in the default build: Riemann.lean does not import it (add import Riemann.Deep.KimEff.Lean.Main to include it). The .lake placeholder was restored afterwards, so a re-check starts with rm .lake and lake exe cache get.

Per-file compiles: 30–250 s each, one at a time. Memory (top, mem column): 330–345 MB during elaboration, peak ~450 MB at the end of a file; the Mathlib import baseline alone is ~300 MB, so splitting files did not lower it (Core was split into Core + Step; no change).

## FORMALIZATIONS.md entry

- Explicit Kim–Lee wedge (kimeff) in Lean 4 + Mathlib, partial: for Riemann's ξ, J^{d,n} hyperbolic for n ≥ 998,616,377,771 and d ≤ D_n, d ≤ 6.45 n²/log² n and all four (c, C) pairs, from six named inputs (the zero-free region (i) = DERIVATION §4–§9 + Arb certificate, Obreschkoff–Kim–Lee Thm 3, Hadamard partial products, Gauss–Lucas for genus 0, γ > 0, Bebel's half-plane); proved: geometry, sector bound, Hurwitz for polynomials, the §10 arithmetic, n_* = ⌈Z(12)⌉, the two-thirds wedge inside for d > 1.04710e8, Θ(z²) = 8ξ(1/2+z), γ real; a = 1/2 variant (3.59) needs no half-plane input (kimEff_main, xi_jensen_hyperbolic_half; sorry-free, standard axioms). [LEAN.md](Riemann/Deep/KimEff/LEAN.md)

