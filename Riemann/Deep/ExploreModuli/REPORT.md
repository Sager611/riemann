NO IDEA: the near-optimal spectral-gap proofs get zero-wise positivity of the trace formula for free from self-adjointness (spectrum on ℝ ∪ iℝ) and spend all their randomness on the averaged geodesic side, whereas for Dirichlet characters mod q the averaged prime side is already free where a positivity argument needs it (orthogonality and the sign of the prime term up to length 2 log q), zero-wise positivity of the same kernels for all supports is GRH up to real zeros, and kernels positive on the whole strip resolve no finer than sech(πx) and amplify an off-line zero by at most tan(π/4 + πα/2), so no random-model input transfers (N13).

# Exploration explore10: moduli spaces, Weil–Petersson random surfaces, spectral gaps → RH tree (2026-11-05)

Field: spectral gaps of random hyperbolic surfaces (Weil–Petersson and random covers) and their arithmetic twins. Brief: `.scratch/orca/explore10-brief.md`. Proofs and formulas: [derivations.md](derivations.md). Notation: a zero ρ = 1/2 + α + ix enters the explicit formula as H(x − iα), H(z) = ∫G(u)e^{izu}du, G supported in [−R, R]; R = A log q (or A log g); DH = Davenport–Heilbronn.

## Survey (state of the art, Nov 2026)

- Weil–Petersson (WP): Mirzakhani 2013, λ₁ ≥ 0.0024 with high probability (Cheeger); Wu–Xue 2102.05581 and Lipnowski–Wright 2103.07496, 3/16 − ε (trace method with simple geodesics, R = 4 log g); Anantharaman–Monk 2304.02678, 2502.12268 (Friedman–Ramanujan (FR) functions), road map 2403.12576, intermediate 2/9 − ε 2604.09792, notes 2601.15157: λ₁ ≥ 1/4 − ε with probability → 1; tangle removal by Moebius inversion 2401.01601; Hide–Macera–Thomas 2508.14874: 1/4 − O(g^{−c}) by the polynomial method; cusps: 2602.08352, and 2611.01477 (1 Nov 2026: gap → 1/4 iff n = o(√g); expected number of eigenvalues in (0, 1/4 − ε] ≍ (n+1)²/g, sharp); uniform volume asymptotics 2610.35492.
- Random covers: Magee–Naud–Puder 2003.10911 (3/16 − ε); Hide–Magee 2107.05292 (cusped, 1/4 − ε, via Bordenave–Collins strong convergence; closed surfaces with λ₁ → 1/4 exist); Chen–Garza-Vargas–Tropp–van Handel 2405.16026 (polynomial method: E tr P(X_N) is rational in 1/N, plus Markov's polynomial inequality); Magee–Puder–van Handel 2504.08988 (closed surfaces, beyond rational settings); rate n^{−b}: 2505.08479; random Belyi (Brooks–Makover) surfaces 1/4 − c/log n: 2511.02517; surveys 2503.21619, Petri's Bourbaki talk 2607.06331.
- Arithmetic surfaces: Selberg 3/16 (Kloosterman sums + Weil's bound); Kim–Sarnak 975/4096 (Sym⁴ functoriality); Sarnak–Xue: N(λ_s, Γ(q)) ≪ D_q^{2(1−s)+ε} for cocompact arithmetic Γ (density hypothesis; as quoted in 2608.18236, theorem thm:SX-Oh); thin groups with power saving 2608.18236; GL₂ over number fields 2402.11761. Euler products of Selberg zeta in the strip 1809.10140 (N12 aside, no zero information).

What the proofs use (Petri 2607.06331, AM 2403.12576; details in derivations D1):

| step | surfaces | source | arithmetic counterpart (primitive χ mod prime q) | status |
|---|---|---|---|---|
| S1 positivity | ĥ_R ≥ 0 on ℝ ∪ iℝ suffices: spectrum r_j ∈ ℝ ∪ i[−1/2, 1/2] | self-adjointness, no randomness | zeros anywhere in the strip: zero-wise positivity for all R forces every non-real zero onto the line (GRH up to real zeros) | missing |
| S2 Markov | P(λ₁ ≤ 1/4 − α²) ≤ E[trace]/ĥ_R(iα) | — | same | free |
| S3 identity term | ≈ g/R² per surface, so R ≥ (1/α) log g | Plancherel | ≈ log q per member, so R ≥ (1/α) log log q | free |
| S4 trivial eigenvalue | λ₀ = 0 gives e^{R/2}, cancelled by the simple-geodesic main term | WP average (Mirzakhani) | the missing χ₀ gives −S₀ ≈ −(4/R)e^{R/2}, cancelled by the AP primes (q−1)S₁ | provable up to R = 2 log q (sign) |
| S5 D = 1/4 − d²/dℓ² | on the test function: kills λ₀, and (after the weight ℓ/(2 sinh(ℓ/2)), D^m with m > deg p + 1) the p(ℓ)e^ℓ (FR) terms | algebra | (1/4 + r²)^m kills the pole terms | free, not needed |
| S6 lengths A log g, every A | 1/g-expansions of V_g^T(ℓ) with FR coefficients, or smoothness in 1/g (polynomial method) | randomness | square-root precision of AP prime sums at x = q^A, A > 2 = GRH on average for the same family | circular beyond support 2 |
| S7 tangles | rare small subsurfaces with short boundary, probability → 0 polynomially in g | randomness | none: orthogonality is exact below length log q (S₁ ≡ 0) | — |
| S8 strong convergence | random permutation reps → regular rep of π₁ (non-amenable, 1/4 = λ₀(ℍ)) | randomness | direct sum of characters mod q → regular rep of ℤ^{π(y)} (amenable): limit norm = trivial value | dead |

Contrast: on WP surfaces the exceptional count is tangle-limited (≍ (n+1)²/g for fixed ε with n ≥ 1 cusps, sharp, 2611.01477; no α-power appears); in the character family there are no tangles and any α-dependence (q^{−2α}) comes from the length cutoff R = 2 log q.

## Transfers, tests, verdicts

### T1. Trace-method first moment over characters mod q (best transfer)

E_χ Σ_ρ H(γ_ρ) = Z = arch − (2/(q−2))·D, D = (q−1)S₁ − S₀, with S₁ the prime sum over n ≡ 1 (mod q) and S₀ over q ∤ n, G = (1 − |u|/R)₊ (H ≥ 0 on ℝ ∪ iℝ, the AM class). Unconditional: S₁ ≥ 0, so Z ≤ Z_up = arch + 2S₀/(q−2). Exact sieve to 1.1·10⁹ (`t1_moment.py`, 71 s, 154 MB):

| q | A | (q−1)S₁/S₀ | D/√q | Z/arch (true) | Z_up/arch (provable) |
|---|---|---|---|---|---|
| 101 | 1 | 0 | −0.52 | 1.043 | 1.043 |
| 101 | 2 | 0.627 | −1.50 | 1.180 | 1.483 |
| 101 | 3 | 0.934 | −1.89 | 1.282 | 5.29 |
| 101 | 4.5 | 0.996 | −2.28 | 1.402 | 109.8 |
| 1009 | 2 | 0.786 | −1.95 | 1.034 | 1.157 |
| 1009 | 3 | 0.980 | −3.96 | 1.073 | 4.56 |
| 10007 | 2.25 | 0.886 | −6.97 | 1.024 | 1.213 |
| 100003 | 1.8 | 0.750 | −4.82 | 1.004 | 1.015 |

The analogue of "simple geodesics cancel λ₀" is visible: (q−1)S₁ → S₀ with residual |D| ≈ 0.4–7 √q (the square-root, Friedman–Ramanujan scale). Provably, Z stays near arch only up to A = 2; beyond it the needed precision is GRH on average for the same family (the explicit formula is an identity). This gap does not matter: if every single zero contributed ≥ 0 (true when all zeros lie on the line or on the real segment, as on a surface), Markov at the provable A = 2 already gives P(zero in [1/2 + α, 1] × [−π/(3R), π/(3R)]) ≤ Z_up/H(iα) ≲ 2α²(log q)² q^{−2α} → 0 for every α > 0 (q = 10007: 0.31, 0.26, 0.073 for α = 0.05, 0.1, 0.25; `t1b_out.txt`). With positivity no long lengths are needed; without it no length helps.
- Randomness test: randomness pays for S4, S6, S7. Arithmetically character orthogonality replaces it, exactly up to length log q; up to 2 log q the sign of the prime term does the rest; beyond that, Brun–Titchmarsh, Bombieri–Vinogradov and Montgomery–Hooley stay inside the classical large-sieve / 1-level-density range (support < 2). Not new.
- DH test: DH-type combinations ω̄Λ(s,χ) + ωΛ(s,χ̄) satisfy the functional equation but have no Euler product, so S₀, S₁ do not exist for them; the Euler product enters only S4/S6, which is not the bottleneck, and the decisive step (positivity) uses nothing DH lacks and fails for DH (T2). FAILS: no RH content.
- Verdict: dead (classical up to support 2, circular beyond it; the positivity step, required for all R, is GRH up to real zeros).

### T2. Positivity versus amplification (why S1 cannot be replaced)

To drop unknown zeros pointwise without GRH, a kernel needs Re H ≥ 0 on the whole closed strip. Exactly the Poitou–Odlyzko kernels G = k/cosh(u/2), k positive definite, have this (Bochner on the boundary line; derivations D3). Two consequences hold for every such H, whatever its support R:
  (resolution) on the line Re H = k̂ * sech(π·), so no window is narrower than sech(πx) (measured half-width 0.43 at R = 160 against 0.419; `t2c_out.txt`);
  (Harnack cap; strip → disc by tanh(πz/2)) Re H(x − iα) ≤ tan(π/4 + πα/2)·Re H(x) for all x and 0 ≤ α < 1/2, sharp for positive harmonic functions.
An off-line zero thus counts as at most 2.414 (α = 1/4) or 3.224 (α = 0.3085) on-line zeros at the same height. For characters mod q the first moment is unconditionally Z ≥ k(0)(log q − 5.4) − O(k(0)R/q), so the Markov bound is ≥ (log q − 5.4 − O(R/q))/(4π tan(π/4 + πα/2)) for every support R ≤ q (main term 0.46 at q = 10⁶, 2.28 at q = 10³⁰ for α = 1/4): it grows with q and never isolates a member. Near Re s = 1 (α = 1/2 − δ) the cap is ≈ 2/(πδ), the classical pole-anchored mechanism of width 1/log q (R30). Numbers (`t2_out.txt`):
- Line-positive kernel (surface type): amplification H(iα)/H(0) = 1.64, 5.86, 220 at R = 10, 20, 40 (α = 1/4), but the worst complex off-line pair weighs −0.65, −6.6, −384 times H(0); at α = 0.3085, R = 40: +1502 against −2756. One complex pair outweighs a detected real zero.
- Full-strip kernel: Re H ≥ 0 on the strip grid (min ≥ 2.8·10⁻⁹); amplification at α = 0.3085 is 1.13, 1.31, 1.51, 1.64 for R = 5, 10, 20, 40 (cap 3.22); at α = 0.45: 1.30 → 3.86 (cap 12.71).
- DH zero 0.808517182457 + 85.699348485378i (`t2b_out.txt`): line-positive pair weight/H(0) changes sign with R and is −1204 at R = 80; full-strip ratio 0.92, 1.10, 0.99, 1.00, 1.00 (R = 5 … 80): no weight.
- Randomness test: this step involves no randomness. Only self-adjointness (a Hilbert–Pólya operator, R51) supplies positivity on ℝ ∪ iℝ.
- DH test: the cap holds for DH too (no Euler product is used). The DH zero is correctly not counted, and nothing excludes it either. PASSES vacuously; it is a barrier, not a mechanism.
- Verdict: barrier (routine: Bochner + Harnack; not novelty-checked; possible COROLLARIES line).

### T3. Strong convergence ↔ amenable limit

Dirichlet characters are abelian covers of Spec ℤ. For fixed y the vectors (χ(p))_{p≤y} equidistribute on the torus T^{π(y)}, so the direct sum of all characters mod q converges strongly to the regular representation of ℤ^{π(y)}; its limit norm is the trivial value because ℤ^{π(y)} is amenable (Brooks: amenable covers have no uniform gap). Strong convergence holds; a gap needs a non-amenable limit. Random permutation representations of the free group converge to Kesten's non-amenable value. Numbers (`t3_out.txt`): fixed y = 30 (k = 10 primes), max over χ ≠ χ₀ of 2Re Σχ(p)/(2k) = 0.553, 0.711, 0.848, 0.893 for q = 10³ … 10⁶ (→ 1). Random 2k-regular Schreier graphs of S_n: 0.4308, 0.4353, 0.4359 for n = 200, 800, 2000 (Kesten 0.4359). For y = q^θ the family maximum tracks the Steinhaus union bound √(log q/k) (69–98% of it) and stays 2.0–2.4 times above Kesten (union-bound ratio √(log q/2) = 2.6 at q = 10⁶; q = 10⁶ + 3, θ = 1: 0.0103 against 0.0133 and 0.0050). Heuristically this is a finite-family large deviation: classical large-sieve mean values.
- Randomness test: random permutations ↔ equidistribution of characters (orthogonality); the replacement is the large sieve. Not new.
- DH test: DH coefficients are not multiplicative, so the torus model does not apply; but the conclusion (the trivial bound) carries no zero information for anything. No RH content.
- Verdict: dead.

### Ranking (RH relevance × chance × gain)

| transfer | node | chance | gain | verdict |
|---|---|---|---|---|
| T1 family first moment | N11 | ≈ 0: needs zero-wise positivity (GRH up to real zeros) | "almost all" near the real axis (q^{−2α}), only if positivity | dead |
| T2 positivity vs amplification | N13, N1 | theorem | barrier only | barrier |
| T3 strong convergence | N12 | 0 | none | dead |

### (iii) Side result

Not attempted. There is no cheap target: the "large genus" thresholds (g_ε, C_ε in 2611.01477; Mirzakhani's 0.0024) rest on Weil–Petersson volume estimates with unspecified constants (2610.35492 gives uniform asymptotics with an unspecified O(1/g)).

## Suggested NOGO line

- explore10 (Orca exploration: random hyperbolic surfaces, Weil–Petersson volumes, random covers, strong convergence → N1, N11, N12, N13) NO IDEA: the trace method's two inputs swap roles. On surfaces positivity is free (self-adjoint, r ∈ ℝ ∪ iℝ) and randomness pays for the averaged geodesic side (Mirzakhani integration, Friedman–Ramanujan cancellation, tangles). For characters mod q the averaged prime side is free up to length 2 log q (orthogonality and sign; with positivity this already gives P(zero near 1/2 + α) ≲ q^{−2α}), and zero-wise positivity of those kernels for all supports is GRH up to real zeros. Strip-positive kernels (exactly k/cosh(u/2), k positive definite) resolve no finer than sech(πx) and amplify an off-line zero by at most tan(π/4 + πα/2) (Harnack; 1.64 against the cap 3.22 at R = 40), so a family first moment never isolates a member (Markov bound ≥ (log q − O(1))/(4π tan(π/4 + πα/2)) for supports R ≤ q); line-positive ones lose to one complex pair (−2756 against +1502). Strong convergence needs a non-amenable limit, and characters converge to the amenable torus (0.55 → 0.89 for q = 10³ → 10⁶ against Kesten's 0.436). Filter: positivity that comes from self-adjointness transfers only with a Hilbert–Pólya operator. [ExploreModuli](Riemann/Deep/ExploreModuli/REPORT.md)

## Files

- `REPORT.md`, `derivations.md` (D1 trace method; D2 family moment and the Markov step; D3 strip-positive kernels: characterisation, resolution, Harnack cap; D4 amenable limit).
- `t1_moment.py` → `t1_out.txt` (segmented sieve to 1.1·10⁹, primes and prime powers; 71 s, 154 MB); `t1b_markov.py` → `t1b_out.txt` (reads `t1_rows.json`, which `t1_moment.py` writes).
- `t2_harnack.py` → `t2_out.txt`; `t2b_dh.py` → `t2b_out.txt`; `t2c_resolution.py` → `t2c_out.txt`.
- `t3_amenable.py` → `t3_out.txt` (FFT over the character group; random Schreier graphs).
- Run: `OMP_NUM_THREADS=1 uv run --no-project --with numpy --with mpmath python3 t1_moment.py 1.1e9`; the others take < 12 s.

RH relevance: N13 (barrier: positivity comes from self-adjointness), with consequences for N11 (families), N12 (Selberg zeta analogue) and N1 (positivity).
