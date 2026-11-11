IDEA: guard-band vertex exchange (node N1) — on an oversampled lattice P = mL Weil's lattice design admits moment-preserving band-pass weight exchanges, impossible at Chekmachev's P = 2L, that remove 91 % of the c = 47 high-frequency forced nodes; the certificate is rigorous at c = 10, but the RH-sensitive low block grows (m − 1)-fold, so the gain is a constant factor (about 6× memory at c = 47, a loss at c ≤ 16).

# ExploreDesign — optimal experimental design → Weil positivity (explore27, 11 Nov 2026)

RH relevance: N1 (Weil positivity on windows; certificate cost of the c = 47 record).

## 1. Which design tool, and why this one
| tool | what it becomes for Weil's window form Q_c | status |
|---|---|---|
| Kiefer–Wolfowitz (D), E-optimal and Elfving duality | LP/SDP duality of the minorant problem: every dual certificate of λ_min is a Krein/minorant certificate; the dual ("worst-direction measure") gives only upper bounds | known (explore19: Krein/Fejér–Riesz duals = Chekmachev's lattice minorant; R8.9, R40) |
| KW variance function at the negative nodes | Chekmachev's Woodbury step (Schur complement Id − BᵀP⁻¹B) | known |
| Dette–Studden canonical moments | Grommer/Stieltjes ladder (explore22), Verblunsky/OPUC for Li (R333) | known |
| Titterington–Silvey / Fedorov–Wynn | active-set node selection; tested as K2 thinning (sec. 4.2): fails below 0.7 × Chekmachev's rule | barrier (numbers) |
| **Böhning vertex exchange** (move weight between support points, moments fixed) | **new: needs an oversampled lattice (sec. 2)** | **tested here** |

Chekmachev's identity is a signed design: Q(f) = pole(f) + (1/2L) Σ_k μ(k)|F(πk/L)|², whose moment constraint is "Fourier transform = Weil's distribution on (−L, L)". Negative weights (forced nodes) become matrix columns at about 1000 bits. At c = 47, 8219 forced nodes lie above the low block (6047 odd + 2172 even, heights 3·10³ to 8·10⁹, from WeilZeta's full scan). They supply 6047 of the 6301 negative odd columns and, through the neighbour rule, most of the 6063-column positive block and of the 371,000-term Cauchy set. A vertex exchange would move their deficit onto neighbouring nodes with surplus.

## 2. Statements (proofs in DERIVATION.md)
- **μ at any frequency.** μ(t) = φ(t) + 2Σₙ c^{−bₙ}[bₙ cos tL − t sin tL]/(bₙ² + t²) − 2Σ_{n<c} Λ(n)n^{−1/2} cos(t log n), φ(t) = −log π + Re ψ(1/4 + it/2). At t = πk/L this is Chekmachev's μ(k); float64 code checked against mpmath at c = 47 to 10⁻¹⁵ (wsym.py).
- **Theorem G (lattice identity with a guard band).** For P ≥ 2L, t_k = 2πk/P and any real even κ ∈ L¹ supported in the guard band L ≤ |x| ≤ P/2: Q(f) = pole(f) + (1/P) Σ_k [μ(t_k) + κ̂(t_k)]|F(t_k)|² for every f ∈ L²(−L/2, L/2), in (−∞, ∞].
- **Lemma R (rigidity).** At P = 2L the only exact weight changes are a(−1)^k (the guard band is the null set {±L}). Chekmachev's weights are canonical, and every forced node must enter the matrix: no vertex-exchange direction exists. At P > 2L the directions are the band-pass sequences κ̂(t_k).
- **Theorem C (certificate on P = mL).** The m sublattices (2π/L)(Z + r/m) are orthogonal bases of L²(I). Sublattice 0 keeps exact weights (diagonal, as Chekmachev's even nodes), and the others are thresholded at τ = m·d_eff/(m − 1). Then Q − ε‖f‖² ≥ Σ_j dl_j⟨f, φ_j⟩² + (2/P)Σ_{k∈K}(w_k − τ)⟨f, cos t_k x⟩² + 2⟨f, cosh(x/2)⟩² (even f), with Chekmachev's Woodbury and Cholesky steps unchanged. m = 2 is his Theorem 5.10.
- **Patches.** κ = Σ_d Σ_i hat_i(|x|)[a cos(t_d x) + b cos(t_d x + (π/2) sgn x)] has the closed-form transform D sinc²((t − t_d)D/2)[a cos((t − t_d)x_i) + b sin((t − t_d)x_i)] plus the mirror at −t_d. Patches change only the weights; the identity stays exact, and the low block (all nodes up to Chekmachev's K2 height) keeps exact weights. The only loss is the same off-K drop that Chekmachev already makes.

## 3. Rigorous test (cert_guard.py, python-flint balls, interval Cholesky; even sector)
| c | lattice | ε certified | margin ν proven (λ_min(Id − W) estimate) | order / U | high forced nodes | time, RSS |
|---|---|---|---|---|---|---|
| 5 | P = 2L (Chekmachev) | 7.447489666e-18 | 2.0e-19 (4.00280318265e-19 = WeilZeta's margin, all 12 digits; order/U 9/94 equal) | 9 / 94 | 0 | 1 s, 0.05 GB |
| 10 | P = 2L | 9.133835606e-44 | 3.2e-45 (6.40e-45; WeilZeta 5.03e-45 with U 132) | 33 / 178 | 21 (= his 5 + 16) | 4 s, 0.11 GB |
| 10 | **P = 3L (new)** | **9.133835606e-44** | **2.16e-45 (4.33e-45)** | 54 / 305 | 29 found, 0 patchable, 27 kept | 11 s, 0.12 GB |

Every row is a rigorous certificate: Arb weights, an Arb recheck of each node flagged by the float pre-scan (flag margin 0.01, against float error < 10⁻⁸), the analytic tail beyond k_scan (Re ψ(1/4 + it/2) increases in |t|), an exact Woodbury solve, and an interval Cholesky of Id − W − νId. Theorem C therefore works, but at c = 10 it costs more than P = 2L: the forced nodes sit where the surplus is small and the dips are as wide as the patch lobe π/L (c = 10: 0 of 29 patchable). One caveat: unflagged nodes rest on the stated float bound, which is argued here but not proved. WeilZeta's RIGOR.md proves such a bound for its own scan; none is written here.

## 4. Tests for a new bound (feasibility numbers)
### 4.1 Patch census at c = 47 (patchlat.py: LP per dip cluster, lattice nodes within ±25, |h| ≤ 0.5 at the window edge; float64)
| height t | Chekmachev forced (P = 2L) | P = 3L forced → unpatched | P = 4L forced → unpatched |
|---|---|---|---|
| 3·10³–10⁴ | 273 | 561 → 232 | 752 → 440 |
| 10⁴–10⁵ | 1306 | 1519 → 232 | 2072 → 453 |
| 10⁵–10⁶ | 2178 | 1476 → 35 | 2195 → 68 |
| 10⁶–10¹⁰ | 4462 | 1729 → 0 | 3402 → 0 |
| total | 8219 | 5285 → 499 (−91 %) | 8421 → 961 (−89 %) |
All dip clusters above t = 7.2·10⁵ (P = 3L) or 8.2·10⁵ (P = 4L) are patched, with slack ≥ 0. The census is limited to windows around the dips that the P = 2L scan found; narrow dips that fall between its nodes are not counted. They lie at large t, where every dip was patchable.
At small c patching helps less: c = 16, P = 3L: 179 high forced → 128 (38 patches); P = 2L: 108 → 108 (Lemma R).

### 4.2 Can the design support be thinned? (Chekmachev's code cert_zeta.py, rigorous, c = 16 even, ε = 1.5128e-75)
| K2 | 427 (his rule 5c²/3) | 301 | 201 | 151 | 101 |
|---|---|---|---|---|---|
| margin | 7.29e-77 ✓ | 6.25e-78 ✓ | −1.75e-76 ✗ | −1.1e-60 ✗ | −9.6e-39 ✗ |
The low block, all nodes up to about 1.2c², is needed: a Fedorov–Wynn-type thinning gains at most 30 %. This block holds the near-null ladder, which is what makes the certificate RH-sensitive.

### 4.3 Cost estimate at c = 47 (P = 3L, from 4.1 and 4.2)
- Low block: 2 × 1842 rank-one columns, against Chekmachev's 1842 odd ones.
- High block: about 500 unpatched nodes, against his 8219.
- Matrix: order ≈ 950, U ≈ 4500, so about 3·10⁷ balls ≈ 6 GB at 1000 bits, against his 1.75·10⁸ ≈ 35 GB. The Cauchy set |G| falls from 371,000 to about 4·10⁴.
- Verdict: a constant factor, roughly 6× in memory and 10× in Cauchy sums, at c = 47. The asymptotics are no better. The low block grows like c², and the unpatched band [t_K2, 10⁶] grows with the prime sum 2Σ_{n<c}Λ(n)/√n. At c ≤ 16 the method is worse than P = 2L.
- The scan cost is unchanged: k_scan grows like e^{2Σ_{n<c}Λ(n)/√n}, which is the real wall beyond c ≈ 60 (weildefl).

## 5. Adversarial test (adversarial.py; float64 Galerkin, J = 40, odd nodes ≤ 4001)
γ₁, γ₂ (14.13, 21.02) are merged into an off-line pair 17.58 ± iδ, keeping the zero count (ξ·(s − ρ₊)(s − ρ₋)/((s − ρ₁)(s − ρ₂)) and conjugates).
| c | ζ even / odd λ_min | moved pair δ = 0.3 even / odd | δ = 0.05 even / odd |
|---|---|---|---|
| 2 | 1.3e-3 / 7.7e-2 | 1.9e-5 / 7.2e-3 | 6.4e-5 / 8.0e-3 |
| 3 | 5.7e-8 / 1.7e-5 | −2.1e-3 / −4.4e-3 | −1.5e-5 / −3.0e-6 |
| 5 | (below float resolution; certified 7.4e-18 / 5.7e-15) | −4.2e-2 / −4.9e-2 | −1.1e-3 / −1.2e-3 |
| 10 | (certified 9.1e-44) | −0.16 / −0.18 | −4.1e-3 / −4.8e-3 |
- From c = 3 on, the moved-pair form has a negative direction, so every sound certificate fails for it, the guard-band one included. At c = 2 it is still positive: small windows do not see the pair (explore19: c* ≈ 0.072m + 0.73 = 2.0).
- The moved pair changes μ by at most 3.5·10⁻⁴ on t ∈ [3000, 3100] at c = 47, so the patches are identical for ζ and for the moved pair. The new step is RH-blind by design. All RH content stays in the exact low block, as in Chekmachev's certificate.

## 6. Honest assessment
- New: Theorem G (lattice identity with a guard band), Lemma R (rigidity of Chekmachev's weights at P = 2L), and Theorem C (an m-sublattice certificate), checked rigorously at c = 5 (m = 2, equal to WeilZeta to 12 digits) and c = 10 (m = 2, 3). These are tools for certificate engineering, not mathematics about ζ.
- No new bound: at c = 10 the same ε is certified with a larger matrix. The c = 47 gain (about 6×) is an estimate and was not certified (memory limit 0.3 GB). Asymptotically nothing changes: the low block and the scan dominate.
- Barrier read-off: design-theoretic dualities add nothing beyond the minorant LP. The one genuinely new freedom, vertex exchange, exists only with oversampling and pays an (m − 1)-fold low block for it.

## Suggested NOGO line
- explore27 (Orca, exploration: optimal experimental design — Kiefer–Wolfowitz/E-optimal/Elfving duality, canonical moments, Titterington–Silvey and Böhning vertex exchange → N1) IDEA (engineering, constant factor): the dualities are the minorant LP; Chekmachev's P = 2L weights are rigid (exact changes only a(−1)^k), but on P = mL a guard band L < |x| < P/2 allows band-pass patches (Theorem G, m-sublattice Theorem C; rigorous at c = 5 = WeilZeta to 12 digits and c = 10 with P = 3L) that remove 91 % of c = 47's 8219 high forced nodes (all above t = 7·10⁵); K2 cannot be thinned below 0.7 × 5c²/3 (c = 16: 301 ✓, 201 ✗), so the low block grows (m − 1)-fold: about 6× less memory at c = 47 (estimate), worse at c ≤ 16, scan wall unchanged [ExploreDesign](Riemann/Deep/ExploreDesign/REPORT.md)

## Files
- wsym.py: μ(t) at any frequency (float64), validated against mpmath.
- DERIVATION.md: μ(t), Theorem G, Lemma R, Theorem C, patches and their loss.
- cert_guard.py: rigorous certificate on P = mL with patches (c m eps prec K2; CG_SCAN_ONLY=1 for the scan and patch stage only).
- patchlat.py: c = 47 patch census; census47.py: height census of the forced nodes.
- adversarial.py: the moved-pair test.
- logs/: cert_guard.log, census47_P3L.txt, census47_P4L.txt, e1_K2_c16_even.log, scan_patch_c16.log, adversarial.log.
- Replay: `uv run --no-project --with numpy --with scipy --with python-flint --with mpmath python3 cert_guard.py 10 3 9.133835606e-44 384 201` (11 s, 0.12 GB). The census needs /tmp/e27_forced47.npy from census47.py, then `patchlat.py 3 0 9000 12 0.5 25 0.5` (about 1 min).
