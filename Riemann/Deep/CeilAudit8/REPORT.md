ROBUST: explore49's fibre obstruction holds for even q as well (three parity classes with the same conic structure, ≤ 2 pairs per generic conic, σ = 2 shape, mixed-parity pairs never needed); explore29's quantizer barrier holds, though Fleischer's hypothesis family is the Laguerre inequality of every derivative Ξ^{(m)} (reach ≥ 3 ḡ against 1.1 ḡ for Ξ itself) and Lloyd–Max uniqueness never sees it; explore39's Poincaré margin is intrinsic (P − R ≤ min(λ₀/W, λ₁/W) for every Barta-admissible ψ), and Barta's λ₀ accuracy binds first.

# ceilaudit8: falsification of our barriers explore49 (N9), explore29 (N8, N6, N5), explore39 (N1, N13); explore45, explore42 skipped (Orca, 20 Nov 2026, 13:06–13:55 UTC)

## B1. explore49 (N9): "the determinant method has no grip on 13/84", even q
**Step fixed.** Bourgain's second-spacing count B1 (1408.5794 §3), pairs of minor arcs a/q with |x_j − x_j'| < 1/(2Y_j); the barrier: generic conics B_δ = n hold ≤ 2 counted pairs, so fibrewise determinant counts give σ = 0. **Class measured:** odd q only (spacing2.c l.8; A = 4a needs gcd(4a, q) = 1).

**Even-q coordinates (derived; checked numerically, gauss_check.py: 87,804 cases, 0 failures).** S_I = (1/q)Σ_k G(a, k; q)·(stationary-phase integral in √(KP) − k), P = aq, K = t/π.
- q ≡ 2 (mod 4), q₁ = q/2: G(a, k; q) = G(2a, k; q₁)(1 + (−1)^{k+1}) vanishes for k even; for k odd G = 2e(−\overline{8a}k²/q₁)G(2a, 0; q₁), |G| = √(2q).
- q ≡ 0 (mod 4): G vanishes for k odd; G(a, 2k′; q) = e(−ā k′²/q)G(a, 0; q), |G| = √(2q).
- So the dual sum runs over one parity class: k = b − 2h′ with b the nearest integer of the allowed parity to √(KP) (η = √(KP) − b ∈ [−1, 1)), h′ ≤ H′ = H/2. In the h′ basis:
  - q ≡ 0 (4): x₂ = −ā/q, x₁ = −b·x₂ (mod 1);
  - q ≡ 2 (4): x₂ = −\overline{2a}/q₁ (inverse mod q₁), x₁ = −b·x₂ (odd q: x₁ = −2b·x₂ in the same true-phase basis);
  - spacing2.c stores −x₁ (its x₁ = 2b·x₂), and spacing3.c does the same for every class: a global sign per coordinate leaves every closeness count unchanged;
  - both: x₃ = 2^{3/2}(2/3)/D, x₄ = √2·η/D (D as for odd q), windows from H′. Equivalently, in the h basis X₃, X₄ are as for odd q and X₁ (mod ½), X₂ (mod ¼) are multivalued.
- Conic identity, same SL₂ algebra with (A, modulus) = (a, q) for q ≡ 0 (4), (2a, q/2) for q ≡ 2 (4): A·mod − A′·mod′ = P − P′ = −B_δ(A, mod) (odd q: 4(P − P′)). Watt's Weyl-differenced form (Bourgain–Watt 1505.04161 TeX l.2259) has x = ({ā/q}, {āb/q}, 1/√(μq³)) for every q: the same (A, mod) = (a, q) structure.
- **Mixed-parity pairs are not forced.** Split the arcs into the three classes and run Bourgain's chain on each: |S| ≤ Σ_class |S_class|, each class obeys the same inequality with its own B1, a constant factor ≤ 3. Counted inside one h-basis sieve (each even arc = two points X, X + (½, 0, 0, 0)), mixed pairs depend on the arbitrary representative of X₂ mod ¼ (measured 0.61× the random model, t = 10¹⁶, 10¹⁸): an artefact.

**Test** (spacing3.c = spacing2.c + class selector; class odd reproduces spacing2.c's output line for line (apart from the added class header) at t = 10¹², 10¹⁶; every run ≤ 87 MB, ≤ 32 s).
| log10 t, α, Q | class | arcs | (x2,x3) pairs | max/generic conic | max/δ | identity fails | true/rand x2x3 | all4/x2x3 · H² |
|---|---|---|---|---|---|---|---|---|
| 16, ½, R | odd / 2(4) / 0(4) | 6027 / 1517 / 1579 | 85250 / 20750 / 23708 | 1 / 1 / 1 | 6 / 3 / 3 | 0 | 1.21 / 1.20 / 1.24 | 0.062 / 0.070 / 0.070 |
| 18, ½, R | same | 24310 / 6077 / 6114 | 555684 / 138352 / 141334 | 2 / 2 / 1 | 14 / 4 / 4 | 0 | 1.23 / 1.24 / 1.22 | 0.064 / 0.065 / 0.062 |
| 20, ½, R | 2(4) / 0(4) | 24248 / 24134 | 881998 / 869714 | 1 / 2 | 6 / 8 | 0 | 1.24 / 1.24 | 0.064 / 0.068 |
| 22, ½, R | odd / 2(4) / 0(4) | 386246 / 96068 / 96631 | 22.2M / 5.52M / 5.57M | 1 / 2 / 1 | 6 / 5 / 4 | 0 | 1.24 / 1.24 / 1.24 | 0.066 / 0.067 / 0.067 |
| 20, 0.45, R | odd / 2(4) / 0(4) | 98087 / 23840 / 25303 | 3.52M / 0.85M / 0.91M | 2 / 1 / 2 | 7 / 7 / 10 | 0 | 1.21 / 1.22 / 1.20 | 0.066 / 0.065 / 0.064 |
| 20, 0.42, R | odd / 2(4) / 0(4) | 90717 / 27297 / 23668 | 5.55M / 1.43M / 1.28M | 1 / 1 / 1 | 2 / 3 / 2 | 0 | 1.28 / 1.22 / 1.23 | 0.064 / 0.056 / 0.059 |
| 18, ½, 4R | odd / 2(4) / 0(4) | 1165 / 289 / 268 | 1382 / 306 / 298 | 1 / 1 / 1 | 1 | 0 | 1.23 / 1.55 / 1.45 | 0.93 / 1.05 / 1.07 |
| 20, ½, R, V-route | 2(4) / 0(4) | 24248 / 24134 | 221814 / 216492 | 1 / 2 | 6 / 8 | 0 | 1.24 / 1.23 | 0.066 / 0.067 |
(Huxley–Watt shape per class, normalized by the class's arc share f² and its x₂ window (4w₂ in the h′ basis; the x₃ condition is unchanged): 0.0348 / 0.0335 / 0.0353 (t = 10¹⁶), 0.0358 / 0.0357 / 0.0360 (10¹⁸), 0.0361 / 0.0362 / 0.0361 (10²²) for odd / 2(4) / 0(4): the same constant. all4/x2x3·H² in the full-h normalization H; constant across t = σ = 2 shape. Degenerate share: g = 0 translations 1.0 % (α = 0.45) and 17–19 % (α = 0.42) in every class; square-discriminant conics ≤ 1 %.)

**Verdict B1: ROBUST.** Even q are a third of the arcs and the same obstruction: within each parity class the (x2, x3) pairs are integer points on conics B_δ(A, mod) = P′ − P, generic conics hold ≤ 2 counted pairs, Huxley–Watt's (x2, x3) shape is the random-model shape (true/rand 1.20–1.24), and the full count has the σ = 2 shape (all4/x2x3·H² = 0.056–0.070). Exact class: fibrewise (conic, δ-strip, hyperbola) determinant counts of Bourgain's second-spacing pairs, every q, split by parity of q; not covered: cross-fibre equidistribution of √(KP) (explore49 §4.6).
Reviewer (claims 1–4 pasted inline: the two Gauss-sum parity formulas, the h′ coordinates with the conic identity, the class split): 1, 2, 4 CORRECT; 3 GAP: x₁′ = −b·x₂′, not b·x₂′, in the stated phase basis (fixed above; the counts use one sign for all arcs, so no number changes). Logs: logs/b1_*.txt.

## B2. explore29 (N8, N6, N5): "the RH-sensitive part of quantizer conditions is Fleischer's log-concavity = the first Laguerre inequality"
**Step fixed.** Fleischer's lemma (log-concave density ⇒ unique Lloyd–Max quantizer), applied to p = |Ξ|^α between consecutive zeros; its hypothesis is L₁(Ξ) = 2(Ξ′² − ΞΞ″) ≥ 0 (exact: (log|Ξ|)″ = −L₁/(2Ξ²) off the zeros). **Class:** p built from Ξ itself, k ≤ 4, α = 1, 2.

**Freedoms the lemma admits, excluded by the class.**
- Fleischer applies to every Laguerre–Pólya function made from Ξ. Under RH every derivative Ξ^{(m)} is LP, so |Ξ^{(m)}| is log-concave between its zeros: hypothesis L₁(Ξ^{(m)}) ≥ 0 (the x-dependent Laguerre inequality of the derivative; the coefficient Turán inequalities at the origin, a different member of this circle, hold unconditionally, Csordas–Norfolk–Varga 1986).
- Higher Laguerre inequalities: |f(x + iy)|² = Σ_n L_n(f)(x) y^{2n}/(2n)!, L_n = Σ_k (−1)^{k+n} C(2n, k) f^{(k)} f^{(2n−k)}; for real zeros |f(x+iy)|²/f(x)² = ∏(1 + y²/(x − γ)²), an off-line pair g₀ ± iδ contributes 1 + 2y²(w² − δ²)/A² + y⁴/A² (w = x − g₀, A = w² + δ²). So L_n(x)/f(x)² = (2n)! [Y^n] of that product: exact, every order, no differentiation.
- Heat flow forward and other multiplier sequences: they move off-line zeros toward the line and merge them (expected to lose reach; not tested).

**Test** (explore29's model: a ζ pair replaced by g₀ ± iδ; Lehmer #6709/6710, generic #6607/6608; 22,491 zeros to 19,999 + a smooth density tail; b2_laguerre.py, b2c_turan_stable.py, b2d_quant.py; peaks ≤ 180 MB by /usr/bin/time; two first drafts were killed by memcap at 2.0 GB and 1.1 GB before chunking).
- Cross-check: max (log|Ξ_δ|)″ḡ² on the merged interval = +48.045, +6.048, +1.615, +0.083 at δ/ḡ = 0.2, 0.5, 0.75, 1.0 (Lehmer) and +0.311 at 1.0 (generic), explore29's table to every digit; Σ_n c_n(g₀)δ^{2n} = 0 to 1e−13 (f(g₀ + iδ) = 0).
- **L₁(Ξ) is blind beyond δ ≈ 1.1 ḡ** (Lehmer: fails at 1.0, holds at 1.1; generic: fails at 1.0, holds at 1.25).
- First failing Laguerre order n* (Lehmer / generic): δ/ḡ = 1.1: 2; 1.25: 2 / 2; 1.5: 3 / 2; 2: 4 / 4; 3: 6 / 7; 4: 9 / 9; 6: 15 / 15; failures form a band of 3–6 orders, centred near πδ/ḡ (real zeros on both sides give [sinh(πy/ḡ)/(πy/ḡ)]², whose coefficient sequence the pair factor (1 − y²/δ²)² can make negative only near 2n ≈ 2πδ/ḡ).
- L₁(Ξ^{(m)}) (stable Taylor form; control with the true zeros: every m ≤ 7 positive, min +0.37 … +0.88): fails at m = 2 (δ = 1.25 ḡ), 2–3 (1.5), 4–5 (2.0), 5–6 (2.5), 7 (3.0), generic alike (also m = 1 at 1.0–1.25): **Fleischer's hypothesis on derivatives reaches δ ≥ 3 ḡ, m ≈ 2.4 δ/ḡ**; the window in m is 1–2 wide (the pair is pulled to the axis, then becomes real).
- **The quantizer's own conclusion sees none of it**: 2-point Lloyd–Max stationary points of p = |Ξ^{(m)}|^α on the dip interval and on 66 sub-intervals, α = 1, 2: unique in all 12 (δ, m) cases, including every case where L₁(Ξ^{(m)}) fails.
- DH (9 off-line pairs below 420, δ/ḡ = 0.06–0.31): L₁ fails at all 9 (min c₁ḡ² from −18.8 to −596).

**Verdict B2: ROBUST, class stated too narrowly.** Quantizer conditions add no RH information beyond N6's inequalities: Fleischer's hypothesis family is {L₁(Ξ^{(m)}) ≥ 0, m ≥ 0} (reach δ ≈ 0.42 m ḡ), the higher Laguerre family {L_n(Ξ) ≥ 0} reaches δ ≈ 0.4 n ḡ, but uniqueness, the only thing a quantizer observes, stayed unique in every case beyond δ = 0.35 ḡ. Since ζ's off-line zeros have δ < ½ and ḡ = 2π/log(t/2π), L₁(Ξ) alone covers δ = ½ only for log(t/2π) < 4π·1.1, t ≲ 6·10⁶ (in this pair model); above, the first failing order grows like 2.5δ/ḡ = log(t/2π)/5 at δ = ½ (band centre log(t/2π)/4), the derivative order like m ≈ 2.4δ/ḡ (2–3 at 1.25–1.5 ḡ).

## B3. explore39 (N1, N13): "beyond the Markov range a certificate needs a Poincaré margin of 1.2e−4 (c = 3), 1.4e−8 (c = 4)"
**Step fixed.** Identity C: for even ψ > 0, f = ψg, Q(f) − ∫V f² = E_ψ(g) − R·Var_ν(g), V = Qψ/ψ, R = 2⟨C, ψ⟩², ν ∝ Cψ, C = cosh(x/2); certificate = Barta (V ≥ 0 pointwise) + Poincaré (P(E_ψ, ν) ≥ R). **Class measured:** ψ = ψ₀ (ground state of Q^ev).

**Freedoms.** (a) any other positive ψ (another ground-state transform, e.g. the Perron vector of the Markov part M = Q^ev − 2CCᵀ, or the pole vector C); (b) the conservative-process certificate of M's h-transform; (c) a non-Markov comparison form in place of E_ψ.

**Lemma (new here, three lines).** If V ≥ 0, then for every f not proportional to ψ, with g = f/ψ, P − R ≤ [E_ψ(g) − R Var_ν(g)]/Var_ν(g) = [Q(f) − ∫V f²]/W_ψ(f) ≤ Q(f)/W_ψ(f). With f = ψ₀ (when ψ is not ∝ ψ₀) and f = φ₁: **P − R ≤ min(λ₀/W_ψ(ψ₀), λ₁/W_ψ(φ₁))** for every admissible ψ. The same step bounds any splitting Q(f) = ∫V f² + B(f) with V ≥ 0: B(φ₁) ≤ λ₁, so a comparison partner B ≥ κB₀ needs B₀(φ₁) ≤ λ₁/κ, i.e. it must reproduce Q's near-kernel.

**Test** (b3_transform.py on explore39's Legendre Galerkin, n = 20; 11 s, 59 MB):
| c | λ₀(Q^ev) | λ₁(Q^ev) | μ₀(M), μ₁(M) | M-Perron margin μ₁/\|μ₀\| | ψ₀: (P−R)/R, cap/R | Perron(M): inf V | C: inf V | ψ₀ + sφ₁ |
|---|---|---|---|---|---|---|---|---|
| 2.5 | 1.03e−5 | 0.0766 | −1.813, 0.0531 | 2.9e−2 | 6.72e−3, 7.61e−3 | −0.108 | −0.334 | Barta fails for s ≥ 1e−4 |
| 3 | 5.58e−8 | 1.46e−3 | −2.170, 9.49e−4 | 4.4e−4 | 1.18e−4, 1.20e−4 | −0.150 | −0.335 | fails for s ≥ 1e−6 |
| 3.5 | 2.51e−10 | 1.88e−5 | −2.544, 1.21e−5 | 4.7e−6 | 1.56e−6, 1.57e−6 | −0.133 | −0.283 | fails for s ≥ 1e−4 |
Barta holds in 5 of the 22 computed (ψ, c) cases (ψ₀, and ψ₀ + 1e−6φ₁ at c = 2.5, 3.5); the cap held in all 5, and in 16 of the other 17; the one violation (c = 3.5, s = 1e−2) has inf V < 0, outside the lemma. The Perron(M) and C rows have (P − R)/R < 0: as certificates they fail outright.

**Verdict B3: ROBUST; the margin is intrinsic to Q^ev, and the binding requirement is finer still.** Every transform's margin is capped by Q's own second even eigenvalue (λ₁/(R·W), and λ₀/(R·W) once ψ ≠ ψ₀); M's Perron vector fails Barta by 0.11–0.15 and its own conservative-process margin is μ₁/|μ₀| ≤ λ₁(Q^ev)/|μ₀| (interlacing). Barta itself needs 0 ≤ inf V ≤ λ₀, an accuracy λ₀/λ₁ = 1.3e−4, 3.8e−5, 1.3e−5 (c = 2.5, 3, 3.5) below the Poincaré margin, so the NOGO line's margins understate the difficulty. Exact class: certificates of Q^ev ≥ 0 through a pointwise potential V = Qψ/ψ ≥ 0 and a Dirichlet-form (or any comparison-form) inequality against the pole, for any positive ψ; not covered: exact high-precision eigenvalue enclosures at fixed c (explore39 cites ε_c certified at c = 47).

## Skipped (lower chance; reasons)
- **explore45 (N8, N7), star function vs Levinson's concavity step.** The three freedoms fold back into the barrier's own optimization (argued, no float test):
  - twisted means ∫|VΨ|²|B|²: Littlewood gives ∫log|VΨ| = ∫log|VΨB| − ∫log|B| with ∫_T^{2T} log|B(a+it)| dt = 2πΣ_{B(ρ)=0, β>a}(β − a) + o(T) ≥ o(T) when B(+∞) = 1, so a twist is the mollifier ΨB, capped by the admissible mollifier length and the Levinson–Conrey–Feng optimization;
  - several pieces: per-piece Jensen with local mean squares C_i gains only where the C_i differ, and on every scale where a mean value is proven C_i = C(1 + o(1)); several mollifiers Ψ_j give log|V| = log|VΨ_j| − log|Ψ_j| exactly, so the best one wins;
  - atoms: explore45 measured them (≤ 17 %).
- **explore42 (N5, N8).** Its remaining claims are an identity and fitted laws, not a method limit with unused freedoms. Under RH, for z_k = ρ_k + ε, ∏_{j≠k}|b_j(z_k)| = ∏_{j≠k}|γ_k − γ_j|/|2ε + i(γ_k − γ_j)| = 2ε|ξ′(ρ_k)|/|ξ(ρ_k + 2ε)| (Hadamard product; checked by hand here). λ_min ≈ 1.14σ^{1.84} is a fit shared by ζ, CUE and DH.

## Suggested NOGO lines
- explore49 amendment (ceilaudit8): even q covered: Gauss sums vanish on one parity of k, so the arcs split into odd / q ≡ 2 / q ≡ 0 (mod 4) with coordinates in h′ = h/2 and conics B_δ(4a, q) = 4(P′ − P), B_δ(2a, q/2) = B_δ(a, q) = P′ − P; in every class generic conics hold ≤ 2 counted pairs, the class-normalized Huxley–Watt ratio is 0.036, true/random 1.20–1.24, all4/x2x3·H² 0.056–0.070 (t = 10¹⁶–10²², α = 0.42–0.5, Q = R, 4R, V-route); mixed-parity pairs are never needed (the class split costs ≤ 3); class: fibrewise determinant counts of Bourgain's second-spacing pairs, every q [CeilAudit8](Riemann/Deep/CeilAudit8/REPORT.md)
- explore29 amendment (ceilaudit8): quantizer conditions add no RH information beyond N6's inequalities: Fleischer's hypotheses are L₁(Ξ^{(m)}) ≥ 0 for every m; in the moved-ζ-pair model L₁(Ξ) is blind beyond δ ≈ 1.1 ḡ while L₁(Ξ^{(m)}) fails at m ≈ 2.4δ/ḡ (to 3 ḡ at m = 7) and L_n(Ξ) at n ≈ 2.5δ/ḡ, yet 2-point Lloyd–Max uniqueness of |Ξ^{(m)}|^α (α = 1, 2; 66 sub-intervals) held in every case with δ ≥ ḡ; class: quantizer conclusions (stationarity, uniqueness) for Laguerre–Pólya functions built from Ξ [CeilAudit8](Riemann/Deep/CeilAudit8/REPORT.md)
- explore39 amendment (ceilaudit8): the Poincaré margin is intrinsic: for every positive ψ with Barta's V = Qψ/ψ ≥ 0, P − R ≤ min(λ₀/W_ψ(ψ₀), λ₁/W_ψ(φ₁)) (W_ψ = Var_ν(·/ψ); a comparison form B ≥ κB₀ needs B₀(φ₁) ≤ λ₁/κ); the Markov part's Perron vector fails Barta (inf V = −0.11 to −0.15), and Barta needs inf V ∈ [0, λ₀] with λ₀/λ₁ = 3.8e−5 at c = 3, finer than the margin; class: potential-plus-Dirichlet-form (or comparison-form) certificates of Q^ev ≥ 0 for any ψ; not covered: fixed-c interval eigenvalue enclosures [CeilAudit8](Riemann/Deep/CeilAudit8/REPORT.md)

RH relevance: N9 (B1; DH obeys the same counts, so no RH information), N6 and N8 (B2: the RH content sits in N6's Laguerre inequalities, not in quantizers), N1 and N13 (B3: the binding quantities are λ₀, λ₁ of Weil's even form on [−a, a]).

## Files
- spacing3.c: spacing2.c plus a parity class selector (argument 10: 0 odd, 1 q ≡ 2 mod 4, 2 q ≡ 0 mod 4, 3 mixed). Compile with `cc -O2 -o sp3 spacing3.c -lm`, run as `./sp3 log10t nu Qmult nshuffle seed Vroute verbose alpha - cls`.
- gauss_check.py: the parity and phase formulas of the Gauss sums (87,804 cases).
- b2_laguerre.py: Laguerre orders via the y-expansion.
- b2c_turan_stable.py: L₁(Ξ^{(m)}) in a stable Taylor form. A first Bell-polynomial version lost precision near zeros (its control failed at m = 4) and was deleted.
- b2d_quant.py: two-point Lloyd–Max uniqueness on |Ξ^{(m)}|^α.
- b3_transform.py: the margin for eight choices of ψ, using explore39's xform.py.
- logs/: every run's output and its /usr/bin/time footprint. The peak was 180 MB. Two drafts were killed by memcap at 2.0 GB and 1.1 GB before chunking.
- Run the Python files with `uv run --no-project --with numpy --with scipy python3 <file>` from the repo root, single-threaded (each takes ≤ 5 min).
- No OPEN item, so no novelty search. The B3 lemma is a three-line routine bound (COROLLARIES level).
