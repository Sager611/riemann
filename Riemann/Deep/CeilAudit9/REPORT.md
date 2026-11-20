ROBUST: explore62 (N8) holds in its exact class (𝒮 + normalized microscopic four-point means, uniform in the shifts + kurtosis + one-point law of Z), but its stated witness is separated by two untwisted inputs, macroscopic K(α) (min K_ζ = 0.29–0.34 < 1, repaired unconditionally by Cox-placed windows, argued) and the one-point law of Z′ (the bulk carries 42–54% of ⟨Z′²⟩, repaired only under a hump-level CLT); explore30 (N1) holds and is stated too weakly: two bumps are optimal for periodic blocks (|g(y)| ≤ g(0)/2 for |y| ≥ a), and LP placement under its own constraints gives admissible integer blocks with W < 0 on c = 13.46 at 3·10¹² in its model (LP thresholds there c ≈ 5.2–11.0 by end and rounding model; anchored LP at t → ∞: W < 0 at c = 13.46, not at 11.02) instead of c ≥ 30 (47 at every height); explore45 (N8, N7) holds: twists remove ≤ 4.1% of the concavity loss (a plain mollifier of the same length 8.1%), splits ≤ 1.6%, atoms add nothing beyond its T3.

# ceilaudit9 — falsification of explore62 (N8), explore30 (N1), explore45 (N8, N7) (Orca, 20 Nov 2026, 13:59–15:20 UTC)
Order explore62, explore30, explore45; none skipped. Compute ≈ 41 min (e30c's 20 min included), one process at a time, memcap 290 MB (peaks ≤ 118 MB, RSS ≤ 149 MB).
Lessons applied: each barrier gets its step, its exact input class, the freedoms its own lemmas admit, one float test.

## A. explore62 (N8): shift-dependent fourth-order means cannot lift the 1/√3 ceiling
Step fixed: an argument whose only inputs on Z are limits of untwisted means forces no positive proportion of line zeros (floor 0, dilution).
Class: 𝒮 (second-order unfolded means, derivatives ≤ J) + normalized microscopic four-point means (fixed shifts, derivatives) + kurtosis
M4 ~ L²/(2π²) + (Thm 3′) one-point law of F below b = L^{1/4}/log L. Witness: independent CUE windows of vanishing length + Gaussian or Φ(Y) bulk.
Freedoms:
1. Uniformity in the shifts: normalized four-point functions agree uniformly: fixed x (Hughes–Young leading term = CUE), 1 ≪ x ≪ L (both
   ≈ 3/(π²x²) decay), x ≍ L (both O(1/L²) normalized, difference O(1/M4)). Inside the class.
2. How fourth-order means enter zero counting: unmollified arguments use them through normalized ratios (Hölder ∫|Z| ≥ (∫Z²)^{3/2}(∫Z⁴)^{−1/2},
   Paley–Zygmund): inside the class. Levinson/Conrey/Selberg use mollified (twisted) second moments: outside (explore62's escape (1)). Compatible.
3. Macroscopic shifts: shifts α ≍ 1/log T are x ≍ 1 (inside); α fixed in t gives the unnormalized K(α) = ⟨Z²(t)Z²(t+α)⟩/⟨Z²⟩²: outside.
4. Fourth moment with its secondary terms: P₄(L)/P₁(L)² is a scalar sequence, matched exactly by tuning p (or λ) at every V; the 1/L shape
   corrections of the four-point function are not limits (outside the class); ζ-windows (explore62's variant, length ≫ T^{1−δ₀}) carry them.
5. Joint one-point law of (Z, Z′): outside (Thm 3′ copies the law of Z only).
Tests (e62.py, e62b.py; 3000 units, 32 points per unit, Riemann–Siegel from ExploreMalliavin/zdil.py):
- Freedom 3: K_ζ(α) min 0.343/0.310/0.294 at α ≈ 3 (t = 10⁴/10⁶/10⁸), < 1 on 57–67% of α ∈ [0.2, 12]; the stated witness has K ≥ ≈ 1 there
  (independent placement; Gaussian background 1 + 2sinc²), so K(α) separates it. Repair (argued, unconditional): place Thm 3′'s energy windows
  (t-length 2√(6p) → 0) as a Cox process with intensity ∝ |Z|² smoothed at t-scale σ, √p ≪ σ → 0: K_w(α) = (1−ε)²⟨Λ_σ(t)Λ_σ(t+α)⟩ → K_ζ(α)
  at fixed α ≠ 0; cross-window microscopic terms ≍ √p/σ → 0; inside one window nothing changes.
- Freedom 5: the bulk {|F| < b} carries 54/46/42% of ⟨F′²⟩, 65/62/60% of ⟨|F′|⟩, 73–74% of ⟨|F′|^{1/4}⟩ (vs 11/8/7% of ⟨F²⟩) at 10⁴/10⁶/10⁸;
  Thm 3′'s bulk Φ(Y) has F′ → 0, so its ⟨|F′|^m⟩ is the windows' share only: the stated witness misses the one-point law of Z′ (separation is
  unconditional only if a lower bound for fractional moments of ζ′ is proven; not checked).
  The joint law alone forces nothing: F̂ = σ_block·|F| (sign constant on blocks of 50 humps) keeps 2.0% of the zeros with TV 0.050/0.048/0.035
  between the 16×16 histograms of (F̂, F̂′) and (F, F′), at the level of ζ's own reflection asymmetry TV((F,F′),(−F,−F′)) = 0.048/0.032/0.039.
  With 𝒮's F″ mean the kinks cost: rearranging only small humps (max|F| < b; large humps keep their zeros) keeps 20.4/19.6/18.3/16.3% of the
  zeros (10⁴/10⁶/10⁸/10¹⁰) with kink budget Ω/⟨F″²⟩ = 0.074/0.059/0.052/0.050, Ω = (1/W)Σ_kinks F′(γ)²; a bounce of width δ costs (3π/8)F′(γ)²/δ
  of F″ energy and moves the joint law on ≈ 2δ per kink, so the repair is admissible in the limit iff Ω = o(⟨F″²⟩) and #large humps = o(N):
  explore62 Prop 4's unproven hump-level input (Selberg CLT at relative extrema). Finite-T trend: slow decay.
Verdict: ROBUST in the exact class. Exact class after this audit: 𝒮 + normalized microscopic four-point means (uniform in the shifts) +
kurtosis + one-point law of Z + macroscopic K(α) (Cox witness); conditionally on the hump CLT also the one-point law of (Z, Z′).

## B. explore30 (N1): window records are not zero verification
Step fixed: a zero-side certificate (zeros verified below 3·10¹², B_2, Trudgian's S and S₁ bounds, the count of 2606.31441) for Weil positivity
on the window [−log c, log c]. Class (adversary): off-line pairs δ = 0.37 above 3·10¹², two-bump f, zeros spread over the negative phases with
full gaps elsewhere; stated thresholds c ≥ 30 (3·10¹²), c ≥ 47 (every height).
Freedoms: (i) other f (Chekmachev's extremals, prolate windows, the window ground state); (ii) placement inside a period; (iii) longer
(multi-harmonic) periods; (iv) δ per zero (δ = δ_max is optimal wherever the kernel is negative: adds nothing).
(i) Lemma A (Cauchy–Schwarz): f ∈ L²(−a, a), g = f * f^* ⇒ |g(y)| ≤ g(0)/2 for |y| ≥ a (the pieces [y−a, a] and [−a, a−y] are disjoint).
   For blocks periodic in γ with period P ∈ (π/a, 2π/a] only the Bragg peaks 0, ±y₁ (y₁ = 2π/P) meet supp g, so
   W_norm ≥ 1 − |Σ_ρ cosh(δ_ρy₁)e^{iy₁γ_ρ}|/n (n = DP), with equality for two narrow bumps at ±y₁/2: no f beats the two-bump f on such blocks,
   and the threshold is a property of the placement. (Reviewer: correct for infinite periodic blocks; finite blocks need boundary control.)
(ii) Test (e30.py, e30b.py): LP over the phase distribution (240 bins; |S| ≤ 5.661 on both sides of every atom; ∫S = 0 per period; |S₁| ≤ 2.762;
   D = 4.280: explore30's values and 1-unit margins). W_norm < 0 for c > 5.15 (free S, S₁ offsets), c > 7.09 (+1 unit for integer rounding),
   c > 7.58 (anchored: block ends at S = S₁ = 0) at 3·10¹²; at t → ∞ (S/D ≤ 0.594, S₁/D ≤ 0.371) c > 10.1 (w → 0: 9.46); with margins
   7.12 (10¹³), 7.38 (10²⁰), 8.41 (10¹⁰⁰), 9.42 (10³⁰⁰); anchored at t → ∞ (e30d.py, 16 phases): W_norm +0.098/−0.035/−0.170/−0.340 at
   c = 11.02/13.46/16.44/20.09 (explore30 at t → ∞: c ≥ 47; narrow-bump limit 35.2).
   Integer blocks (e30v.py: 824 periods, second-order sigma-delta rounding, clusters split symmetrically, exact two-bump kernel K = Re[f̂(γ−iδ)
   conj f̂(γ+iδ)], w = 0.03a, on-line lattice outside, S and S₁ exact including the extrema between zeros), explore30's constraints:
   c = 13.46 (a = 1.3): W_norm −0.101, max|S| 5.38 ≤ 5.66, S₁ ∈ [−1.56, 2.58] ⊂ ±2.76: admissible;
   c = 20.09 (a = 1.5): W_norm −0.394, max|S| 5.00, S₁ ∈ [−1.11, 2.63]: admissible (explore30's family: not admissible below c = 30).
(iii) Not tested (P > 2π/a lets g(y₁) approach g(0) but doubles the zeros per period, which the S₁ bound penalizes).
Prime side: a block of period P is invisible to every window with log c′ < 2π/P (its Bragg peak lies outside supp g) up to its ends, so it agrees
with the explicit formula on those windows as well as the uniform zeros do; on its own window the prime side (W_true ≈ 1 + O(Σ_{n≤c}Λ(n)n^{−½}/log t₀))
excludes it (for the c = 13.46 block's f only n = 13 meets supp g, |prime term| ≤ Λ(13)13^{−½}/(2πD) = 0.026, so W_true ≥ 0.97 > −0.10).
So explore30's converse extends: certified windows c ≥ 13.46 already carry Euler-product information beyond zero verification.
Reviewer (verdict GAP, logs in this section): Lemma A and the kernel identity correct; gaps: finite-block optimality (only infinite periodic
blocks are covered by Lemma A), interior S₁ extrema (covered by the exact integer check), and the splice to ζ's true zeros (D·L ∉ ℤ, actual S,
S₁ at the block ends, the exterior count): explore30's own 1-unit-margin model, so c = 13.46 is a result in its model, not an unconditional block.
Verdict: ROBUST and stated too weakly: the barrier holds on c ≥ 13.46 at 3·10¹² in explore30's model (integer floats; LP c ≈ 5.2–11.0) and on
c ≥ 13.46 in the t → ∞ limit (anchored LP, floats; intermediate heights not rechecked);
30 and 47 are thresholds of the full-gap placement, not of the test-function class. Certificate steps for the stronger line: splice to ζ's
computed zeros at a block start (S_true, S₁_true), integer closure of S (one adjusted period), Arb evaluation of W and of S, S₁.

## C. explore45 (N8, N7): symmetrization and star-function inequalities cannot beat Levinson's concavity step
Step fixed: Jensen ∫log|VΨ| ≤ (T/2) log mean|VΨ|² on the Levinson line. Class: inequalities valid for every subharmonic function + the line mean
squares C(x) of one fixed VΨ; T3: analytic competitors with the same 42 line mean squares (zero fraction 0.200 vs smeared 0.2407).
Freedoms (ceilaudit8 argued them, untested): twisted means ∫|VΨB|², several pieces, the zeros' atoms with the star function's freedom.
Tests (e45.py; T₀ = 10⁴, H = 400, dt = 0.02, explore45's θ = ½, P(x) = x, b = 1.03217, R = 1.24618): concavity loss J = 0.1727 (explore45: 0.179).
- Twist B = 1 + Σ_{2≤k≤K} β_k k^{−s}, β complex, optimized in sample: U_B = ½ log mean|VΨB|² − mean log|B| removes 0.76/1.69/2.37/4.08% of J
  (K = 2/3/4/6); the plain Levinson mollifier of length Ky removes 5.60/7.76/8.06% (K = 2/3/4): by ∫log|VΨ| = ∫log|VΨB| − ∫log|B| the twist is
  the mollifier ΨB, and a worse one than the plain mollifier of the same length. First-order gain ∝ c₂ = 0.147 (coefficient of 2^{−s} in VΨ),
  → 0 as T → ∞ at fixed θ.
- Two pieces, per-piece Jensen: E = {cos(t log p + φ₀) > 0} and complement (best of 16 φ₀): 0.59/0.74/1.57% of J (p = 2/3/5); these pieces are
  twisted weights (Fourier series in p^{it}), same order as the twists; two halves in t: 0.032%.
- Atoms: T3's competitor is log|f| for an analytic f, so every inequality valid for subharmonic functions (star function, symmetrization) holds
  for it with the same data: the star function's freedom plus atomicity proves nothing below 0.200, i.e. adds nothing beyond T3 (≤ 17%).
Verdict: ROBUST in the class; ceilaudit8's two arguments hold in floats; the remaining gains are mollifier changes, capped by the admissible length.

## Suggested NOGO lines
- explore62 → ceilaudit9: ROBUST in 𝒮 + normalized microscopic four-point means (uniform in the shifts) + kurtosis + one-point law of Z; its
  witness fails macroscopic K(α) (K_ζ min 0.29–0.34 at α ≈ 3, t = 10⁴–10⁸; Cox-placed energy windows repair it, argued) and the one-point law of Z′
  (the bulk carries 42–54% of ⟨Z′²⟩; sign-rearranged |Z| matches the (Z, Z′) law at 2% of the zeros but needs a hump-level CLT to keep 𝒮's
  F″ mean: kink budget 0.074 → 0.050 of ⟨F″²⟩ from 10⁴ to 10¹⁰) [CeilAudit9](Riemann/Deep/CeilAudit9/REPORT.md)
- explore30 → ceilaudit9: stated too weakly: |g(y)| ≤ g(0)/2 for |y| ≥ a makes two bumps optimal for periodic blocks; LP placement under
  explore30's own constraints and margins gives admissible integer blocks with W/(D·2π‖f‖²) = −0.10 on c = 13.46 and −0.39 on c = 20.1 at
  3·10¹² (LP thresholds c ≈ 5.2–11.0 by end and rounding model; anchored LP at t → ∞ negative from c = 13.46): zero verification certifies no window c ≥ 13.46 in this model
  [CeilAudit9](Riemann/Deep/CeilAudit9/REPORT.md)
- explore45 → ceilaudit9: ROBUST for subharmonic inequalities + line mean squares of a fixed VΨ: at T = 10⁴ twists B = 1 + Σ_{k≤6}β_k k^{−s}
  remove ≤ 4.1% of the concavity loss (plain Levinson mollifier of the same length: 8.1%), two-piece splits ≤ 1.6%, time halves 0.03%; atoms
  add nothing beyond T3 [CeilAudit9](Riemann/Deep/CeilAudit9/REPORT.md)

RH relevance: N8 (critical-line proportion: explore62's class sharpened, two witness gaps named; explore45's concavity step confirmed), N1 (Weil
windows: zero verification reaches no window c ≥ 13.46 in explore30's model; certified windows from c = 13.46 up carry Euler-product
information), N7 (the same Littlewood + Jensen step in zero-density detection). No discovery; no NOVELTY search (nothing OPEN).

## Files
Riemann/Deep/CeilAudit9/: e62.py (macroscopic K, bulk shares, sign-rearranged |Z|), e62b.py (small-hump rearrangement, kink budget),
e30.py (LP adversary, anchored variant), e30b.py (bisection thresholds), e30c.py (anchored bisection; killed at 1200 s after its first line),
e30d.py (anchored LP at t → ∞), e30v.py (integer blocks, exact kernel), e45.py (twists, splits, longer mollifiers);
logs/: e62.log, e62b.log, e30.log, e30b.log, e30c.log, e30d.log, e30v.log (first integer runs, rounding/origin bugs), e30v2.log, e30v3.log
(final integer check), e45_small.log (memcap peaks and /usr/bin/time -l inside each log). .scratch/ceilaudit9/: STATE.md, notes.md, dbg30.py.

## Working notes (chronological, kept for the record)
- 13:59 UTC start. Read NOGO 748 (explore62), 608 (explore30), 650 (explore45), 749 (ceilaudit8), ExploreMalliavin REPORT §1–5.
- explore62 class (Thm 3/3′): inputs = 𝒮 (untwisted second-order unfolded means, derivatives ≤ J) + (I4) normalized microscopic
  four-point means (fixed shifts x, derivatives) + (Ik) kurtosis M4 ~ L²/(2π²) + (I1/Thm 3′) one-point law below |Z| = L^{3/4}/log L.
  Witness: independent CUE(N) windows (length fraction p → 0, energy η → 0 or 1−ε) + Gaussian / dead-block / low-frequency Φ(Y) bulk.
  Step fixed: "a criterion whose inputs on Z are exactly these limits forces ≥ c·N(T) zeros on the line": none can (floor 0).
- 14:50 E62 done (e62.py → logs/e62.log, e62b.py → logs/e62b.log; W = 3000 units, 32 points/unit, Riemann–Siegel; peaks 118/109 MB memcap, RSS 142/131 MB).
  (a) macroscopic K_ζ(α) (α fixed in t): min 0.343/0.310/0.294 at α ≈ 3 (t = 10⁴/10⁶/10⁸), K < 1 on 57–67% of α ∈ [0.2, 12];
      stated witness: independent windows + Gaussian/low-frequency bulk give K ≈ (1−η)²(1 + 2sinc²) + window part ≥ ≈ 1 → separated (explore62 escape (4)).
  (b) bulk {|F| < b}, b = L^{1/4}/log L (0.83/0.75/0.72): carries 54/46/42% of ⟨F′²⟩ and 65/62/60% of ⟨|F′|⟩, 73% of ⟨|F′|^{1/4}⟩ (vs 11/8/7% of ⟨F²⟩);
      Thm 3′'s bulk Φ(Y) has F′ → 0, so its ⟨|F′|^m⟩ (m < 2) is the windows' share only → the stated witness misses the derivative's one-point law.
  (c) repair 1 (all humps): F̂ = σ_block·|F| (σ constant on blocks of 50 humps): zeros/unit 0.020 vs 1.000, TV of joint (F,F′) histogram (16×16)
      0.050/0.048/0.035 vs ζ's own asymmetry floor TV((F,F′),(−F,−F′)) 0.048/0.032/0.039 → the one-point law of (Z, Z′) alone forces nothing;
      cost: kinks at every former zero; ⟨F̂″²⟩/⟨F″²⟩ = 7.3–7.7 (grid kinks), smoothing sqrt(F² + d²F′²): 3.8 at d = 0.05, TV 0.09–0.12.
      repair 2 (small humps only, max|F| < b, take the preceding large hump's sign; large humps keep their zeros): zeros/unit 0.204/0.196/0.183/0.163
      (t = 10⁴/10⁶/10⁸/10¹⁰); kink budget Ω = (1/W)Σ_kinks F′(γ)²: Ω/⟨F″²⟩ = 0.074/0.059/0.052/0.050 (all humps: 0.156–0.163);
      smoothing width δ costs (3π/8)Ω/δ of ⟨F″²⟩ (0.87→0.59 at δ = 0.1) and moves the joint law on ≈ 2δ per kink.
  Reading: (F, F′) one-point law + F″ second-order mean (in 𝒮) are compatible with o(N) zeros iff Ω_small = o(⟨F″²⟩) and #large humps = o(N):
  the same hump-level input as explore62 Prop 4 (unproven: Selberg CLT at relative extrema); finite-T trend is slow decay (0.074 → 0.050 over 6 decades).
  Macroscopic K repair (argument): place the energy windows (Thm 3′, energy 1−ε, t-length 2√(6p) → 0) as a Cox process with intensity ∝ |Z|² smoothed
  at t-scale σ (√p ≪ σ → 0): macroscopic K_w(α) = (1−ε)²⟨Λ_σ(t)Λ_σ(t+α)⟩ → K_ζ(α); cross-window microscopic terms ≍ √p/σ → 0; unconditional.
- 15:08 E30 (e30.py LP; e30v.py integer check; e30b.py bisection; logs/e30*.log; peaks ≤ 117 MB memcap, RSS ≤ 149 MB; e30c.py anchored bisection killed at 1200 s after its first line).
  Lemma A (Cauchy–Schwarz): f ∈ L²(−a, a), g = f * f^* ⇒ |g(y)| ≤ g(0)/2 for |y| ≥ a (the overlap pieces [y−a, a], [−a, a−y] are disjoint).
  ⇒ for a block periodic in γ with period P ∈ (π/a, 2π/a] only the Bragg peaks 0, ±y₁ (y₁ = 2π/P) meet supp g, W_norm ≥ 1 − |Σ_ρ cosh(δ_ρy₁)e^{iy₁γ_ρ}|/n
    (n = DP zeros per period), equality for two narrow bumps at ±y₁/2: explore30's two-bump f is OPTIMAL for single-harmonic blocks; Chekmachev extremals,
    prolate windows and the window ground state can only do worse there. The threshold is a property of the zero placement.
  LP over the phase distribution (240 bins, off-line pairs δ = 0.37 where the kernel is negative, |S| ≤ 5.661, |S₁| ≤ 2.762, D = 4.280 at 3·10¹²: explore30's values):
    unanchored (free S, S₁ offsets): W_norm < 0 for a > 0.820 (c > 5.15); with 1 extra unit for integer rounding a > 0.979 (c > 7.09); t → ∞ (S/D ≤ 0.594,
    S₁/D ≤ 0.371): a > 1.158 (c > 10.1; w → 0: 9.46); a*(t) with margins: 7.12 (10¹³), 7.38 (10²⁰), 8.41 (10¹⁰⁰), 9.42 (10³⁰⁰).
    anchored (block ends where S = S₁ = 0, as explore30's ends), margin 0: a* = 1.013 (c* = 7.58) at 3·10¹².
  Integer blocks (824 periods, sigma-delta rounding, clusters split symmetrically, exact two-bump kernel with w = 0.03a, on-line lattice outside):
    a = 1.3 (c = 13.46): W_norm −0.101, max|S| 5.38 ≤ 5.66, S₁ ∈ [−1.56, 2.58] ⊂ ±2.76: ADMISSIBLE; a = 1.5 (c = 20.09): W_norm −0.394, max|S| 5.00,
    S₁ ∈ [−1.11, 2.63]: ADMISSIBLE (explore30: not admissible below a = 1.7). Earlier integer runs (logs/e30v.log, e30v2.log) had S₁ drift from rounding/origin bugs, fixed.
- 15:08 E45 (e45.py → logs/e45_small.log; T0 = 10⁴, H = 400, dt = 0.02, 2·10⁴ points, explore45's θ = ½, P(x) = x, b, R; peak 68 MB, RSS 94 MB, 113 s):
  U(Ψ_y) = 0.30763, m = 0.13496, concavity loss J = 0.17267 (explore45 on H = 4000: 0.179).
  twist B = 1 + Σ_{k≤K} β_k k^{−s} (β complex, Nelder–Mead, in-sample): gain 0.76/1.69/2.37/4.08% of J for K = 2/3/4/6;
  plain Levinson mollifier of length Ky (P(x) = x, no twist): 5.60/7.76/8.06% of J for K = 2/3/4 → the twist is a worse mollifier of the same length.
  split by sign of cos(t log p + φ₀) (best of 16 φ₀): 0.59/0.74/1.57% of J (p = 2/3/5); two halves in t: 0.032% of J.
  first-order twisted gain ∝ c₂ = coefficient of 2^{−s} in VΨ = 0.147 here, → 0 as T → ∞ at fixed θ (2^{a−½}P(1 − log 2/log y) → 1, b log 2/L → 0).
