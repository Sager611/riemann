# Smale's mean value constant, arXiv 2610.21472 (Jatar–Ng): audit and hunt for an explicit counterexample

Date 2026-10-22. Brief: `.scratch/orca/smale-brief.md`. Verdict: **NO DISCOVERY** (see REPORT.md).
Notation as in the paper: P(0) = 0, P'(0) = 1, S(P) = min_{P'(c)=0} |P(c)/c|, 𝓜 = sup S, 𝓑 = sup 1/B(0) over admissible
univalent Blaschke pairs (W, B), 𝓐 = sup inf_{∂Ω} |f(w)/w| over admissible (f, Ω).
Scripts (this folder): `smale_exact.py` (sympy, 4 s), `smale_numerics.py` (numpy + scipy; modes dieudonne, cshape, nmin,
slp, sym; each < 4 s and < 120 MB).

## 0. The earlier disproof ("a new disproof")

The paper's reference [Lean] is github.com/tadamcz/mean-value-problem (created 2026-09-03): a Lean 4 / Mathlib proof of
`MeanValueProblem.mean_value_problem.disproof` (negation of the Formal Conjectures statement, constant 1), written by a
pre-release GPT-6 Astra in Epoch AI's LeanOpenProblems run `wikipedia-vega-1000usd`, checked by Comparator (only standard
axioms). With z = 0 and p'(0) = 1 the negation says exactly: some p has |p(c)/c| > 1 at every critical point, i.e. 𝓜 > 1.

Construction (file `MeanValueProblem/Resolutions/MeanValueProblem.lean`): q = (z⁴ + 64i)/√4097, K = {|q| ≤ 1} (star-convex
about 0), explicit conformal map W: K → D̄, base point z₀ = 33/32 + 2i, real Blaschke factor B_r, model
F = (z − z₀)·B_r(W z)/B_r(W z₀) (`modelF`), then P = (X − z₀)g + ∫q^{2n}, translated and normalised.
- r comes from `disk_exclusion` (an ∃ statement), g from `exists_base_polynomial` (compactness + C¹ approximation, ∃),
  n from `∀ᶠ n in atTop`. **No explicit polynomial**; the README says "very large, unspecified degree" and
  "violates the bound by a small margin". The paper says this route gives at best about 1 + 2⁻³⁵.
- The approximation is also very slow in practice. W has branch points where 1 + 64i z⁴ = 0, i.e. z⁴ = i/64, just
  outside K. There the Green function of K is (1/4)·log(|z⁴ + 64i|/√4097) = (1/8)·log(4097/4096) = 3.05·10⁻⁵, so the
  Bernstein–Walsh rate is e^{−3.05·10⁻⁵ N}: about 10⁵ degrees for three digits in g alone.

Other sources searched: the arXiv API ("mean value conjecture", 2025–26) finds only Jatar–Ng 2608.27047 (a survey) and
2510.16875 (dual conjecture, odd polynomials). github.com/openai/math: 749 preprint folders and CONTENTS.md, none on
Smale's mean value problem. Zenodo: a non-Archimedean variant (10.5281/zenodo.22226153) and a degree-6 computational
record (10.5281/zenodo.23068020). None of these is a disproof, and none gives an explicit polynomial.

## 1. Audit of 𝓜 ≥ 𝓑 and of the pair (9/16, 3/8)

Chain: 𝓜 ≤ 𝓑 (Prop 3.1), 𝓑 ≤ 𝓐 ≤ 𝓜 (Props 4.2/4.3, Thm 4.4). Each step was checked line by line.
- Lemma 2.1 (regularisation): J_t = {|B| < t} is an analytic Jordan domain once t > max |B(crit B)|. The Riemann map
  ψ_t gives (W∘ψ_t/ψ_t'(0), B∘ψ_t/t), which is regular with B_t(0) = B(0)/t. Implicitly t > B(0) is also needed (0 ∈ J_t);
  the lemma's "t close to 1" covers this.
- 𝓑 ≤ 𝓐: P̃(w) = w·B(W⁻¹w)/B(0) is holomorphic near W(D̄) with P̃'(W z)·W'(z) = (WB)'(z)/B(0) ≠ 0, and
  |P̃/w| = 1/B(0) on W(T). Correct.
- 𝓐 ≤ 𝓜 (Prop 4.3): Runge gives p ≈ f with p' ≠ 0 inside and |p/w| > A + c on a thin annulus N. Hilbert's lemniscate
  {|q| = 1} is a single Jordan curve in N; this is standard via Faber polynomials of the exterior map (the paper cites
  Nagy–Totik). For H_m = p + w q^m:
  - inside {|q| ≤ 1−r}: H_m' → p' ≠ 0;
  - outside {|q| ≥ 1+r}: |q^{m−1}(q + m w q')| ≥ (1+r)^{m−1}(m|wq'| − |q|), which dominates p' (and the case b_m → ∞
    is handled);
  - so all critical points lie in the collar, where q^m(b) = −p'(b)q(b)/(q(b) + m b q'(b)) = O(1/m) uniformly (|bq'| is
    bounded below there);
  - finally divide by H_m'(0) = 1 + q(0)^m → 1.
  **No false or unjustified step.**
- Cor 4.5 (𝓜_d > 1.1018 for all large d) also uses 𝓜_{d+1} ≥ 𝓜_d (Crane). That holds: perturbing by εw^{d+1} adds one
  critical point near ∞ with a huge |P/w|.

The pair x = 9/16, y = 3/8 (Thm 6.1). The exact recomputation is `smale_exact.py`:
- u = −x − iy, Q = (z−1)⁴ − u⁴, a = 1+u, b = 1−iu, H = (1 − āz)Q/(z − a). H matches the paper's a₀..a₄, and
  1/B(0)² = 1 + (1−|a|²)(1−|b|²)/|a−b|² = **8081/6656** exactly.
- Dieudonné polynomial P_c(w) = Σ_{k=1}^4 a_k U_{k−1}(c) w^{k−1}. Schur–Cohn transforms with c symbolic:
  - δ₁, δ₂, δ₃ have degrees 6, 12, 24 (the paper's D₁, D₂, D₃), and δ₁·2³² = D₁ exactly;
  - |lead T₁|² has degree 8 (L₁) and |lead T₂|² has degree 22 (c²M);
  - Sturm counts on [−1,1]: 0, 0, 0, 0 roots, and only c = 0 for |lead T₂|²; all values are positive at c = 1.
- Degenerate values of c:
  - c = 0: |a₁|² − |a₃|² = 79390996785/2³² > 0;
  - c = ±1/√2: δ₁ is the same, the lead is ≠ 0, and δ₂ = 11066090771698524513/2⁶⁴.
  All three equal the paper's numbers.
- Hence H is univalent on D and the certificate is correct. Numerically, min_c min |root P_c| = 1.00221 (at c = 0.7327),
  and H' has zeros at |z| = 1.00387, 1.38843, 1.60486, so H is univalent on D with very little room.
- (WB)' = const·Q'(φ)φ' with Q' = 4(z−1)³ ≠ 0 in D, so the pair is admissible. It is not regular: the model
  f = const·Q∘H⁻¹(N·) has a critical point at w* = H(1)/N ∈ ∂Ω, where |f(w*)/w*| = 1/B(0) exactly (numerically
  1.1018586788). Lemma 2.1 is what makes this harmless.
- The auxiliary quartic normalised at b has S = 1/4: its only critical point is z = 1, and |Q(1)|/(|1−b|·|Q'(b)|) =
  |u|⁴/(4|u|⁴).

**Independent confirmation of 𝓜 > 1:** the Lean proof of §0.

## 2. Searches for an explicit P with S(P) > 1 (all negative)

(a) **Critical-point parametrisation.** P'(w) = ∏_{j=1}^{d−1}(1 − w/c_j) and P(c_j)/c_j = ∫₀¹ ∏_i(1 − τc_j/c_i)dτ.
Gauss–Legendre with ⌈d/2⌉ + 4 nodes is exact. The objective is max min_j log|P(c_j)/c_j|, solved by an LP trust-region
method (SLP), with analytic gradients checked against finite differences to 2·10⁻¹⁰. Results from random starts:
- d = 4, 5, 6, 8, 12 (6 starts each): best S = (d−1)/d, which is Smale's w − w^d/d;
- d = 16, 24, 32, 48 (4 starts each): 0.9375, 0.9500, 0.9526, 0.9655.

So nothing beats (d−1)/d. The mode is `smale_numerics.py slp D TRIALS`.

(b) **Exact-certificate family** P(w) = w·R(w^k).
- Here P'(w) = T(w^k) with T = R + k u R'. So the critical points are exactly the k-th roots of the zeros t_l of T, and
  S(P) = min_l |R(t_l)|.
- Writing T = ∏(1 − u/t_l) gives R_j = T_j/(1 + jk) and R(t) = (1/k)∫₀¹ τ^{1/k−1} T(τt) dτ. A hit could therefore be
  certified in exact Q(i) arithmetic with no root finding.
- As k → ∞, R(t) = 1 + L(t)/k + O(k⁻²) with L(u) = ∫₀^u (T(s)−1)/s ds. The limit problem is Smale's problem for the
  entire function G = u·e^{L(u)}: G' = e^L·T, and |G(c)/c| = e^{Re L(c)}.
- deg L = 2 is impossible. With x_i = λ₁c_i one gets (1+x₁)(1+x₂) = 1, which is incompatible with Re x_i > 1.
- SLP over N ≤ 20 orbits, k ∈ {1, 3, 4, ∞}: the optima are Smale's extremal (−1/N in the limit) or worse.
  The mode is `sym N K TRIALS`.

(c) **Relaxed pairs on a disk.**
- Setup: polynomial f of degree k ≤ 16, a disk D(w₀, ρ) ∋ 0 with |w₀|/ρ ≤ 0.85, |f/w| > A on the circle, and critical
  points of f inside allowed if |f(c)/c| > A. By Prop 4.3 with q = (w−w₀)/ρ, A > 1 would give an explicit
  P = (f + w q^m)/(1 + q(0)^m).
- Best value 1.000 (f = w), checked on 4·10⁵ circle points, from 10 starts each for k = 4, 6, 8, 10, 12, 16.
- An earlier apparent 1.0035 was undersampling near w = 0 and is discarded.

(d) **Astra-type lemniscates** K = {|zⁿ − 1| < ρ}, using the exact conformal map V(s) = s((ρ²−1)/(ρ − sⁿ))^{1/n}.
For n = 2, ..., 6 and ρ ∈ [1.02, 1.3], 6000 random (ζ₀, α) plus hill climbing found no admissible pair with A ≥ 1.003.
Zeros were counted by an adaptive argument principle. Earlier "A = 1.047" hits came from aliasing in a fixed-grid
winding count.

(e) **Polynomial 𝓐-pairs.** R*(p) is the level at which the component of {|p/w| < R} containing 0 first meets a
critical point of p; it is computed with a merge tree of saddles via Newton flow. For degrees 4–6 (300 random starts
plus hill climbing), R* = 1 always.

(f) **Side result (numerical):** a new degree-one admissible pair.
- W = e^{3.1z} is univalent because |Im 3.1z| < π. Take base b = −0.1814 + 0.6221i and zero a = −0.0555 − 0.8659i.
- ((W − W(b))B_a)' has 0 zeros in D: robust winding, with min |·| on T = 8.4·10⁻⁴. So 1/B(0) = 1/|B_a(b)| = **1.03164**.
- This beats the paper's discarded examples (odd Koebe 1.02257, bent slit 1.02553) but not 1.1019.
- Mode: `cshape`.

## 3. Why an explicit P is expensive (degree estimates for the indirect route)

Suppose all critical points lie on one curve Γ around 0 and P' ∝ f' inside Γ, where f is the model, locally univalent
inside. Then the counting measure ν (mass n) satisfies U^ν = log|f'| + C on Γ. Write ν = σ + nω, where σ is the
zero-mass single layer reproducing log|f'| and ω is the equilibrium density. Positivity gives n ≥ n_min = max_Γ(−σ/ω).
- For a circle |w − w₀| = ρ, n_min = 2·max Re((w−w₀) f''/f'), so disks cost nothing here. They fail in §2(c) for a
  different reason.
- Fjords are what cost: ω decays exponentially in them. Computed with Kress log-quadrature (validated: ellipse
  equilibrium density exact to 10⁻¹³, circle n_min = 2|λ| exact), on Γ_t = W({|B_a| = t}):

| pair | min_Γ \|f/w\| | n_min | min ω (relative) |
|---|---|---|---|
| paper (9/16, 3/8) | 1.01 | 9.3·10⁶ | 1.9·10⁻⁶ |
| paper (9/16, 3/8) | 1.02 | 4.9·10⁷ | 4·10⁻⁷ |
| paper (9/16, 3/8) | 1.05 | 2.5·10¹¹ | 1.2·10⁻¹⁰ |
| e^{3.1z} pair of §2(f) | 1.003 | 4.7·10⁴ | 4.7·10⁻⁵ |
| e^{3.1z} pair of §2(f) | 1.01 | 3.9·10⁵ | 6.3·10⁻⁶ |
| second e^{3.1z} pair (A = 1.0197) | 1.003 | 1.4·10⁴ | |

The paper's base point sits at the centre of a large spiral (Ω = H(D)/N has diameter about 60, while dist(0, ∂Ω) is
about 0.5). The construction additionally needs the collar error |p'(b)||q(b)|/|q(b) + m b q'(b)| ≈ |f'|/(2πnω|w|) to
be below the margin min_Γ|f/w| − A. That pushes n higher still, so every explicit instance from the known pairs has
degree ≳ 10⁴. It also needs polynomial approximants p of models with singularities near the fjord, as with Astra's
3·10⁻⁵ in §0.

**Answer to (3):** the indirect argument does yield polynomials, for all large degrees, and the Lean proof confirms
𝓜 > 1. Neither yields an explicit one.

## 4. A degree-independent certificate (for future use)

For P = (p + w q^m)/(1 + q(0)^m) with p, q in Q(i)[w], p(0) = 0 and p'(0) = 1:
1. **Inner curve** C_in ⊂ {|q| < 1}: if p' has no zero inside C_in and |p'| > |q|^{m−1}|q + m w q'| on C_in, then
   by Rouché P' has no zero inside C_in.
2. **Outer curve** C_out ⊂ {|q| > 1}: if |q|^{m−1}|q + m w q'| > |p'| on C_out, and all zeros of q and of q + m w q'
   lie inside it, then all m·deg q critical points lie inside C_out.
3. **Collar bound:** at any critical point b between the two curves,
   |P(b)/b|·|1 + q(0)^m| ≥ |p(b)/b| − |p'(b)|·|q(b)|/|q(b) + m b q'(b)|.
   So an interval-arithmetic lower bound of the right side on boxes covering the collar proves S(P) > A.

The cost of this certificate depends only on deg p and deg q, not on m. What is missing is a pair (p, q) of moderate
degree with positive margins. A k-symmetric variant (§2b) would make even the last step exact.
