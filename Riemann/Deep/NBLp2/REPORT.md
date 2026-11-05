THEOREM: For every Dirichlet character χ mod q and every 1 < p < 2: L(s,χ) ≠ 0 on Re s > 1/p ⟺ 1_(0,1) lies in the L^p(0,1)-closure of span{κ_χ(1/(kx)) − δ_χ k^{−1} κ_χ(1/x) : k ≥ 1}, where κ_χ(y) = δ_χ (φ(q)/q) y − Σ_{n≤y} χ(n) and δ_χ = 1 iff χ is principal (⟺ that span is dense in the step-function space S_p); hence unconditionally for all χ and 1 < p ≤ 8/7 (Lean-checked 7/8), and 1 < p ≤ 1.1449744 with Bebel's half-plane (node N5).

Node N5 (Nyman–Beurling / Báez-Duarte), Dirichlet L-function analogue. Agent `nblp2`, 2026-11-05. Paper proof, no Lean.
Full proof: [PROOF.md](PROOF.md). It adapts [NBLp/PROOF.md](../NBLp/PROOF.md) (the case q = 1) lemma by lemma.

## 1. The brief's premise was wrong: this does not settle Dimitrov–Oliveira's Problem 1

Problem 1 of arXiv 1608.07887 (TeX lines 247–251, quoted exactly):

> For any fixed n ∈ ℕ and p ∈ (1,2], determine the best approximation of 1_(−∞,∞) by products of the form
> L(1/p+it) A_n(1/p+it) in L²(ℝ, ω), where ω(t) = 1/(1/p² + t²) is the weight function.

That is a fixed-n, weighted-L² extremal problem (its answer is an orthogonal projection, already written there as a
quotient of Gram determinants). It is not an L^p natural-number problem, and nothing here solves it. Round R312's
reading ("ask for a true L^p natural-number endpoint theorem") was a misreading.

What the paper does contain:
* Their functions: κ(x) = β x^α − Σ_{k≤x} χ(k) k^{1/2−1/p}, with α = 3/2 − 1/p and β = φ(q)/(αq) for principal χ,
  α = β = 0 otherwise (lines 203–207). The weights k^{1/2−1/p} shift the line Re s = 1/p to Re s = 1/2.
* Thm GBNew/GBLp (p ∈ (1,2]): L ≠ 0 on Re s > 1/p ⟺ d_n(L,p) → 0 in L²(0,∞). In the proof (lines 787–822) this is the
  p = 2 Báez-Duarte criterion for the shifted series η(s) = L(s + 1/p − 1/2, χ): a Hilbert-space statement.
* Thm GNB (line 654): for their class 𝒟 and every p > 1, η ≠ 0 on Re s > 1/p ⟺ 1 lies in the L^p(0,1)-closure C^p_η of
  continuous dilations κ_η(1/(θx)), θ ≥ 1 (side condition β Σ b_k/θ_k^α = 0). The natural-number closure C^p_{η,ℕ} is
  defined (line 330) but used only for p = 2.

The theorem above is the natural-number (Báez-Duarte) strengthening of their Thm GNB in L^p(0,1), 1 < p < 2, for
η = L(·,χ) itself: in their notation, 1 ∈ C^p_{η,ℕ} ⟺ L(s,χ) ≠ 0 on Re s > 1/p. Its family is the unweighted κ_χ
(their p = 2 normalisation, α = δ_χ, β = δ_χ φ(q)/q); their p-weighted κ would detect zeros only on Re s > 2/p − 1/2
(PROOF.md, Remark 9.3).

## 2. Statement (PROOF.md §1)

Notation: S_χ(y) = Σ_{n≤y} χ(n); κ_χ as above (for q = 1, κ(y) = {y}); r_k(x) = κ_χ(1/(kx)); λ_k = r_k − δ_χ k^{−1} r_1
(λ_k vanishes on (1,∞); for q = 1, λ_k = ρ_{1/k} of NBLp). On I_n = (1/(n+1), 1/n]: λ_k = δ_χ S_χ(n)/k − S_χ(⌊n/k⌋).
* **Theorem A_χ.** For 1 < p < 2 and any χ: (i) L(s,χ) ≠ 0 on Re s > 1/p ⟺ (ii) 1_(0,1) ∈ closure of span{λ_k} ⟺
  (iii) the span is dense in S_p. Under (i), for each ε ∈ (0, ε_0] and N ≥ N_0(p,χ,ε):
  ‖1 + Σ_{k≤N} μ(k)χ(k)k^{−ε}(1 − k/N)λ_k‖_p ≤ A_{p,χ} ε^{1−1/p} + (1 + qδ_χ)ε.
* **Theorem B_χ (open range).** If L(s,χ) ≠ 0 on Re s > θ ≥ 1/2 and p < 1/θ: ‖1 + Σ_{k≤N} μ(k)χ(k)(1 − k/N)λ_k‖_p ≪ N^{−η}
  for every η < 1/p − θ.
* **Necessity** holds for all 1 < p < ∞; so for p > 2 the criterion fails for every χ (unconditionally).
* **Corollaries.** Unconditional for every χ and 1 < p ≤ 8/7 (OpenAI's 7/8 half-plane for all Dirichlet L, Lean port
  `OAI.DirichletCharacter.LFunction_ne_zero_of_seven_eighths_lt_re` in Riemann/Deep/QRH, sorry-free) and for 1 < p ≤ 1/0.87338196661320752
  = 1.1449744 (Bebel, QuasiRiemannTracker, 31 Oct 2026: ζ and every Dirichlet L; Lean-kernel-checked per the tracker,
  axiom list not inspected by us). GRH for χ ⟺ the criterion for all p ∈ (1,2). For principal χ mod q the criterion is
  equivalent to ζ ≠ 0 on Re s > 1/p (a family of ζ-criteria, one per modulus), and the side condition can be dropped
  (Remark 9.1: a conjugate-pair functional).

## 3. Changes against NBLp

| NBLp item | Status | Change for χ mod q |
|---|---|---|
| §0 notation | changed | ρ_{1/k} → λ_k = κ_χ(1/(kx)) − δk^{−1}κ_χ(1/x); κ_χ(y) = δ(φ(q)/q)y − S_χ(y), bounded by φ(q) |
| Prop. C (structure) | changed (values) | λ_k on I_n = δS_χ(n)/k − S_χ(⌊n/k⌋); D_mλ_k = λ_{km} − δk^{−1}λ_m; (b), (d) unchanged |
| Prop. N (necessity) | unchanged | ∫_0^1 λ_k x^{w−1} = −L(w,χ)(k^{−w} − δ/k)/w = 0 at a zero w |
| Lemma T (Mellin) | changed | 𝓜r_k = −k^{−s}L(s,χ)/s by Abel summation; for non-principal χ valid for all σ > 0 (r_k = 0 on x > 1/k) |
| Lemma M (Hölder), Lemma PW (line shift) | unchanged | — |
| Lemma Γ | changed for χ* ≠ 1 | γ = (q*/π)^{(s+a)/2}Γ((s+a)/2): bounds gain a factor (q*)^{ε/2}; no exceptional point s = 1; χ* = 1 unchanged |
| Lemma E (Euler factors) | new | imprimitive χ: L(s,χ) = L(s,χ*)E_χ(s); zeros of E_χ on Re s = 0; \|E_χ(s)/E_χ(s+ε)\| ≤ q^{3ε} on Re s ≥ 1/2 |
| Lemma H (Hadamard midline) | changed | product for ξ(s,χ*) (Davenport Ch. 12; only Re B = −Σ Re 1/ρ enters, so complex χ is fine), times Lemma E |
| Lemma A (Lindelöf in zero-free region) | changed | convexity \|L\| ≤ (2 + 5C_κ)(\|t\|+2)^{1/2} via Abel summation; no slit for non-principal χ; t < 0 via χ̄ |
| Lemma P (Riesz–Perron) | changed | μ(k) → μ(k)χ(k) |
| Thm B | changed | the r_1 correction appears only for principal χ |
| A1 (functions) | changed | E_{ε,χ}(y) = Σ_{j≤y}χ(j)φ_ε(j) − δ(φ(q)/q)y/L(1+ε,χ); H_ε = −Σ_{j≤m}χ(j)φ_ε(j) + δS_χ(m)/L(1+ε,χ) on I_m |
| A2, A3 | unchanged argument | bounds through Lemmas H_χ, Γ_χ; q-dependent constants |
| A4 (key estimate B(ε) ≤ Kε) | changed | pole term only for principal χ (residue φ(q)/q, and (φ(q)/q)/L(1+ε,χ_0) ≤ ε); **new step** for non-principal χ: R(2σ') may be complex, \|1 − R(2σ')\| ≤ ε\|ζ'(1/2+1/p)\|ζ(1/2+1/p); C_Γ → C_{Γ,q} |
| A5 | changed | correction δr_1/L(1+ε,χ), at most qε |
| A6, A7 | unchanged | constants depend on χ; N → ∞ at fixed ε, then ε → 0 |
| Corollaries | changed | inputs: 7/8 for all Dirichlet L (Lean port), Bebel's 0.87338… |
| Remark 9.1 | new | principal χ (and q = 1): the side condition Σ b_k/k = 0 is removable in L^p(0,1) |

Verdict on difficulty: routine, as R312 predicted. The only step that is not a change of constants is A4 for
non-principal χ, where R(2σ') ≥ 1 (true for ζ and principal χ) fails and is replaced by |1 − R(2σ')| = O_p(ε).

## 4. Self-review

* **Pass 1** (5 parallel reviewer agents, one per section: setup and Mellin; Γ/E/H; A/P/Thm B; A1–A4; A5–A7 and
  corollaries): all SOUND.
* **Pass 2** (4 adversarial reviewer agents: order of limits and hidden use of the line Re s = 1/p; independent
  re-derivation of the residues; gamma and Euler-factor bounds uniformly on Re s = σ', σ'' and the horizontal sides,
  odd and imprimitive χ, principal exceptional points; setup attacks, Siegel and multiple zeros): all SOUND. One formula
  fixed: Lemma T's bound for ‖r_1‖ is now stated separately for δ = 0 (it was formally 0/0 at σ = 1).
* **Own checks.** Order of limits: N_0(p,χ,ε) is chosen at fixed ε (A7); K_{p,χ} does not depend on ε; d_{N,p} is
  non-increasing in N. The hypothesis is only used on Re s ≥ 1/p + ε/4, never on the line 1/p. Residues: the rectangle is
  clockwise for the right-hand region, J = 2π[R(2σ')/(2σ') − δ(φ(q)/q)/((2σ'−1)L(1+ε,χ))]. Gamma factors: Re ψ(z) ≤ log|z| + 4
  for Re z ≥ 1/4 covers both parities a = 0, 1 on Re s ≥ 1/2.
* **Numeric residue check** (checks2.py part 2; p = 8/7, ε = 0.1, σ' = 0.825; Simpson to |t| = 100 plus tail 2/T):

| χ | J numeric | J formula | difference |
|---|---|---|---|
| χ₋₃ | 3.7269922 | 3.7269979 | 5.7·10⁻⁶ |
| χ₀ mod 3 (with the pole term) | 3.2418085 | 3.2418035 | 5.1·10⁻⁶ |
| complex χ mod 5 (χ(2) = i) | 3.7842012 + 0.0295347i | 3.7841977 + 0.0295297i | 6.2·10⁻⁶ |

With the sign of the s = 2σ' residue flipped the formula would give −3.727, −4.978 and −3.784 − 0.030i. The complex
case confirms that J is genuinely complex and that only Re J enters. Numerically B(ε)/ε = 0.0133, 0.386, 0.0133, far
below the proven K_{p,χ} ≈ 8·10⁴–6·10⁵ (the constants are qualitative, as in NBLp).

## 5. Numerics (nbchi.c, exact integer sieve, M = 10⁹ breakpoints, one process; log run1e9.txt)

‖1 + Σ_{k≤N} μ(k)χ(k)(1 − k/N) λ_k‖_{L^p(0,1)} = (Σ_n |v(n)|^p/(n(n+1)))^{1/p}, v(n) = 1 + δS_χ(n)P_N(1) − Σ_{j≤n} χ(j)A_N(j),
A_N(j) = Σ_{k|j,k≤N} μ(k)(1 − k/N). Intervals: partial sum to M plus the rigorous tail Vmax^p/(M+1).

| χ | N = 10 (p = 1.1) | N = 100 | N = 1000 | p = 2 (N = 10, 100, 1000) |
|---|---|---|---|---|
| χ₋₃ (mod 3, odd) | 0.1975366999 | [0.0501847712, 0.0501848121] | [0.0177404, 0.0177410] | 0.2798, 0.1809, 0.1785 |
| χ₋₄ (mod 4, odd) | 0.1902035455 | [0.0566491282, 0.0566491640] | [0.0206820, 0.0206824] | 0.2702, 0.2174, 0.2216 |
| χ₀ mod 3 | 0.6564771811 | [0.2007612123, 0.2007612834] | [0.0455231, 0.0455242] | 0.6973, 0.2923, 0.1459 |
| χ₀ mod 4 | 0.7521779228 | [0.2338089364, 0.2338089817] | [0.0518532, 0.0518539] | 0.7765, 0.3207, 0.1500 |
| ζ (q = 1, NBLp) | 0.3459480980 | [0.0592261752, 0.0592262243] | [0.0150506, 0.0150513] | 0.3646, 0.1326, 0.1207 |

* N = 10 column: exact infinite sums over the period q·lcm(1..10) (7560 or 10080) with digamma weights (checks2.py part 1);
  each lies inside the M = 10⁹ interval [0.1975366979, 0.1975367021], [0.1902035435, 0.1902035476],
  [0.6564771751, 0.6564771797], [0.7521779169, 0.7521779203]; periodicity was checked over two periods.
* Regression: the code with χ ≡ 1 reproduces NBLp's run digit for digit (M = 10⁶, 10⁷).
* At p = 1.1 the norms fall by 2.7–4.5× per decade (local exponents 0.44–0.65), consistent with the GRH-conditional
  rate N^{−η}, η < 1/p − 1/2 = 0.409 (Theorem B_χ; L(s,χ₋₃), L(s,χ₋₄) satisfy GRH numerically far beyond the heights that
  matter here). Principal characters start higher (their Möbius sums skip k sharing a factor with q) but then fall
  faster. At p = 2 the non-principal norms stall, as for ζ.
* Endpoint family (illustration, M = 2·10⁶, tail not bounded rigorously): for χ₋₃,
  ‖1 + H_ε‖_{1.1} = 0.138, 0.082, 0.046, 0.024 at ε = 0.4, 0.2, 0.1, 0.05, i.e. about (0.35–0.49)ε (ζ in NBLp: 0.360, 0.205,
  0.110, 0.057). The proven rate is ε^{1−1/p}; the true zeros sit on Re s = 1/2 < 1/p, so the decay is faster.

## 6. Novelty pass

Question: has anyone proved the L^p (p ≠ 2) natural-number criterion for Dirichlet L-functions? Not found.
* Citing works of 1608.07887: OpenAlex (W2514882434) lists 2, Semantic Scholar 4: arXiv 1704.01234 (zeros of Dirichlet
  *polynomials* via a density criterion), and three papers on orthogonal Dirichlet polynomials (2019, 2019, 2020). None
  treats L^p closures of natural-number dilations.
* Delaunay–Fricain–Mosaki–Robert (arXiv 1101.1199, Trans. AMS 2013; 1112.0166, Math. Z. 2013): zero-free discs and a
  Beurling–Nyman criterion for L-functions in weighted L²((0,1), dx/x) (Hilbert space, Nikolski-type distances).
* de Roton (C. R. 2005; Bull. SMF 2006; Trans. AMS 2007; JNT 2009): Beurling–Nyman and Báez-Duarte criteria for the
  Selberg class, p = 2 only (as D–O state, line 202). Semantic Scholar: 3 citing works of the 2005 note (1608.07887,
  1704.01234 and one whose title has no L^p, Dirichlet or Nyman–Beurling keyword).
* Báez-Duarte 2003/2005 (ζ, p = 2); NBLp's pass already scanned their 101 citing works for L^p (none).
* arXiv listings math.NT 2026-10 (854 entries) and 2026-11 (140 so far), plus math.CA and math.FA for both months: no
  title on Nyman–Beurling, Báez-Duarte, or L^p closures for L-functions (keyword scan of all titles).
* Zenodo (4 queries): Nyman–Beurling records exist (e.g. 6982337, 17957671, 22647730, 22650745, 23129278, 23370435),
  none mentions L^p or Dirichlet characters.
* Our own files: NBLp (ζ only) and R312 (the idea).
Limits of the search: OpenAlex and the arXiv API throttled (HTTP 429) part-way; the Evidence API was not used.

Assessment: new to our search, but a routine extension of our own NBLp theorem (one genuinely new step, A4 for
non-principal χ). It belongs in COROLLARIES.md (or as a generalisation line on the NBLp discovery page), not as a separate
DISCOVERIES entry; the orchestrator decides.

## 7. Files

* [PROOF.md](PROOF.md): full proof, lemma-by-lemma status against NBLp.
* [nbchi.c](nbchi.c): L^p norms of the explicit approximants (types: ζ, χ₋₃, χ₋₄, χ₀ mod 3, χ₀ mod 4); [run1e9.txt](run1e9.txt) its log.
* [checks2.py](checks2.py): exact N = 10 period sums; residue identity of A4(b) incl. a complex character mod 5; endpoint
  illustration. Run: `uv run --no-project --with mpmath --with numpy python3 checks2.py`; output in [checks2.txt](checks2.txt).

RH relevance: N5
