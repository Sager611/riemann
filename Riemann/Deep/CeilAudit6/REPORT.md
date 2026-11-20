ROBUST: (N5, explore38) exact for zero-kernel test vectors, since P(N) = dist(P_C1, V_N)²; ζ keeps δ = L²(d_N² − B) ≈ .068 flat to N = 4000 and χ₋₄ plateaus at ≈ .085 to 3000, but DH's "first order" is not established: a plateau δ_∞ ≈ .39 fits equally well to N = 2800, and the split between P and the continuous part is unmeasured. (N13, aclgh2) holds against every freedom of Theorem condpos (any block set, any absorption with κ♮ → 0, weights; S is forced), but ACLG's own window model contradicts its semiclassical law: thin pockets give κ ≈ .16 ≫ √(A/2π), and one window with 8 pockets has Σκ² = 1.11 against Σ A/2π = .28. So the stated exponent is not the growth law; a larger global count is a conjecture. (N12, explore48) holds in every rank: rank ≥ 3 lattices with a ζ factor have T_s < γ₁, rank 3's ζ(2s − 1) is only formal, and class-group sums that keep their phases are arithmetic, not geometry. (N1/N13, explore32) step one is RH-blind at the statement level, and the c = 2 barrier is exact for gamma-factor-only bounds and pointwise-sign Perron–Frobenius. (N10, explore33) exact for level-weight and influence statistics whose top-level weight is 2^{o(n)} relative to the total. Float tests were done for N5, N13 and N12; N1 and N10 were argued without one.

# CeilAudit6 — barrier audit of N5, N13, N1, N12, N10 NOGO lines (Orca ceilaudit6, 20 Nov 2026)


Brief: .scratch/orca/ceilaudit6-brief.md. Work dir .scratch/ceilaudit6/ (code, logs). Lessons 191, 235, 309, 323, 329 applied:
each barrier gets (1) the step it fixes and its exact input class, (2) the freedoms of the method's own lemmas outside that class,
(3) a float test of the most promising freedom, (4) ROBUST/OPEN.

## 1. N5 — explore38 (integer-dilation penalty P(N) = d_N² − D(1/N)²)

### 1.1 Step fixed and class
- Step: the duality lower bound d_N² ≥ D(1/N)² ≥ B_R(N) = c_Rᴴ H_R⁻¹ c_R with test vectors Y_ρ ⊥ C_{1/N} (continuous span).
  Class: test vectors in C_{1/N}^⊥ (Christoffel–Darboux kernels at zeros, interior vectors at off-line zeros). Every such bound is
  ≤ D(1/N)², so "no kernel at the zeros controls P" is exact for this class: a tautology of the class, not a measured limit.
- Exact identity (orthogonal decomposition V_N ⊂ C_λ ⇒ V_N^⊥ = C_λ^⊥ ⊕ (C_λ ⊖ V_N)):
  P(N) = ‖P_{C_λ ⊖ V_N} 1‖² = dist(P_{C_λ}1, V_N)², i.e. the cost of replacing the continuous optimiser f·G* by f·(Dirichlet
  polynomial of length N). So P(N) is a second Nyman–Beurling problem with target P_{C}1 instead of 1.
- The scaling claims ("ζ, χ₋₄ second order; DH first order 0.020/log N") are measured on N ≤ 1280 (ζ), 600 (χ₋₄), 1000 (DH),
  not derived; explore59 §6 makes P(N) = o(1/log N) the whole discrete content of RH ⇒ d_N² ~ C/log N.

### 1.2 Freedoms outside the class
- (F1) test vectors in C_λ ⊖ V_N (components of non-integer dilations ρ(θ/x), θ ∉ 1/ℕ, orthogonal to V_N): the only vectors that see P.
- (F2) exact evaluation: Vasyunin's closed-form Gram (cotangent sums) and, for periodic coefficients, the digamma formula —
  this gives d_N exactly, hence the penalty's measured scaling at larger N (test below).
- (F3) extra mean values the cited papers prove: BCF (1211.5191) RH + J₋₁ ⇒ d_N² ~ C/log N, hence (with Burnol) P = o(1/log N)
  for ζ conditionally — so ζ's "second order" is consistent with a conditional theorem, and an Euler-product L-function with GRH
  (χ₋₄) must also have P = o(1/log N) if BCF extends (its V_N = Σ μχ(n)(1 − log n/log N)n^{−s} has bounded coefficients).
- (F4) several pieces: none needed; P is one quadratic form.
- Most promising falsification: (F2) at the largest N under the cap — does χ₋₄ (Euler product, GRH true numerically) show the
  same first-order gap as DH (then the separation is quantitative, not an order), and does DH's gap·log N stay at 0.020?

### 1.3 Test (exact d_N by the closed-form Grams; Burnol/Christoffel B as in explore38: K = 400 (ζ), 300 (χ₋₄, DH) zeros + diagonal tail)
- Code: .scratch/ceilaudit6/gram.c (C, double: Vasyunin for ζ, digamma breakpoint formula for period-q coefficients), dN.py (one in-place
  Cholesky gives every d_N, N ≤ NMAX), b_zeta.py / b_dh.py / b_chi.py (explore38's s3/s7/s9 with d_N read from file), ana1.py.
- Validation: ζ d²(1280) = 0.006361102637571 (nbrate certified …562); DH d²(400, 1000), χ₋₄ d²(100, 600) equal explore38's logs to
  8 digits; B·L(10, 100, 1280) = 0.026772, 0.032671, 0.036029 (explore38: 0.0268, 0.0326, 0.0360).
- Precision: entrywise relative noise ε on the Gram moves d² by ≤ 8·10⁻¹⁴ (ε = 10⁻¹³) and ≤ 9·10⁻¹² (ε = 10⁻¹¹) for ζ at every
  N ≤ 4000 (χ₋₄, N ≤ 2000: ≤ 2·10⁻¹⁵); optimal coefficients ‖x‖₂ ≤ 15.5, ‖x‖_∞ ≤ 1.1, min pivot² 2.8·10⁻⁷: well conditioned.
- New sizes: ζ to N = 4000, χ₋₄ to 3000, DH to 2800 (explore38: 1280, 600, 1000). δ = L²(d_N² − B), gap·L = L(d_N² − B):

| N window | ζ δ | ζ gap·L | χ₋₄ δ | χ₋₄ gap·L | DH δ | DH gap·L |
|---|---|---|---|---|---|---|
| 600–1000 | .0658 | .0101 | .0870 | .0133 | .1265 | .0194 |
| 1000–1500 | .0678 | .0095 | .0927 | .0131 | .1396 | .0197 |
| 1500–2000 | .0697 | .0094 | .0837 | .0112 | .1508 | .0202 |
| 2000–2500 | .0711 | .0092 | .0834 | .0108 | .1592 | .0207 |
| 2500–3000 | .0703 | .0089 | .0846 | .0107 | .1639 (≤ 2800) | .0208 |
| 3000–4000 | .0670 | .0082 | .0860 (3000) | .0107 | – | – |

- Fits on N ≥ 1000 (δ = a + bL | a + b/L, residual sd): ζ .0825 − .0018L | .0572 + .091/L (sd .0017 both); χ₋₄ .157 − .0094L |
  .015 + .54/L (.0028, .0027); DH −.0844 + .0316L (.00101) | .3853 − 1.745/L (.00108).
- Reading (δ = L²P + L²(D² − B), both parts ≥ 0, so δ is an upper bound for L²P; all "orders" below are finite-range plateaus).
  ζ: flat, δ = .068 ± .002 to N = 4000 (an oscillation, no trend) — explore38's "flat 0.064 ± 0.003" holds, at a level ≈ .004 higher;
  hence P_ζ ≤ .07/log²N on N ≤ 4000. χ₋₄: δ peaks at .093 (N ≈ 1300) and plateaus at .084–.086 on [1500, 3000], which settles
  explore38's open choice (0.040 + 0.007L vs 0.1134 − 0.189/L) in favour of a plateau on this range. DH: gap·L = .0194 → .0208,
  rising, no decay to N = 2800; the saturating law δ_∞ = .39 (gap·L peaking at L = 9.1, N ≈ 8600) fits as well as the linear one, so
  "first order" is not established; nor is it shown to sit in P rather than in D² − B (explore38 split it only at N ≤ 40, ≥ 89 % P).
  The separation at reachable N is quantitative: δ_DH/δ_ζ = 2.4 at N = 2800 and growing. Under these two fitted laws (a model-
  dependent forecast; residual sd is not a predictive error) they differ by ≥ .003 from N ≈ 4300–4500 (Gram 150–165 MB; ≈ 30 min
  single-threaded with this O(N³) breakpoint code; a cotangent-sum reciprocity for the periodic Gram would make it cheap).
- External check: C₂(2560) = L²d² − CL = +.0032 here (exact to 10⁻¹¹ by the noise test) vs 2611.03924's V2 +.0028 and NB2nd's
  inferred V1 ≈ +.0025: the inferred V2 − V1 ≈ +.0003 (or the transcribed value) is off by ≈ .0007; flag for NB2nd.
- Caveat: B uses explore38's zero truncation (K = 400/300 + diagonal tail); a K-doubling check at the largest N is in §1.5.

### 1.4 Verdict N5: ROBUST (exact class), scaling claim rescoped
- "No kernel at the zeros controls P" is exact for test vectors in C_{1/N}^⊥ (they are orthogonal to P_C1, so they give 0 for P).
- Scaling: ζ (δ ≈ .068 flat, N ≤ 4000) and χ₋₄ (δ plateau ≈ .085, N ≤ 3000) keep bounded δ at 3–5× the old sizes, i.e.
  P = O(1/log²N) on the range (an upper bound, since D² − B ≥ 0); DH's "first order 0.020/log N" is confirmed only as a rising
  gap·L ≈ .02 to N = 2800: its order is undetermined (a plateau δ_∞ ≈ .39 fits equally), and its split between P and D² − B is
  unmeasured above N = 40. explore59's reduction (RH ⇒ d_N² ~ C/log N ⟺ P = o(1/log N), given its Lemma R′) stands: ζ's data
  satisfy P ≤ .07/log²N to 4000; DH's do not decide.
- Reviewer (A1 identity / exact class): CORRECT; (A2 scaling reading): GAP — "d² − B = P + (D² − B), plateaus are finite-range,
  the forecast is model-dependent": accepted and written in above.


### 1.5 K-doubling check (logs/run5.log): δ changes by ≤ .0007 at the largest N (ζ K 400→800 at 4000: .06735→.06723; DH 300→600 at
  2800: .1665→.1658; χ₋₄ 300→600 at 3000: .0860→.0856). Conclusions unchanged.
- Process note (cap): one dN_check run (ζ, N = 4000) reached maxrss 370 MB for ≈ 8 s (a 128 MB slice copy in solve_triangular),
  over the 300 MB cap. memcap did not catch it: it under-reports `uv run` children (a 150 MB test: 18 MB reported under uv, 157 MB
  for plain python3) and misses runs shorter than its first top sample. Runs 3–8 used /usr/bin/time -l only. All other peaks
  ≤ 248 MB maxrss (logs/run*.log, run8.time). Fix suggested for memcap: follow uv's python child (pgrep -f) or cap via ulimit.

## 2. N13 — aclgh2 (onset of ACLG's (H2ii) localisation)

### 2.1 Step fixed and class (ACLG TeX prolate-weil-fine15.tex, arXiv 2611.13129)
- Theorem condpos (l. 1318–1330): W_λ(σ, n₀) has codimension ≤ 2 + 2(n₀+1)|M|, M = {m : κ_m(λσ) > κ_{n₀}(σ)}; n₀ is fixed by
  ε_σ(n₀) ≤ min(1, log σ/(C₀λ)); the λ comes from Lemma archdom (l. 1309) Λ(f) ≤ QW_λ(f,f) + C₀λ‖f‖², Λ(f) = (1/2π)∫|Mf(is)|²
  log(2+|s|)ds (l. 1223), i.e. the sup bound of the prime sum Σ_{p≤λ²}W_p. (H2′)/(H2ii) cut the pair at an angle index n₀.
- Barrier: n₀ ≥ N_S(κ♮) ≥ Σ_{κ>κ♮}κ² ≈ 2F_S with κ♮ ≈ τ/A_R ~ 1/λ, resting on the semiclassical law κ² ≈ A/2π per coupling pocket.
  Class: conditioned spaces built from angle-ordered Halmos blocks of the semilocal pair with a threshold κ♮ → 0 as λ → ∞.

### 2.2 Freedoms outside the class
- (a) Any finite block set F instead of {n ≤ n₀}: Prop. pairexp is a sum over orthogonal Halmos blocks, so it holds for any F; but
  (H2′) for the excluded blocks needs Σ_{n∉F}κ_n‖θ(f)a_n‖² ≤ εΛ(f), and Λ's weight is only log(2+|s|): blocks with κ_n ≳ ε log(2+τ_n)
  must stay in F. Gain: log factors only.
- (b) Weighted κ: the codimension counts blocks; Σκ² is already the barrier's lower bound for the count.
- (c) Different S: forced; QW_λ = W_{0,2} − W_R − Σ_{p≤λ²}W_p (l. 1315).
- (d) Different absorption inequality (e.g. Λ ≤ QW + C log λ‖f‖²): κ♮ still → 0 and every pocket still counts (test below); only a
  λ-independent threshold κ_* would drop pockets, which needs Λ(f) ≤ C·QW + O(‖f‖²), stronger than the positivity being proved.
- Most promising falsification: the semiclassical law itself on THIN pockets (the bulk: median A/2π ≈ .02 for |S| = 6), which
  ACLGH2 validated only on thick ones (A/2π = .21–.44, τ ≤ 2·10⁴).

### 2.3 Test (ACLGH2's band-limited window model, code/window_count.py, tmax 100, depth 16; pockets from a dump of pockets.py)
| S | τ_c | A/2π | √(A/2π) | κ_top | # κ > .01 |
|---|---|---|---|---|---|
| {∞,2,…,11} | 1 046 476.56 | 2.1·10⁻⁵ | .0045 | .1645 | 4 |
| | 1 026 625.04 | 2.5·10⁻⁴ | .016 | .1644 | 7 |
| | 1 044 028.94 | 1.0·10⁻³ | .032 | .1689 | 6 |
| | 1 040 929.10 | .0186 | .136 | .2010 | 4 |
| | 3 077 708.30 | 3.3·10⁻⁵ | .0058 | .1614 | 4 |
| | no pocket: 1 033 777, 1 042 479, 10⁶, 1.5·10⁶, 4.5·10⁶ (> horizon 4.1·10⁶) | 0 | 0 | .038, .083, .021, .042, .031 | 8, 8, 4, –, – |
| {∞,2,3,5,7} | 35 053, 28 608, 26 125, 41 979 | 1.7–5.3·10⁻⁴ | .013–.023 | .151–.171 | 19–24 |
- Stability (logs/run4.log): at τ_c = 1 046 476.56, κ_top = .1645 for tmax 50/100/200, depth 16/32, centre ±40: not an edge artefact.
- Quasimode check (wc_vec.py: top singular vectors, spectral profile; logs/run6.log, run7.log): at every tested thin pocket the
  κ ≈ .16 mode is pinned to the pocket's global position (windows shifted by ±40 move its spectral peak to τ′ = ∓40.0), holds 98 % of
  its spectral mass within ±10 of the pocket and leaks 1.1–1.8·10⁻³ beyond τ₁ (where window and global symbols differ). Since the
  window symbol equals the shifted global one on |τ′| ≤ τ₁, the global operator has an approximate singular pair with residual
  ≤ 2(leakage)^{1/2} ≈ .07, i.e. a global angle ≳ .09 per pocket; pockets ≥ 150 apart give nearly orthogonal quasimodes. Some
  window-top angles are closure artefacts (spectral mass > 80 % beyond τ₁, e.g. κ = .024 at τ_c = 1.5·10⁶); the pocket modes are not.
- One operator, many pockets (run8.log; no gluing): S = {∞,2,3,5,7}, τ_c = 23 150, τ₁ = 1200 (tmax 1600, depth 6; depth 6 vs 16
  changes κ_top by .0009), 8 pockets inside (A/2π = 1.4·10⁻³ … 7.2·10⁻², Σ A/2π = .282 = N_sc): Σκ² = 1.11, 8 angles ≥ .156
  (.266, .262, .259, .232, .215, .165, .161, .156), 27 ≥ .1, 70 ≥ .05, 377 ≥ .01 (248 MB maxrss).
- Reading: on the 8 thin pockets tested (A/2π ≤ 10⁻³; 5 with the quasimode check), κ does not follow √(A/2π); it stays near .16 as A → 0 (near-resonance, Landau–Zener-like),
  and genuine pocket-free modes (inside τ₁) carry κ ≈ .02–.08 > κ♮ (≤ .007). If typical, every pocket counts: N_S(κ♮) ≳ #pockets
  (970 for |S| = 6, against 2F_S = 97) + background.

### 2.4 Verdict N13: ROBUST against (a)–(d); its semiclassical mechanism is contradicted in the model
- None of (a)–(d) lowers the count (structural, from the TeX). The stated exponent exp((½ + o(1)) log² λ²) rests on κ² ≈ A/2π,
  which the window model contradicts on thin pockets (κ ≈ .16 at A/2π = 2·10⁻⁵) and in bulk (one 8-pocket operator: Σκ² = 1.11
  against Σ A/2π = .28). So the semiclassical count is not the growth law; the window data suggest a larger onset, but global pocket
  counting is a numerical conjecture: it needs well-conditioned Gram matrices of the input and output quasimodes and a collective
  residual bound ‖Σ_j a_j(Γv_j − κ_ju_j)‖ ≤ ε‖a‖ (individual residuals .066–.085 give only ‖Γv_j‖ ≳ .076 each).
- Reviewer, pass 1 (B): GAP — "window-top angles need not be pocket modes; separate windows need not give independent global modes".
  Pass 2 (B2, with the quasimode data): GAP — "individual residuals plus near-orthogonal inputs do not give one global singular value
  per pocket; a numerical conjecture with an identified proof obligation". Both accepted; the text above states only that.


## 3. N12 (N5, N1) — explore48 ("Euler products (h = 1) have T_s ≤ π√163, so geometry says nothing about ζ above height 40")

### 3.1 Step fixed and class
- Step: Lemma K (rigorous Bessel bound |K_{s−½}(u)| small for |t| < u) on the dual-frequency tail of the Chowla–Selberg expansion
  along the shortest vector; the Proposition (sketch) then gives no off-line zeros below (1 − δ)T_s, T_s = 2π/λ₁² (y ≥ y₀(δ)).
  The height is fixed by u_min = 2π·μ₁(L′)/λ₁ (shortest dual frequency), and ζ enters only via E_Q = ζ·L(s, χ_D), i.e. h(D) = 1.
- Class: single Epstein zetas with absolute-value (phase-blind) tail bounds; ζ as an analytic factor.

### 3.2 Freedoms and tests (.scratch/ceilaudit6/lat.py, logs/run3.log)
- (a) Higher rank. Rank-n Chowla–Selberg has T_s = 2π μ₁(L′)/λ₁ and can be large for skewed lattices. Test: all 158 diagonal forms
  x² + a₂y² + a₃z² + a₄w² (a_i ≤ 24, a₂a₃a₄ a square: trivial character), r_Q(n) for n ≤ 3000; ζ(s)ζ(s−1) is a factor iff
  r_Q(mk) = f(m)σ₁(k) for k prime to the level. Only (1,1,1,1), (1,1,2,2), (1,1,3,3), (1,2,2,4) pass, with T_s = 2π√a₂ = 6.3, 6.3,
  6.3, 8.9; none of the 58 forms with T_s > γ₁ = 14.13 passes (r/σ relative sd ≥ 3 %). Structural reason: in rank n ≥ 3 the Eisenstein
  factors L(s, ψ)L(s − n/2 + 1, φ) have their zeros on Re s = ½ and n/2 − ½, never on the Epstein line Re s = n/4; a rank-n analogue
  of the Proposition would therefore force T_s < γ₁/(1 − δ) for every lattice with ζ | E_L. Higher rank cannot reach ζ's zeros.
- (b) Rank 3 puts ζ's zeros on its line: Shimura gives E_{Z³}(s) = ζ(2s)ζ(2s − 1)·Σ_t r₃(t)t^{−s}/L(2s, χ_{−4t}), zeros of ζ(2s − 1) at
  (1 + ρ)/2 on Re s = ¾, height γ/2. Only formal: E_{Z³}(¾ + iγ₁/2) = 8.10 + 11.17i (theta integral, 30 digits; E_{Z³}(2) = 16.5323
  agrees with the partial sum + tail), so the quotient has compensating poles and ζ(2s − 1) is no analytic factor.
- (c) Hecke characters / class-group sums L_K(s, ψ) = Σ_C ψ(C)E_C(s)/w: Lemma K controls class C only below T_s(C) = π√|D|/a_C,
  and for h ≥ 2 some class has a_C ≥ 2; keeping the phases ψ(C) turns the summed tails into the Dirichlet coefficients of L_K
  (divisor sums with χ_D): arithmetic, not geometry. (d) Exact evaluation of the tail instead of Lemma K: same — the exact tail is
  the arithmetic part.

### 3.3 Verdict N12: ROBUST
- Exact class: Chowla–Selberg expansions of one Epstein zeta of any rank with Lemma-K absolute bounds on the dual-frequency tails;
  ζ appears as an analytic factor only in rank 2 (h(D) = 1, T_s ≤ π√163) and, in rank ≥ 3, only for lattices with T_s < γ₁ (found:
  T_s ≤ 8.9); rank 3's ζ(2s − 1) is formal. Keeping phases across classes is the arithmetic of L_K.


## 4. N1, N13 — explore32 ("step one carries no RH information"; "no prime-blind bound passes c = 2") — reasoning, no float test

### 4.1 Steps fixed and classes
- Claim 1 (statement-level): CCM's step one = "the smallest eigenvalue of the window form is simple, even, nodeless". It holds for
  ζ, L(s, χ₅), DH (also at its failing window c = 31: λ_min = −6.08·10⁻³⁰ even / 7.46·10⁻²⁶ odd, arb at ≥ 128 bits) and 14 moved-pair
  models on every positive window (84/84). Class: qualitative ground-state properties (simplicity, parity, sign) on positive windows.
- Claim 2 (method-level): window positivity at c = 2 has margin 1.3·10⁻³ against the p = 2 weight w₂ = 0.49. Class: bounds that use
  the gamma factor and conductor only (Burnol's conductor operator, Beckner's ψ(¼) − ln π, Hirschman–Beckner, Maassen–Uffink) or a
  pointwise kernel sign (Perron–Frobenius: κ(z) ≥ 2cosh(z/2) ⟺ λ² ≤ ρ = 1.3247).

### 4.2 Freedoms outside the classes
- (a) A proof of step one that uses ζ's Markov sign (prime kernel −Λ(n)n^{−½}δ(x − y ± log n) ≤ 0, which helps Perron–Frobenius):
  RH content is a property of the statement, not of the proof; DH satisfies the statement without the sign, so Claim 1 is
  unaffected. And the prime kernel sits at |x − y| = log n ≥ log 2 > log ρ, so it cannot repair the archimedean kernel's sign on
  (log ρ, log 2): pointwise PF stops at λ² = ρ for every L-function with factor Γ(s/2).
- (b) Exact evaluation of finitely many prime terms (the p = 2 shift is an explicit operator): this is the computer-assisted window
  route (Botha L ≤ 2, Kramarenko–Byrd); records, not a node move; outside Claim 2's class by definition.
- (c) Quantitative ground-state data (λ_min itself; DH's 4.5·10⁻¹⁸ at c = 20 vs L(s, χ₅)'s 8.6·10⁻⁶): RH-sensitive, but that is
  positivity, i.e. not step one. (d) Parity after failure: explore32's own log has an odd ground state only after failure
  (moved pair k = 5, δ = 0.3, c = 3.53), so parity flags failure late, never on a positive window.

### 4.3 Verdict: ROBUST (both claims, in the classes above). Skipped a float test: nothing in (a)–(d) is a step-one bound that could
  separate ζ from DH on a positive window.

## 5. N10 — explore33 ("every Boolean statistic within 3.1 % of a planted-zero Beurling system") — reasoning, no float test

### 5.1 Step fixed and class
- Encoding: f_x(S) = 1 iff ∏_{p∈S} p ≤ x on {0,1}^n, n = π(x); μ(f) = Q(x)/2ⁿ, top coefficient f̂([n]) = 2⁻ⁿM(x).
- The cited theorems (level-k inequality, hypercontractivity, KKL, Friedgut, Peres, noise sensitivity at fixed ρ) bound level weights
  W^k = Σ_{|T|=k} f̂(T)² and influences in terms of μ; none is sharper than Parseval (Σ W^k = μ) at the top level.
  The top weight is (2⁻ⁿM)²; relative to Parseval's total it is M²/(Q·2ⁿ) = 2^{−n(1−o(1))}.
  Class: statistics Σ_k w_k W^k and influence functionals whose top-level weight is 2^{o(n)} relative to the total μ.

### 5.2 Freedoms
- (a) Fourier weight at high levels W^{≥n−j} at its own scale 4⁻ⁿ: these are the sign-twisted squarefree sums (Möbius except on j
  primes), so they are RH-sensitive, but no cited theorem bounds them below Parseval: the level-k inequality
  W^{≤k} ≤ μ²(2e ln(1/μ)/k)^k reaches k = n (2 ln(1/μ) ≈ 1.39n) yet gives μ²(3.77)ⁿ ≫ μ there.
- (b) Noise sensitivity at ρ → −1 or at measure 2^{−n(1−o(1))}: ⟨f, T_ρ f⟩ = Σ ρ^{|T|}f̂(T)² weights the top by |ρ|ⁿ ≤ 1 among a total
  Q/2ⁿ: invisible at relative scale 2⁻ⁿ. (c) The top coefficient alone is M(x): circular.

### 5.3 Verdict: ROBUST in the class of 5.1. Skipped a float test: freedom (a) is a restatement of M-type sums, not a Boolean bound.


## Suggested NOGO lines
- ceilaudit6 (Orca, N5 audit of explore38): "no kernel at the zeros controls the integer-dilation penalty" is exact for test vectors
  ⊥ C_λ (P = dist(P_C1, V_N)², and such vectors are ⊥ P_C1); exact Grams to N = 4000 (ζ), 3000 (χ₋₄), 2800 (DH): ζ's δ = L²(d_N² − B)
  = .068 ± .002 flat, χ₋₄'s saturates at .085 (second order), DH's gap·log N rises .0194 → .0208, first vs second order (δ_∞ ≈ .39)
  undecidable below N ≈ 4500 [CeilAudit6](Riemann/Deep/CeilAudit6/REPORT.md)
- ceilaudit6 (N13 audit of aclgh2): the onset barrier survives every freedom of Theorem condpos (any finite Halmos block set, any
  absorption with κ♮ → 0, weighted counts; S ⊇ {p ≤ λ²} forced; only a λ-free threshold, i.e. Λ ≤ C·QW + O(‖f‖²), would drop
  pockets), but its law κ² ≈ A/2π fails in ACLG's window model (thin pockets κ ≈ .16 at A/2π = 2·10⁻⁵; one 8-pocket window
  Σκ² = 1.11 vs Σ A/2π = .28): exp(½ log²λ²) is not the growth law; a larger global count is unproved (needs collective quasimode
  residuals) [CeilAudit6]
- ceilaudit6 (N12 audit of explore48): the height-40 barrier covers every rank: ζ is an analytic factor of E_L only for h(D) = 1 in
  rank 2 and, in rank ≥ 3, only for lattices with T_s < γ₁ (Eisenstein factors sit off Re s = n/4; of 158 diagonal quaternary forms
  4 have the factor, T_s ≤ 8.9); rank 3's ζ(2s − 1) is formal (E_{Z³}(¾ + iγ₁/2) = 8.10 + 11.17i); class sums with phases are the
  arithmetic of L_K [CeilAudit6]
- ceilaudit6 (N1, N13 audit of explore32): step one is RH-blind at the statement level (simple even nodeless ground state on positive
  windows for DH and moved pairs, whatever the proof); "no prime-blind bound passes c = 2" is exact for gamma-factor-only bounds and
  pointwise kernel-sign Perron–Frobenius, which stops at λ² = ρ for every Γ(s/2) L-function (prime kernels sit at |x − y| ≥ log 2);
  exact prime terms are the computer-assisted window records [CeilAudit6]
- ceilaudit6 (N10 audit of explore33): the Boolean barrier is exact for level-weight and influence statistics with 2^{o(n)} relative
  top weight (the top weight is 2^{−n(1−o(1))} of Parseval); high-level weights at scale 4⁻ⁿ are M-type sums that no cited
  inequality bounds below Parseval (level-k gives μ²·3.77ⁿ at k = n) [CeilAudit6]

RH relevance: N5 (load-bearing for explore59's reduction: ζ's P(N) = O(1/log²N) holds in floats to N = 4000; the penalty separates
ζ, χ₋₄ from DH only quantitatively at reachable N), N13 (ACLG route: the codimension obstruction stands and is larger than stated),
N12, N1, N10 (barriers stand with exact classes).

## Files
- This folder: REPORT.md; code/ (gram.c, dN.py, dN_check.py, b_zeta.py, b_dh.py, b_chi.py, b_*_k.py, ana1.py, lat.py, pk_dump.py,
  pk_select.py, pk_in.py, wc_vec.py, run1–8.sh); logs/ (run1–8.log, run8.time, ana1.txt, pocket dumps). Work copies in .scratch/ceilaudit6/.
- Reviewer sessions: pass 1 (A1 CORRECT, A2 GAP, B GAP) and pass 2 (B2 GAP), both accepted in §1.4 and §2.4.
- d_N arrays (d2z4000.npy, d2c3000.npy, d2dh2800.npy) in logs/ (exact d_N² for every N ≤ NMAX; 32 KB each).
- Inputs read: explore38 (ExploreChristoffel code s3/s7/s9), ExplorePotential DH zeros, WangAudit zeros1250, ExploreHelio zeros,
  ACLGH2 code (window_count.py, pockets.py), ACLG TeX (arXiv 2611.13129, deleted after use).

