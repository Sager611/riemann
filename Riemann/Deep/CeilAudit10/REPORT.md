OPEN: e1rate, freedom 'keep the Bessel phase on the Deshouillers–Iwaniec type-II terms' (exponent bookkeeping in e1rate's own exponents.py, conditional on a phase-kept Kuznetsov bound nobody has): ν < (2/7)(3θ−1) = (8/7)·E1 for θ ≤ 0.8, a loss up to X0^{1/8} tolerated at θ = 0.6, short-window LB 0.1408 → 0.1550 at θ = 0.6 (gain only for θ < 0.6176). ROBUST: explore63 — C1 fails for every λ tested (floats to 10⁵), every damping is capped by a sharp Newton budget (Lemma D), the continuous-Hahn route dies the same way, and Jensen's optimality extends to every preserver faithful on one exponential (Lemma M′); gl2ceil — the certified one-piece cap covers optimiser, degree, R and Feng's first piece (μ_f∗Λ_f = −μ_f·log); skipped: the Feng k = 2 two-piece float test (no implementable GL(2) main term).

# ceilaudit10 — falsification of OUR barriers explore63 (N6), e1rate (N8), gl2ceil (N8/N12)

Status: started 20 Nov 2026 15:32 UTC; order of work: explore63 (best chance) → e1rate (exponent bookkeeping) → gl2ceil.
Logs: Riemann/Deep/CeilAudit10/code/ (scripts and logs); every run under bash .scratch/memcap.sh 290 <s>.

## 1. explore63 (N6) — orthogonal-basis Pólya–Turán transfers; Lemma M
### 1.1 Steps fixed and input classes
- C1 fixes the pre-image test: G_N(z) = Σ_{n≤N} (−1)^n r^{2n} c_{2n}^{(λ)} β_n^{(λ)} z^{2n}, β_n = (2λ)_{2n}/(2n)!, must be real-rooted (Romik's lemma
  for MP_λ); Newton at k = 1 needs R(λ) := c₂²/(c₀c₄) ≥ ν(λ)·2N/(N−1), ν = (2λ+2)(2λ+3)/(12λ(2λ+1)). Class: every positive Φ, fixed λ, dilation,
  Poisson factor r, undamped truncation. Proof for λ < λ* = (1+√13)/4 = 1.1514 (2ν > 1 ≥ R); ζ grid λ = 1.25…30; asymptotic λ → ∞.
- C3 fixes the approximant itself (undamped truncation in MP/Hermite/continuous-Hahn bases): non-real zeros (Arb) — an observation, no preserver.
- Lemma M fixes the damping of a degree-d DIAGONAL test (finite multiplier sequence, γ₀ = γ₁ = 1): γ_k ≤ (d)_k/d^k (Maclaurin). Class: diagonal
  operators x^k ↦ γ_k x^k at one centre. Reading "reaches no further than √d" is an argument (damping ⇒ reach), not a theorem.
### 1.2 Freedoms the cited lemmas admit but the classes exclude
- F1 damping μ_n of the pre-image (Fejér/Jackson/Jensen/any) — C1 is stated for undamped truncations only.
- F2 λ = λ(N) growing with N (C1 is per fixed λ; the gaps (1.1514, 1.25) and (30, ∞) are not certified for ζ).
- F3 damped truncations of the approximant (Fejér, Jackson) — C3 tested only undamped ones.
- F4 non-diagonal real-rootedness preservers (Borcea–Brändén symbols), translation of the centre, the t-variable vs the x = −t² variable.
### 1.3 Tests (scripts and logs in Riemann/Deep/CeilAudit10/code/; peaks ≤ 59 MB)
- T1 λ-gaps of C1 (lam_gap.py/.log; ζ, mpmath dps 30, adaptive Gauss–Legendre). R(λ) < 2ν(λ) everywhere tested:
  [1.10, 1.30], 26 points: R = 0.3842…0.3832 against 2ν = 1.034…0.917 (relative margin ≥ 58 %; the "gap" (1.1514, 1.25) is not tight).
  [25, 10⁵], 60 points (ratio 1.15): λ(2ν − R) = 0.344, 0.325 (min, λ ≈ 58), 0.326, 0.329, …, 0.33332 → 1/3 (the λ → ∞ law); min relative
  margin 1.05·10⁻⁵ at λ = 9.5·10⁴ (that is (1/3)/λ ÷ 1/3, as predicted). Reproduces explore63 at λ = 1.25 (0.383437) and 30 (0.344524).
  Consequence: Newton k = 1 fails for every λ on the grid, every N ≥ 2, undamped AND Jensen-damped (the damped test needs R ≥ 2ν exactly).
- Lemma D (proof; every positive Φ, every λ > 0, N, r, dilation): if a damping μ ∈ (0,∞)^{N+1} makes the MP_λ pre-image
  Σ_{n≤N} (−1)^n μ_n c_{2n}β_n z^{2n} real-rooted, then for 2 ≤ m ≤ N
    μ_m/μ_{m−1} ≤ (μ₁/μ₀) · 2λ(2λ+1)(2m−1)/((2λ+2m−2)(2λ+2m−1)) · (N−m+1)/N  (≈ (μ₁/μ₀)λ(2λ+1)/m for m ≫ λ).
  Proof: Newton at n = 1…m−1 for the degree-N polynomial in w = z² with positive roots, minus the log-convexity of the Hausdorff moments c_{2n}
  (c² ≤ c₋c₊), summed: second differences of log μ telescope to Δ₁ − Δ_m, those of log β to log β₁/β₀ − log β_m/β_{m−1}, and
  Σ log(1+1/n) + Σ log(1+1/(N−n)) = log m + log(N/(N−m+1)). So every real-rooted damped pre-image keeps ≲ e·λ(2λ+1)·μ₁/μ₀ effective terms
  (λ = 3/4, μ₀ = μ₁: μ₃ ≤ 0.375, μ₅ ≤ 0.0522, μ₁₀ ≤ 2.9·10⁻⁵, every N), while seeing γ₁ needs ≈ 0.049·14.13² ≈ 10 MP terms: damping (Fejér, Jackson, Jensen, any) cannot rescue the route at fixed λ.
  λ = λ(N) needs λ ≳ √(N/2); the undamped MP_λ terms do not decay for n ≲ λ (Gaussian regime: c_{2n}|P_{2n}| ≍ n^{−1/2}) and decay as
  e^{−4√(πn)} beyond, so the best case is N ≍ max(λ², 0.049t²): a tie with the T² law, never better (heuristic exponent count).
- T2a damped MP truncations of Ξ itself (damped_mp.py/.log; Arb, λ = 3/4, c_{2n} from explore63): non-real zeros
  none: from N = 2 (4/4); Fejér: from N = 2; Jackson: N = 2, 5 real, from N = 10 (6/20; 22/80 at N = 40); Jensen μ = (N)_n/N^n: real for
  N = 2, 5, 10, non-real from N = 20 (4/40; N = 40: 4/80); Gaussian e^{−n²/N}: from N = 5. Damping also shifts the real zeros (Jensen N = 40:
  13.565, 21.087 for 14.135, 21.022). No damping keeps the approximants real-rooted, and none is certified by a preserver.
- T2b NEW FREEDOM — a Pólya–Turán–Romik lemma for the continuous-Hahn basis (hahn_lemma.py, hahn_lemma2.py/.log): z^k ↦ κ_k g_k, with
  explore63's g_n ((2n+3)g_{n+1} = 8x g_n − (2n+1)g_{n−1}), maps real-rooted p (degree 2–10, random roots, scales 0.1–30, 300 trials each,
  all-sign and negative roots) to real-rooted images for κ_k = 1 (Romik's normalisation), k!^{1/4}, k!^{1/2}, (3/2)_k^{1/2}; it fails for
  κ_k = k!^{3/4} (6/300, 35/300), k! (54, 73), (3/2)_k = monic (11, 42). explore63 covered Hahn only by C3; its C1 (MP) did not reach this basis.
- T2c the Hahn pre-image (hahn_pre.py/.log): Romik's d_{2n} (n ≤ 100) are log-convex (d²/(d₋d₊) = 0.522, 0.715, 0.801, …, 0.975 (n = 20),
  0.9968 (n = 99) < 1), so with κ = 1 Newton k = 1 reads 0.522 ≥ 2N/(N−1): false; Arb: the pre-images q_N(w) have only non-real roots for
  N = 2…40 under all four κ (1, k!/(3/2)_k, 4^k k!/(3/2)_k², 1/(3/2)_k). With the largest passing κ ≈ (3/2)_k^{1/2} the basis supplies
  log-concavity (1 + 1/n + O(n⁻²)), the log-convexity deficit of d is ≈ 1.4 n^{−4/3} (n = 20–50), and Lemma D (general form, §1.4,
  with b_n = 1/κ_{2n} and the measured ρ_n = d²/(d₋d₊)) gives μ_m/μ_{m−1} ≤ 0.535, 0.384, 0.308, 0.262 (m = 2…5), 0.165 (10), 0.058 (99) × μ₁/μ₀
  (budget.py/.log): μ₅ ≤ 0.0166 for μ₀ = μ₁, while seeing γ₁ needs ≈ 0.2·14.13^{3/2} ≈ 11 Hahn terms.
  Only κ_k ≍ k!^s with s > 1/2 (lemma fails in floats at s = 3/4; threshold in (1/2, 3/4) unknown) would give the basis enough log-concavity.
- Lemma M′ (proof; non-diagonal, every preserver): if P is real-rooted of degree d, positive on I = [x₀−X, x₀+X], and |log P − (ax+b)| ≤ ε on I,
  then d ≥ (aX − 4ε)₊²/(16ε). Proof: h = log P − ax − b is concave (h'' = −Σ(x−x_j)⁻²), |h'| ≤ 4ε/X on the middle half J, Cauchy–Schwarz
  (Σ(x−x_j)⁻¹)² ≤ dΣ(x−x_j)⁻² gives h'' ≤ −(a − 4ε/X)²/d on J, and a concave h with |h| ≤ ε has h'' ≥ −16ε/X² somewhere on J.
  Hence every linear real-rootedness preserver T (any Borcea–Brändén symbol, any centre) that is ε-close to the identity on one exponential
  e^{ax} ∈ LP⁺ over a window with aX = A has degree ≥ (A − 4ε)²/(16ε) (Jensen: A²/(2ε)); reproducing R_n[x] on [0, X] forces A ≍ n, so
  d ≳ n²: Lemma M's √d reach is optimal up to a constant in this non-diagonal class. Not covered: preservers that are NOT close to the identity on
  exponentials (rank ≤ 2 Hermite–Kakeya maps are real-rooted for every input, hence vacuous as tests; no non-vacuous candidate found). Applied to
  F(x) = Ξ(√−x) itself, Cauchy–Schwarz gives only d ≳ (log F)'²/|(log F)''| ≍ t log t/2 (the truncated product is real-rooted): the T² law is a
  property of universal (input-independent) preservers, not of real-rooted approximation of Ξ.

### 1.4 Lemma D in general form; verdict
- Lemma D (general; reviewer: symbolic statement and proof CORRECT, sharp). Basis with a Pólya-type lemma z^k ↦ κ_k p_k, expansion coefficients
  e_n log-convex (ρ_n := e_n²/(e_{n−1}e_{n+1}) ≤ 1), pre-image a_n = μ_n e_n b_n, b_n = 1/κ_{2n} (MP_λ: b_n = (2λ)_{2n}/(2n)!). If Σ_{n≤N} a_n(−w)^n
  has N positive roots, then μ_m/μ_{m−1} ≤ (μ₁/μ₀)(b₁/b₀)(b_{m−1}/b_m)·Π_{n<m}ρ_n·(N−m+1)/(Nm). Attained for e ≡ 1, q = (1 + cv/N)^N (reviewer).
  The factor (N−m+1)/N is Jensen's (N)_m/N^m; the 1/m must be paid by the basis. Taylor basis (Pólya–Jensen): Ξ's γ(k) are log-concave (Turán
  inequalities), only Jensen's factor remains (√N faithful terms). MP_λ: b supplies ≈ 1 + (2λ−1)/m, so ≲ e·λ(2λ+1) terms; Hahn: ≲ 3 terms.
  Reviewer's GAP was rounding only: μ₅ ≤ 0.0521924 (not 0.052) at λ = 3/4; corrected above.
- Verdict explore63: ROBUST, class extended. C1 holds for every tested λ (floats on [1.10, 1.30] ∪ [25, 10⁵], λ(2ν − R) → 1/3), undamped and
  Jensen-damped; every damping is capped by Lemma D (proof, every positive Φ); the continuous-Hahn route (new Pólya lemma in floats) dies the same
  way (log-convex d_{2n}); damped truncations are non-real from N ≤ 20; Lemma M′ extends "Jensen is optimal" from diagonal multipliers to every
  preserver ε-faithful on one exponential. Stated too strongly in explore63: "so N6's T² law stands" was an inference from coefficient
  domination; it is now a theorem only for the ε-faithful class (Lemma M′). Not covered: preservers not close to the identity on exponentials,
  bases (other than Taylor) in which Ξ's coefficients are log-concave.

## 2. e1rate (N8) — the separated Bessel phase
### 2.1 Step fixed and class
- DERIVATION §2: Mellin separation of Φ(x) = e(−x)𝒲(x) in the single variable x = n/(pq), cost Λ = √X0 ≍ (T/H)^{1/2} at N ≍ N0, then a θ = 1
  bound uniform in the Mellin variable (DI Lemma 1/Thm 12, BC trilinear, Weil type I with conductor pX0), arbitrary outer coefficients; §8 transfer
  principle ⇒ ν < (2θ − 1)ν₁; beats E1's (3θ − 1)/4 iff θ > (4ν₁ − 1)/(8ν₁ − 3).
### 2.2 Freedoms (cited lemmas admit, class excludes)
- (a) keep the t-average |t| ≍ X0 inside Cauchy–Schwarz/Kuznetsov (hybrid Kloosterman–archimedean large sieve); (b) type I sums averaged over n
  and the second mollifier variable (Lemma I numerics: |S| ≈ √D); (c) smoother or band-limited window: X0 at N ≍ N0 is fixed by H (uncertainty),
  smoothness only shortens N > N0 (E1 Step 4 already uses (N0/N)^{K/2}) — no exponent moves; (d) stationary phase/Poisson in n inside the
  Kloosterman sum: dual length ≍ p²X0²/N0 ≍ TP/Q, the pre-Voronoi shifted convolution (involutive) — no exponent moves; (e) lift SL's cap ν < 1/2.
### 2.3 Test: exponent bookkeeping against h(θ) (e1_book.py → e1_book.log, e1_book2.log; e1_book_exact.py/.log; peaks ≤ 46 MB)
- Exact criterion (reviewer CORRECT): h(θ;κ) = (3+κ−√((1−κ)(9−κ−8c)))/4 has h(θ; c) = c and ∂h/∂κ > 0, so for c(θ) ≥ 0
  LB = max{0, c, (c+2κ)/3, h} > c(θ) ⟺ κ(ν) > c(θ). e1rate's dominance check (numeric, κ = 0.409) is an identity.
- κ(ν) (one-piece SL class; kappa_opt.opt deg 4 = lower bounds; = e1rate's deg 8 within 10⁻⁵ for ν ≥ 0.36, up to 0.009 below for ν ≤ 0.25):
  0.1426 (0.2), 0.2186 (0.3), 0.2936 (0.4), 0.36350 (1/2), 0.40806 (4/7), 0.4968 (3/4), 0.5736 (1).
- c(θ) = κ(1/2) at θ = 0.71399 (c(5/7) = 0.36398); = κ(4/7) at 0.74252; = κ(1) at 0.87936. Needed ν*(θ) = κ⁻¹(c(θ)): 0.190 (θ = 0.6), 0.266 (0.625),
  0.337 (0.65), 0.402 (0.675), 0.466 (0.7), 0.528 (0.725), 0.592 (0.75), 0.804 (0.825), 0.98 (0.875), none ≤ 1 for θ ≥ 0.88.
- Route (a) in e1rate's own bookkeeping (e1_routes.py → e1_routes_a.log: exponents.py imported read-only, only eII patched to etr + gDI, i.e. the DI
  type-II terms keep the Bessel phase at no X0 cost; type I by Lemma I, the trivial bound, both orientations and the Vaughan-cut minimax unchanged —
  reviewer: patch correct): ν_K = (2/7)(3θ−1) = (8/7)·E1 at θ = 0.55, 0.6, 0.65, 0.7, 0.75, 0.8 (0.18571, 0.22857, 0.27143, 0.31428, 0.35714,
  0.40000; binding block P = Q = ν, N = N0). Partial keeping (e1_delta.py/.log, type II charged X0^δ): θ = 0.6 gives 0.21428 (δ = 1/16), 0.20571 (0.1),
  0.2 = E1 (1/8, 1/4, 3/8): ν = 0.22857 − 0.2286δ, so a loss below X0^{1/8} suffices there; δ = 1/2 is e1rate's (4/7)(2θ−1) (crossover 9/11);
  θ = 0.75, δ = 0.2: 0.33123 (E1 0.3125). e1rate's heuristic (DERIVATION §5 (a): "≤ (2/3)(2θ−1), beats E1 only for θ > 5/7") fixed one decomposition;
  it is not a ceiling of the method's own minimax (reviewer).
- Falsification of route (a) on finer grids (e1_fine.py/.log; G = 24, 13 values of n, 31 cuts, K = 24): θ = 0.6 worst block exponent −0.001875 at
  ν = 0.2275 and +0.001625 at ν = 0.2295 (zero at 0.22857 by interpolation; slope 7/4; binding corner P = Q = ν, N = N0): the value is grid-stable.
- LB with full keeping (e1_lb.py/.log; κ = kappa_opt.opt deg 6 in both columns, LB = kappa_opt.lb): 0.02231 → 0.03267 (θ = 0.51), 0.07827 → 0.09014
  (0.55), 0.11662 → 0.12984 (0.58), 0.14079 → 0.15496 (0.60), 0.15861 → 0.16716 (0.61); no gain for θ ≥ 0.6176 (κ((2/7)(3θ−1)) = c(θ)), none above 5/7.
- What a certificate needs: (i) a Deshouillers–Iwaniec Lemma 1 / Kuznetsov bound for Σ β_p γ_q α_n e(nq̄/p)Φ(n/(pq)) with the Bessel weight Φ (phase
  X0 ≍ T/H) inside, uniform over the blocks and coefficient classes, losing at most X0^δ with δ < 1/8 at θ = 0.6 — a hybrid Kloosterman–archimedean large
  sieve; none is known (e1rate's own escape route); (ii) the exponent minimax certified on the continuum (here sampled: G = 16 grid, 21 cuts, 17 bisections;
  reviewer: a grid-model prediction, the closed form (2/7)(3θ−1) is read off six θ); (iii) SL's remaining steps at ν < 1/2 (holds); (iv) κ(ν), LB in Arb.
- Route (b) dropped: setting type I to −∞ also voids piece 2 of the decomposition (a different model, not "type I averaged").
- Numeric slip in e1rate: "c(θ) ≥ c(9/11) = 0.5103": c(9/11) = 0.508194 (0.5103 is c(0.82)); its margin 0.0648 at 9/11 is computed correctly.
### 2.4 Verdict e1rate: OPEN (bookkeeping, conditional). The barrier is exact for its separated class (type II charged X0^{1/2}); the freedom "keep
  the Bessel phase on type II" passes it in the method's own bookkeeping for every θ ∈ (1/2, 0.8] (8/7·E1), tolerates X0^{1/8} at θ = 0.6, and would
  raise h(θ)/LB for θ < 0.6176. For θ ≥ 0.71399 nothing can change h with SL's one-piece detector (c(θ) ≥ κ(1/2)). Lesson 329's price: exactly the
  √X0 the method pays; the lead is conditional on an input that does not exist. Novelty: arXiv math.NT 2026-10/11 titles 0 hits; Zenodo
  '"short intervals" mollifier Kloosterman' nothing relevant; our files: only e1rate's escape-route sentence (DERIVATION l. 153).

## 3. gl2ceil (N8, N12) — one-piece Levinson–Conrey cap for GL(2)
### 3.1 Step fixed and class
- sup of c = (1 + w(1)²)/2 + (C_P/ξ)∫w² + ξB_P∫w′² (w = e^{Ry}Q, B_P = ∫P², C_P = ∫P′², ξ = ν/2 ≤ 25/256) over R > 0, real P with P(0) = 0,
  P(1) = 1, Q in Bernard's class; Arb dual certificate over (R, λ) ∈ [1/2, 28] × [0, 3] plus hand tails (R ≤ 0.5367, R ≥ 28, λ ≥ 3); inner P-minimum
  exact (sinh family); any degree (H¹ closure); complex Q with real P proved useless.
### 3.2 Freedoms checked
- Optimiser, degree, R: covered by construction (a dual bound over every cell, not a numerical optimum). Constraints P(0) = 0, P(1) = 1: the method's
  own. Own float check (gl2_xcheck.py/.log, e1rate's independent optimiser, Q(1) = 0 class): ξ = 25/256: 0.06294 (deg 6), 0.06922 (deg 10)
  ≤ 0.07005633; ξ = 1/8: 0.08691, 0.08954 ≤ 0.08971026 (no float above the caps; slow convergence at R ≈ 9.9, not a tightness test).
- Feng's first piece is INSIDE the certified class (new observation, exact): Σ(μ_f∗Λ_f)(n)n^{−s} = −L′/L² = Σ(−μ_f(n) log n)n^{−s}, so
  μ_f∗Λ_f = −μ_f·log (likewise μ∗Λ = −μ·log for ζ). A second piece Σ_{h≤y₂}(μ_f∗Λ_f)(h)h^{σ₀−1/2−s}P₂(log(y₂/h)/log y₂)/log y₂ is the one-piece
  mollifier with P ↦ P + G, G(x) = −((1−x)/τ)P₂((x−1+τ)/τ) on [1−τ, 1], 0 below (y₂ = y^τ, x = log(y/h)/log y): at τ = 1 exactly P − (1−x)P₂;
  G is continuous (P₂(0) = 0), G(1) = 0, G(0) = 0, so P(0) = 0, P(1) = 1 are kept and P + G lies in the H¹ closure. So the cap already covers it: no gain. Genuinely new are Feng's k ≥ 2 pieces (μ_f∗Λ_f∗Λ_f ≠ μ_f·log², the
  difference carrying Σ_{p|n} λ_f(p)² log²p-type weights). GL2Two (l. 115) and CeilAudit (l. 95) record no GL(2) mean value for them; GL2Two also
  reports AT21's off-diagonal bound as coefficient-free per pair, which would admit them at ν < 25/128 — unverified here.
- Two pieces of KRZ/BCY type: ≤ 2.1·10⁻¹¹ (gl2two, float). Complex (P, Q): AT21 l.996's √(AC − Im B²) ≤ the real value; GL2Ceil notes Bernard's
  mean value does not cover complex P; Degree2Mol cplx.py: random starts return to the real optimum (float).
- Skipped: the float test of Feng's k = 2 piece. PRZZ (arXiv 1802.10521, Thm meanvalueintegral, general-d contour combinatorics) gives the twisted
  moment for these coefficients but not a closed main-term functional I could implement and validate in the budget (download deleted).
### 3.3 Verdict gl2ceil: ROBUST in its class (certified), class now explicitly including Feng's k = 1 piece at any length ≤ y; open only for
  k ≥ 2 Feng pieces (needs their GL(2) main term) and complex (P, Q) (needs a mean value for complex P).

## Suggested NOGO lines
- explore63 (append, ceilaudit10): class extended — R(λ) < 2ν(λ) in floats on [1.10, 1.30] ∪ [25, 10⁵] (λ(2ν − R) → 1/3), so undamped and Jensen-damped
  MP_λ pre-images fail Newton k = 1 for every λ; any damping μ of a real-rooted pre-image in a basis with log-convex expansion coefficients obeys
  μ_m/μ_{m−1} ≤ (μ₁/μ₀)(b₁/b₀)(b_{m−1}/b_m)(N−m+1)/(Nm) (Lemma D, sharp; MP_λ every positive Φ: μ₅ ≤ 0.0522 at λ = 3/4); continuous Hahn has a Pólya
  lemma in floats (κ_k ≤ k!^{1/2}) but log-convex d_{2n} (μ₅ ≤ 0.0166); Fejér/Jackson/Jensen/Gaussian-damped MP approximants non-real from N ≤ 20;
  every real-rootedness preserver ε-faithful on one exponential e^{ax} over a window with aX = A has degree ≥ (A − 4ε)²/(16ε) (Lemma M′; Jensen A²/(2ε));
  not covered: preservers not close to the identity on exponentials [CeilAudit10](Riemann/Deep/CeilAudit10/REPORT.md)
- e1rate (rescope, ceilaudit10): the barrier is the separated class (DI type II charged √X0); in e1rate's own exponents.py a type II bound that keeps
  the Bessel phase at loss X0^δ gives ν < (2/7)(3θ−1) at δ = 0 (θ ≤ 0.8, grid) and still beats E1 at θ = 0.6 for δ < 1/8, raising LB for θ < 0.6176
  (0.1408 → 0.1550 at θ = 0.6); no such bound is known; LB > c(θ) ⟺ κ(ν) > c(θ), so nothing changes for θ ≥ 0.71399 (c(θ) ≥ κ(1/2)); c(9/11) = 0.5082,
  not 0.5103 [CeilAudit10](Riemann/Deep/CeilAudit10/REPORT.md)
- gl2ceil (append, ceilaudit10): class = real P (any degree, H¹ closure), every R > 0, Q in Bernard's class, and Feng's k = 1 piece (μ_f∗Λ_f = −μ_f log,
  absorbed by P) at any length ≤ y; outside: Feng k ≥ 2 pieces (no GL(2) main term), complex (P, Q) (float only), KRZ two-piece (≤ 2.1·10⁻¹¹, gl2two)
  [CeilAudit10](Riemann/Deep/CeilAudit10/REPORT.md)

RH relevance: N6 (Jensen/Pólya hyperbolicity tests: the T² law now holds for every pre-image route with log-convex coefficients and every ε-faithful
preserver), N8 (short-window and GL(2) critical-line proportions), N12 (GL(2) zero proportions). No record moved; no discovery claimed.

## Review and novelty
- Reviewer 1 (inline formulas): Lemma D symbolic CORRECT (sharp; one rounding fixed: μ₅ ≤ 0.0522), Lemma M′ CORRECT, the h(θ; c) = c identity CORRECT.
- Reviewer 2 (e1rate route (a)): monkeypatch CORRECT (only eII changed; Bdec reads module globals); the conditional reading is legitimate; heuristic (a)
  is invalid as a ceiling; GAP: the closed form (2/7)(3θ−1) and δ < 1/8 are grid-model predictions, LB values float only — stated so above.
- Novelty (OPEN item): python3 .scratch/lst.py on math.NT 2026-10/11 titles (short window|mollif × Kloosterman|Bessel|zeros, Levinson short): 0 hits;
  Zenodo '"short intervals" mollifier Kloosterman': nothing relevant (478 scanned); rg over our .md: only e1rate DERIVATION l. 153 (escape route).
- Lemma M′ and Lemma D are short Cauchy–Schwarz/Newton arguments (likely folklore); the continuous-Hahn Pólya lemma is a float observation only
  (300 random tests per normalisation), not checked against the literature.

## Files (Riemann/Deep/CeilAudit10/code/)
- lam_gap.py/.log (T1), damped_mp.py/.log (T2a), hahn_lemma.py/.log, hahn_lemma2.py/.log (T2b), hahn_pre.py/.log (T2c), budget.py/.log (Lemma D caps);
  e1_book.py → e1_book.log (κ table, ν*), e1_book2.log (superseded heuristic routes), e1_book_exact.py/.log; e1_routes.py → e1_routes_a.log (route a),
  e1_delta.py/.log (X0^δ scan), e1_lb.py/.log (LB with full keeping); gl2_xcheck.py/.log; STATE.md. Peaks ≤ 59 MB (memcap lines); compute ≈ 38 min in total, one process at a time.
- Read-only imports: Riemann/Deep/ExplorePositivity/mp_core.py and logs/c2n_zeta.txt, d2n_zeta.txt; Riemann/Deep/E1Rate/kappa_opt.py (__pycache__ removed).

