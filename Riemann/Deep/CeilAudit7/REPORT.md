ROBUST: explore58 — adding Gonek-type off-diagonal shifted means D(a,b) = sinc(a−b) − sinc(a)sinc(b) (sign-carrying; exact CUE limit) and CGG derivative rows leaves the δ = 9/10 hard-core cone feasible (float LP), and the sign criteria fail only at shifts ≥ 1; explore60 — Gaussian block mixtures matching 𝒮 and ζ's one-point fourth-order (⟨Z⁴⟩/⟨Z²⟩² = 30, ⟨Z²Z′²⟩/⟨Z⁴⟩ = π²/15) and L¹ moments have 10.1% zero density, so these inputs cannot pass 57.735% (shift-dependent fourth-order means not covered); explore44 — 7/11 = 2/(2 + κ) is fixed by Theorem U′'s step 6e′ (hairpin majorant ≍ X^{2d_ρ} against density slope κ), and Lemma M, Hölder exponents and the AM–GM step do not move it (robust for that bounding scheme; reviewer GAP on universality); explore54 and explore36 skipped (reasoned notes: ROBUST as stated).

# ceilaudit7 — falsification of our N8/N9/N7/N1 barriers explore58, explore60, explore54, explore44, explore36 (Orca, 20 Nov 2026)

Brief: .scratch/orca/ceilaudit7-brief.md. Scratch: .scratch/ceilaudit7/. Start 11:50 UTC.

## 0. Choice
- Chosen (best chance first): A = explore58 (sign-carrying off-diagonal shifted means are a cheap exact test), B = explore60 (4th-order/L¹ inputs vs the Gaussian), D = explore44 (which step fixes 7/11).
- Skipped tests (reasoned verdict only, §C, §E): explore54 (its structure lemma is algebra for every q, a_j, λ, X), explore36 (stated numerical; positivity reaches less, 13.9c ≪ 2πc²).

## A. explore58 (N8 small gaps): hard-core CUE mixtures vs untwisted inputs
- Step fixed: the cone test E(δ; I): ζ's input vector v_Z(I) lies in the closed convex cone of v(f), f real-rooted type-π (trig polys), all gaps ≥ δ.
  Class: I = (m_0..m_10, d = Gonek's Σ|ζ′(ρ)|², D(α) = Gonek's Σ|ζ(ρ+iα)|² at α ∈ {¼, ½, ¾, 1, 3/2, 2}); criteria linear or concave 1-homogeneous in v.
- Every listed input is a square: m_j = ⟨Z^{(j)2}⟩, d = ⟨Σ_γ Z′(γ)²⟩, D(α) = ⟨Σ_γ Z(γ+α)²⟩. None carries a sign.
- Freedom A1 (sign information Gonek's lemma keeps): the off-diagonal shifted means D(a, b) := ⟨Σ_γ Z(γ+a)Z(γ+b)⟩/⟨Z²⟩ (normalised units),
  i.e. Σ_γ ζ(ρ+iα)ζ(1−ρ−iβ) times the phase e^{i(α−β)L/2}; same Gonek-lemma computation as the diagonal (under RH).
- Sign criterion (exact, elementary): if every gap is ≥ δ (so all zeros simple), then for 0 < a, b < δ the interval (γ − b, γ + a) holds exactly one
  zero (γ), so Z(γ+a)Z(γ−b) < 0 for every γ: D(a, −b) ≤ 0 on the whole hard-core class (any real function, type π not needed).
  Linear in the extended input vector. Hence D_ζ(a, −b) > 0 for some a, b < δ ⟹ (under RH) infinitely many gaps < max(a, b), i.e. μ ≤ max(a, b).
  Also: Z′(γ)Z(γ+a) > 0 for 0 < a < δ (no sign change before the next zero): E(a) := ⟨Σ_γ Z′(γ)Z(γ+a)⟩ ≥ 0 on the class.
- Exact CUE value (Palm + Heine): S(u, v) := E Σ_j Z_U(θ_j+u)Z_U(θ_j+v) = h(0)·D_{N−1}(w h), w = |1 − e^{iθ}|², h(θ) = 4 sin((u−θ)/2) sin((v−θ)/2),
  D_n = Toeplitz determinant (pentadiagonal symbol); D(a, b) = S/(N(N+1)), u = 2πa/N, v = 2πb/N. Checks: a = b → 1 − sinc²(a); D(a,a)/a² → π²/3 = d.
- Result A1 (dab.py, logs/anti_N800.log): the exact CUE values converge (O(1/N)) to
  **D(a, b) = sinc(a − b) − sinc(a)·sinc(b)** (sinc x = sin πx/πx), the covariance of the flat-spectrum Gaussian conditioned on Y(0) = 0.
  N = 800: D(a,a) = 1.001248 at a = 1 (limit 1); D(a,−a) = −0.406028 at a = ½ (sinc 1 − sinc²½ = −0.405285); +0.101926 at 1.2 (+0.10186).
  Consequences (exact in the limit): D(a, −b) = sinc(a + b) − sinc(a)sinc(b) ≤ 0 on (0, 1]², = 0 only at (1, 1); D(a, −a) > 0 iff a ∈ (1, 1.42);
  E(a) = −sinc′(a) > 0 on (0, 1.4303). So the sign criteria are violated by ζ's inputs only at shifts ≥ 1 = the mean spacing: they give μ ≤ 1, no gain.
- Freedom A2 (CGG-type higher derivatives at zeros, Taylor coefficients of D): Σ Z″(γ)² = 4π⁴/45, Σ Z′(γ)Z‴(γ) = −π⁴/5 (normalised, per zero; d = π²/3).
  Exact local identity: f(γ+t) = f′(γ)t·exp(At − Bt²/2 + …), B = Σ_k (x_k − γ)^{−2}, so (3/4)f″(γ)² − f′(γ)f‴(γ) = 3B f′(γ)², and gaps ≥ δ give B ≤ π²/(3δ²).
  Linear criterion ⟨Σ_γ ((3/4)Z″² − Z′Z‴)⟩ ≤ (π²/δ²)·d on the class; ζ: 4π⁴/15 = (4π²/5)·d, so it fails only if δ > √(5/4) = 1.118: no gain.
  (Energy weighting again: the Z′²-weighted mean of B is 4π²/15 = 2.63 < the lattice's π²/3 = 3.29; unweighted sine-kernel mean is 2π²/3.)
- Test A3 (cueAB.py = explore58's lk_cueD columns + off-diagonal rows; logs/A2_*.log), N = 48, δ = 9/10, rows m_0..m_6, d, D(α) at explore58's six α, plus
  - 6 pairs (¼,−¼), (½,−½), (¾,−¾), (½,−¼), (0.8,−0.4), (¼,½), scales s ∈ [½, 1], S = 400: feasible, residual 0 (20 basic columns); same pool without pairs: 0 (14).
  - 10 pairs (adds (0.85,−0.85), (0.85,−0.6), (0.4,0.8), (0.1,0.85)), density-1 columns only (s = 1), S = 600: mean relative residual 1.2·10⁻⁷ (21 columns), i.e. feasible to LP tolerance.
  S = 60 is not enough (2.4·10⁻³): pool size, as in explore58.
- Non-concave ratio criteria: a criterion valid on each local window reaches ζ's averaged inputs only if it survives mixtures; Cauchy–Schwarz-type
  ratio tests (√(v_i v_k) − |v_j| ≥ 0) are concave 1-homogeneous, hence inside the cone statement. Criteria not closed under mixtures do not apply to window averages.

  Off-axis check (logs/A_dabchk.log): Richardson (N = 400, 800) equals the formula to ≤ 6·10⁻⁶ at (½,−¼), (0.8,−0.4), (¼,½), (1.3,0.2), (0.7,−1.1), (2,−½).
- ζ itself (zcheck.py, logs/A_zcheck.log; 250 zeros at t ≈ 1.5·10⁴, L = 7.8, Z by mpmath, ⟨Z²⟩ = L + 2γ_E): diagonal ratios ζ/limit 1.38 (a = ½), 1.25 (a = 1),
  off-diagonal D(¼,−¼) = −0.246 (limit −0.174), D(½,−½) = −0.492 (−0.405), D(0.9,−0.9) = +0.042 (−0.116), D(1.2,−1.2) = +0.079 (+0.102).
  Lower-order terms of relative size 0.2–0.4 at this height (explore58 saw 1.33 for d/m_1 below 237); they vanish as L → ∞, so they cannot enter a statement about μ.
  The identity D_ζ = D_CUE at leading order is Gonek's lemma under RH (diagonal: Gonek 1984); not re-derived here.
- Verdict A: ROBUST, class extended. Exact class now covered (float): energy-weighted untwisted second-order inputs = m_0..m_6, d, D(α, α), and off-diagonal
  D(a, b) on finite grids (sign-carrying), hence every finite set of mixed means Σ_γ Z^{(j)}(γ)Z^{(k)}(γ + α) up to finite differences; linear or
  concave 1-homogeneous criteria; δ = 9/10. Reason: ζ's energy-weighted Palm 2-point function is the Gaussian bridge sinc(a−b) − sinc(a)sinc(b); its sign
  structure is that of a function with no zero within distance 1 of each zero, and its derivative data weight isolated zeros (⟨B⟩_w = 4π²/15 < π²/3).
  Not covered: energy-weighted inputs of order ≥ 4 at zeros; Arb certificate for the new rows (float LP only); unweighted zero statistics (pair correlation) are outside the class by design.

## B. explore60 (N8 critical line): the Gaussian ceiling 1/√3 against fourth-order and L¹ inputs
- Step fixed: (C3) the Gaussian X(t) = √(2θ′)Y(θ(t)/π), Y flat on [−π, π], lies in 𝒮 and has Rice density (1/π)√(λ₂/λ₀) = 1/√3. Class 𝒮: all
  untwisted second-order unfolded means ∫F^{(i)}(v)F^{(j)}(v + u)dv → c(−1)^i sinc^{(i+j)}(u). Criteria: any theorem "every f ∈ 𝒮 has ≥ c·N zeros".
- What the Gaussian fails to match (untwisted, outside 𝒮), leading order:
  - ⟨Z⁴⟩/⟨Z²⟩² = 3 for the Gaussian, ~ L²/(2π²) for ζ (Ingham; CUE (N+2)²(N+3)/(12(N+1))) → ∞;
  - R := ⟨Z²Z′²⟩/⟨Z⁴⟩ (unfolded units): Gaussian λ₂/3 = π²/9 = 1.0966; ζ/CUE (b4.py ratio, Toeplitz D_N((2 − 2cos θ)(2 − 2cos(θ − φ))), logs/B_ratio.log):
    0.68429, 0.67113, 0.66446 at N = 100, 200, 400 → **π²/15 = 0.65797** (1/N extrapolation 0.6578): ζ's large values are smoother than Gaussian ones;
  - ⟨|Z|⟩/⟨Z²⟩^{1/2} = √(2/π) = 0.798 for the Gaussian, ≍ L^{−1/4} → 0 for ζ (Heath-Brown, Heap–Radziwiłł–Soundararajan); ≈ 0.45 at T = 3·10¹² (Keating–Snaith constant, rough).
  - Finite-T second-order data (Rice .682/.629/.601 at T = 1500/10⁶/3·10¹²) are inside 𝒮 at that T; an asymptotic proportion theorem sees only the limit 1/√3.
- Freedom tested (amplitude/length decoupling): every input above is homogeneous of degree ≥ 1 in Z, i.e. energy-weighted, while the zero count is length-weighted.
  Construction: block mixture (long blocks, switching) of stationary Gaussians G_k with spectra uniform on bands B_j ⊂ [0, π] (J bands), variances A_k, length fractions p_k,
  plus a dead remainder (amplitude → 0). Exact constraints, linear in p: Σ_{k∈B_j} p_kA_k = 1/J (energy-weighted spectrum exactly flat: second order = 𝒮 for all i, j, u);
  3Σp_kA_k² = M4; Σp_kA_k²ℓ_k = R·M4 (ℓ_k = λ₂ of band); √(2/π)Σp_k√A_k = M1; Σp_k ≤ 1. Zero density Σp_k√ℓ_k/π (Rice per block; block ends O(1/length)). LP over 61 amplitudes × J bands.
- Test (b4.py lp, logs/B_lp.log), R = π²/15 throughout:
  | M4 | M1 | J | min zero density (1/√3 = 0.5774) | live length |
  |---|---|---|---|---|
  | 30 (≈ ζ at T = 3·10¹²) | 0.45 | 8 | 0.1074 | 0.45 |
  | 30 | 0.45 | 16 | 0.1008 | 0.56 |
  | 10 | 0.60 | 8 | 0.2579 | 0.60 |
  | 100 | 0.30 | 8 | 0.0432 | 0.25 |
  | 3.2 | 0.78 | 8 | infeasible (R = π²/15 needs non-Gaussian M4) | – |
  Matching ζ's second-order class 𝒮 and its fourth-order and L¹ one-point moments leaves the zero density far below 1/√3 (10.1% at the T = 3·10¹² values), → 0 as M4 → ∞ (live length ≥ 3/M4).
- Verdict B: ROBUST. Exact class: 𝒮 plus finitely many one-point moments ⟨|Z|⟩, ⟨Z⁴⟩, ⟨Z²Z′²⟩ (more such rows add columns, not obstructions, as long as each is homogeneous); ceiling ≤ 1/√3 (in fact ≤ 0.11 at ζ's T = 3·10¹² moments).
  Not covered: shift-dependent fourth-order means and twisted moments. Measured (b4g.py, logs/B_G4.log): ζ/CUE's G(x) = ⟨Z(0)²Z(x)²⟩/⟨Z⁴⟩ = 0.722, 0.302, 0.0756, 0.0189, 0.0047
  at x = ½, 1, 2, 4, 8 (≈ 0.30/x²), while every Gaussian block mixture has G ≥ 1/3 (Wick): these means separate ζ from all such mixtures.
  Heuristic only: X = A·Y with A > 0 independent of Y keeps Y's zeros and r_X = R_A·r_Y; r_X = sinc forces λ₂(Y) = π²/3 − |R_A″(0)| < π²/3, i.e. density below 1/√3,
  if sinc/R_A is positive definite (true to first order for small amplitude variance; unchecked for ζ-sized fluctuations).


## D. explore44 (N9/N7 → N8 side): "Theorem U′'s 6e′ needs 2θ < κ by itself, so 7/11 stays"
- Where 7/11 comes from: Theorem U′ runs DPC's method with mollifier R_T of length X = T^θ, θ ∈ (1/φ − 1, min{κ/2, 1}); λ = 1/(1 + θ) > 1/(1 + κ/2) = 2/(2 + κ) = 7/11 at Conrey's κ = 8/7.
  Two steps impose θ < κ/2 (UPAudit/FULLPROOF.md): M3/M4c (moments of Off, pθ < κ with p = 2u > 2) and 6e′(d)–(e).
- Step fixed (6e′): hairpins from the cut points x_ρ = γ_ρ ∓ h/2 down to depth d_ρ = β − ½ ∈ [a/L, (a+1)/L) (set S_a, |S_a| ≤ 2σ_aTL, σ_a ≤ C₅e^{−κa} by (D)):
  𝓗_{T,a} ≤ 4C_P Σ_±Σ_{S_a}∫_0^{d_ρ}|R_TR*|(x_ρ − iy)·Z dy, with (Z+) |R_T(t − iy)| ≤ C(1 + W_{−y}), (Z−) |R*(t − iy)| ≤ CX^y(1 + W),
  W_{−y} ≤ On + X^{−y}Off, W ≤ On + Off, Off(t) = Σ_{d_ρ>0}X^{d_ρ}κ_A(L|γ − t|). Hence |R_TR*| ≤ C²(1 + Off)(1 + On)^{2k}[X^y + Off], and
  (b) (1 + Off)(E + Off) ≤ 3(E² + Off²), E = e^{θ(a+1)}; (d) Σ_{S_a}Off²Y via Hölder (v, u, w) needs Σ Off^{2u} ≪ TL (2uθ < κ); (e) Σ_a e^{2θ(a+1)}σ_a^{1−1/w} < ∞ needs 2θ < κ(1 − 1/w).
  Class: bounds of the hairpin terms through |R_TR*| (absolute values) and a zero-density input of slope κ uniform in a.
- Freedoms the cited lemmas admit, and why none moves 7/11:
  1. Lemma M for every q (explore44's point): irrelevant to (e); the e^{2θa} weight is there with q = 2.
  2. AM–GM in (b) is lossy only in form: for ρ ∈ S_a the zero's own term gives Off(x_ρ) ≥ X^{d_ρ}κ_A(Lh/2) ≍ e^{θa}, so the unexpanded cross term E·Off and the term Off²
     are both ≍ e^{2θa}: Σ_{S_a}(E·Off + Off²)Y ≳ e^{2θa}σ_aTL. No choice of Hölder exponents (v, u, w) removes it; they only cost factors (1 − 1/w).
  3. Is the weight X^{2d} an artefact of Lemma Z's upper bounds? On the hairpin the zero's own residue gives |M(s₊)| ≍ X^{d−y}, |M(s₋)| ≍ X^{d+y}: product X^{2d} for every y ∈ (0, d)
     (no phase rotation along the hairpin: the factors X^{±i(γ−t)} cancel at t = x_ρ fixed). Float test below: attained up to log-factors.
  4. A weighted density estimate in place of (D): Σ_ρ X^{2d_ρ} ≪ TL for θ ≥ κ/2? By Littlewood's lemma with a mollifier of length T^{θ′} (the source of κ = 2θ′, Conrey θ′ = 4/7),
     Σ_ρ[X^{2d} − 1 − 2θLd] = ∫(2θL)²X^{2(σ−½)}Σ_{β>σ}(β − σ)dσ ≲ T(2θL)²∫_0^{½}e^{2(θ−θ′)Lu}du: bounded by ≍ TL iff θ < θ′ = κ/2. Same threshold: 7/11 = 1/(1 + θ′), i.e. DPC's mollifier must be shorter than Conrey's.
  5. The hairpin shape is fixed by the cut of V from x_ρ to the zero; the y-direction (lower vs upper half-plane) is symmetric (s₊ ↔ s₋).
  So the floor moves only with κ (a longer density mollifier, N7/N9) or with cancellation between hairpins (signs of the jumps of P(V) across different zeros), which no absolute-value bound sees.
- Test (dh_hairpin.py, logs/D_hairpin.log; Davenport–Heilbronn, a periodic mod 5, b = Dirichlet inverse, M_X(s) = Σ_{n<X} b_n(1 − log n/log X)³n^{−s}, X = 10³…2·10⁵, 0.3 s, 19 MB):
  | DH zero | d = β − ½ | y | slope of log|M(s₊)M(s₋)| in log X (X ≥ 10⁴) | 2d |
  |---|---|---|---|---|
  | 0.8085 + 85.699i | 0.3085 | d/2, 0.9d | 0.384, 0.438 | 0.617 |
  | 0.6508 + 114.163i | 0.1508 | d/2, 0.9d | 0.264, 0.275 | 0.302 |
  | 0.5744 + 166.479i | 0.0744 | d/2, 0.9d | 0.228, 0.231 | 0.149 |
  | t = 94.5 (no zero), y = 0.1 | – | – | −0.075 | 0 |
  |M₊M₋| grows like a power X^{c(d)}, c ≈ 2d up to the kernel's log X factors (deepest zero: (log X)^{−O(1)} still visible at X ≤ 2·10⁵; shallow zero: other zeros' terms): the weight is real, not a bound artefact.
- Verdict D: ROBUST, and explore44's sentence is right but should name the step: 6e′(b)–(e), hairpin weight |R_TR*| ≍ X^{2d_ρ} (zero's own residue in (Z±)) against density slope κ:
  converges iff 2θ < κ, i.e. λ > 2/(2 + κ). Exact class: absolute-value hairpin bounds + any density input of slope κ (uniform in a); escapes: κ > 8/7, or hairpin cancellation.

- Review (reviewer subagent, verdict-only, claims (1)–(4) of this section pasted inline): **GAP** — (1) Off(x_ρ) ≳ e^{θa} (not ≍): the e^{2θa} cost survives in the
  majorant but does not lower-bound the true hairpin contribution; (2) one explicit-formula summand does not prove |M(s±)| ≍ X^{d∓y} (cancellation with other terms possible),
  DH numerics are evidence only; (3) CORRECT as a majorant calculation (O(TL) iff θ < κ/2, ≍ TL² at equality), but it does not constrain the true weighted sum;
  (4) the threshold limits this majorant/density/Hölder scheme, not every absolute-value method: sharper mollifier estimates or zero–mollifier correlations could improve it
  without a larger κ or hairpin cancellation. "The 7/11 threshold is justified for this bounding scheme; an unavoidable floor is not proved."
- Verdict D after review: ROBUST for the exact class "6e′ as written: Lemma Z's pointwise majorants (Z±) through On/Off, a density input of slope κ uniform in a, Hölder in y and over S_a";
  7/11 = 2/(2 + κ) is fixed there (θ < κ/2). Escapes outside the class: κ > 8/7; cancellation between hairpins; hairpin estimates sharper than Lemma Z near off-line zeros
  (zero–mollifier correlations; DH shows the majorant's growth X^{c(d)} with c(d) > 0 is real on that example, not that it is forced).

## C. explore54 (skipped: no test; reasoned note)
- Step fixed: the structure lemma: the SMT on R_X (every q, a_j, λ, X, every meromorphic g) = Littlewood's identities + a pointwise inequality in (w, w′) = (g, λg′);
  its only analytic content is the log-derivative lemma, an upper bound for ∫log⁺|g′/(g − a)|, i.e. a lower bound for a-point counts. Class: Nevanlinna/Ahlfors/Milloux-type inequalities for g, g^{(k)} on the half-strip.
- Freedoms in the brief: q ≥ 3 targets are already in explore54's grid (q = 3, 4, 6: −482 … −1715 at T = 10⁴); Milloux, T(r,f) ≤ N̄(r,∞) + N(r,0) + N̄(r,1;f^{(k)}) − N₀(r,0;f^{(k+1)}) + S(r,f),
  again puts the zero counts of f on the large side (a lower bound for N(r,0;f)); Ahlfors' covering inequality is void (L/S = 3.88 per period in the RH model).
- Reasoned verdict (untested): ROBUST. In this theory the only upper bound for an a-point count is the FMT, N(r,a) ≤ T(r,g) + O(1) = Littlewood + Jensen on a base line (Levinson's step).

## E. explore36 (skipped: no test; reasoned note)
- Step fixed: Le Cam's two-point lemma for arguments from D_X = {Λ(n) : n < X} with a positive margin η (relative Weil-window operator norm, windows c ≤ X).
  Class: alternatives = ζ's zeros with one pair merged off the line (δ ≤ 0.2) and ≤ 800 neighbours moved on the line by < 1 spacing; distances in float64 (Levenberg–Marquardt on the HS norm, Nyström operator norm).
- Certified? No: float numerics on a finite family; the 2πX² threshold is read off M-scans (ρ ∝ M^{−1.5} below ½log(t/2π), flat above); explore36 itself calls exact invisibility conjectural (D3).
- Freedoms: (i) Weil-window positivity (sign information on the zero side) is margin-free, hence outside Le Cam, but reaches less (13.9c, explore19) than 2πc²;
  (ii) higher moments or any nonlinear functional of the window data are functions of D_X: covered whenever they have a margin;
  (iii) outside the class: the sign and support of the prime side beyond X (Λ(n) ≥ 0 on log ℕ). Every finite compensation violates it (ΔG is a smooth exponential sum on |y| > log X,
  not a discrete measure), and a conductor-1 functional equation leaves only ζ (Hamburger). For near-line pairs, one-sided arguments using these tail signs must beat the dropped tail mass
  with the pair's gain cosh(δy) − 1 ≈ δ²y²/2, which is why such arguments give zero-free regions only within c/log t of Re s = 1 (classical).
- Reasoned verdict (untested): ROBUST as stated (numerical, positive-margin, data below X); the line should add "local on-line compensation, float64".

## Suggested NOGO lines
- ceilaudit7 → explore58 (N8 small gaps) ROBUST, class extended: Gonek-type off-diagonal means D(a,b) = ⟨Σ_γZ(γ+a)Z(γ+b)⟩/⟨Z²⟩ = sinc(a−b) − sinc(a)sinc(b) (exact CUE limit, Toeplitz) carry signs, but violate the hard-core sign pattern (D(a,−b) ≤ 0, Σ Z′(γ)Z(γ+a) ≥ 0 for a, b < δ) only at max(a, b) ≥ 1; ten such rows keep the δ = 9/10 cone feasible (float LP, density-1 columns); the CGG derivative rows give (3/4)Z″² − Z′Z‴ = 3BZ′² with ⟨B⟩ = 4π²/15, useless below δ = 1.118 [CeilAudit7](Riemann/Deep/CeilAudit7/REPORT.md)
- ceilaudit7 → explore60 (N8 critical line) ROBUST: ζ's untwisted inputs outside 𝒮 that the Gaussian misses (⟨Z⁴⟩/⟨Z²⟩² ~ L²/2π², ⟨Z²Z′²⟩/⟨Z⁴⟩ = π²/15 vs π²/9, ⟨|Z|⟩/⟨Z²⟩^{½} → 0) are energy-weighted; block mixtures of band-limited Gaussians matching 𝒮 and these one-point moments (M4 = 30, M1 = 0.45) have zero density 10.1% (LP), so they cannot lift 57.735%; shift-dependent fourth-order means not covered [CeilAudit7](Riemann/Deep/CeilAudit7/REPORT.md)
- ceilaudit7 → explore44 (N8 side, Theorem U′) rescoped: 7/11 = 2/(2 + κ) is fixed by step 6e′: hairpin majorant |R_TR*| ≲ (1 + Off)(1 + On)^{2k}(X^y + Off) with Off(x_ρ) ≳ X^{d_ρ} against density slope κ, Σ_a e^{2θa}e^{−κa} < ∞ iff θ < κ/2; Lemma M for every q, Hölder exponents and the AM–GM step do not move it; DH's mollifier grows like X^{0.38–0.44} on the hairpin of its d = 0.31 zero; escapes: κ > 8/7, hairpin cancellation, estimates sharper than Lemma Z near off-line zeros (reviewer: floor proved for this scheme only) [CeilAudit7](Riemann/Deep/CeilAudit7/REPORT.md)
- (reasoned, untested) explore54: scope sentence to add: "Nevanlinna/Ahlfors/Milloux-type inequalities for g, g^{(k)} on the half-strip bound a-point counts only from below; the only upper bound is the FMT = Littlewood + Jensen (Levinson's step)" [CeilAudit7](Riemann/Deep/CeilAudit7/REPORT.md)
- (reasoned, untested) explore36: scope sentence to add: "float64, local on-line compensation by ≤ 800 neighbours; positivity (13.9c) and functionals of the window data are inside or weaker; the prime-side sign/support beyond X (Λ ≥ 0 on log ℕ) is outside and kills every finite compensation" [CeilAudit7](Riemann/Deep/CeilAudit7/REPORT.md)

RH relevance: N8 (small gaps μ: explore58 extended; critical-line proportion: explore60 extended; Theorem U′'s floor 7/11: explore44 rescoped to step 6e′); N7/N9 (κ is the only lever on 7/11 inside the scheme); N1, N13 (reasoned notes on explore36, explore54).

## Files
- Riemann/Deep/CeilAudit7/code/: dab.py (exact CUE D(a,b) via Toeplitz), dabchk.py (off-axis check), zcheck.py (ζ zeros at 1.5·10⁴), cueAB.py (explore58's cone LP + off-diagonal rows; imports ExploreLK/lk_cue.py),
  b4.py (CUE ⟨Z²Z′²⟩/⟨Z⁴⟩; Gaussian block-mixture LP), b4g.py (G(x)), dh_hairpin.py (DH hairpin growth).
- Riemann/Deep/CeilAudit7/logs/: anti_N800.log, A_dabchk.log, A_zcheck.log, A2_*.log, B_ratio.log, B_lp.log, B_G4.log, D_hairpin.log.
- Compute: ≈ 22 min wall in total (15 of them a first ζ run killed by memcap at 900 s with no output: 1500 zeros plus 30000 Z² samples on a background core); every run under memcap 290 MB, peaks ≤ 48 MB (top), footprint ≤ 28 MB (/usr/bin/time -l).
- Novelty: no OPEN item, so no novelty scan (the D(a,b) identity is a Gonek-type mean value; not checked against the literature).

