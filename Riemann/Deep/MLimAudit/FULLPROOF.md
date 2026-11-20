# mlimaudit FULLPROOF: the chain behind lim sup M(x)/√x > 2.0739001417, lim inf < −2.0764178466

Sources (local copies in .scratch/mlimaudit/src, deleted after the audit): OtR = Odlyzko–te Riele 1985 scan
(Odlyzko's site, mertens.disproof.pdf; PDF page = journal page − 136); BT = Best–Trudgian arXiv 1209.3843 TeX
(LinearIndependence-v14.tex); Hu = Hurst arXiv 1610.08551 TeX (ComputingMertensFunction.tex); W = Weisz
Zenodo 23651459 (THEOREM.md, data/README.md); C = Riemann/Deep/MLimsup/certify.py; L = Riemann/Deep/MLimAudit/lean_certify.py.

## S1. Theorem (OtR p. 144, read from the scan)
K ∈ C²(ℝ), K ≥ 0, K(−y) = K(y), K(y) = O((1+y²)⁻¹); k(t) = ∫K(y)e^{−ity}dy (OtR (2.10), p. 143) with
k(t) = 0 for |t| ≥ T and k(0) = 1. If the zeros ρ = β+iγ with |γ| < T have β = 1/2 and are simple, then for
every y₀: lim sup m(y) ≥ h_K(y₀) ≥ lim inf m(y), h_K(y) = Σ_ρ k(γ)e^{iγy}/(ρζ′(ρ)); m(y) = M(e^y)e^{−y/2}
(OtR's §2 notation; p. 142 not re-read here; Hu TeX l. 352–354 writes the same sum for q(x) = M(x)/√x with y = log x).
- The zero hypothesis is only on |γ| < T. For T = 18292.5 it is a published fact (van de Lune–te Riele–Winter
  1986: first 1.5·10⁹ zeros simple, on the line), independent of Weisz.
- OtR p. 141: RH and simplicity of all zeros follow from each one-sided bound lim sup M(x)x^{−1/2} ≤ A or
  lim inf ≥ −A (refs [18] Ingham, [44] Titchmarsh). So if RH fails or a zero is multiple, lim sup = +∞ AND
  lim inf = −∞. BT TeX lines 103–107 state the same (Ingham 1942 pp. 318–319). The page's sentence holds as written.
- Hu TeX lines 359–362 state the theorem only for the Jurkat–Peyerimhoff kernel; BT line 243 states the
  admissibility rule informally (nonnegative Fourier transform, support [−T, T]). The general-kernel form
  used by mlimsup is OtR's p. 144 statement above: VERIFIED against the scan.
- Real-pair form: k real and even, ρ̄ is a zero with ζ′(ρ̄) = conj ζ′(ρ), so h_K(y) = 2Σ_{0<γ<T} k(γ)Re(c e^{iγy}),
  c = 1/(ρζ′(ρ)). W THEOREM.md lines 32–39 use the same Re(c e^{iγy}); C and L compute cre·cos θ − cim·sin θ.

## S2. Kernel admissibility (C l. 63–76 load_kernel, l. 89–98 kernel_at; kernel_*.txt headers)
- g = Σ_j b_j 1[−T₁/2+jh, −T₁/2+(j+1)h), M = 1350, h = 1355/100, T₁ = Mh = 18292.5 (exact in binary: 36585/2).
  b_j integers, 40105 ≤ b_j ≤ 10⁶ (plus), 36239 ≤ b_j ≤ 10⁶ (minus): all > 0 (checked by L).
- (g∗g̃)(t) = ∫g(s)g(s−t)ds = Σ_{j,l} b_j b_l Λ(t−(j−l)h), Λ(s) = (h−|s|)₊. For t = (d+τ)h, 0 ≤ τ ≤ 1, only
  j−l ∈ {d, d+1} contribute: h[R(d)(1−τ) + R(d+1)τ], R(d) = Σ_j b_j b_{j+d}. ‖g‖² = hR(0). So
  k(t) = [R(d)(1−τ)+R(d+1)τ]/R(0) with R(M) = R(M+1) = 0: exactly C's kernel_at (R list padded with [0, 0]).
- k(0) = 1; k(t) = 0 for |t| ≥ Mh = T₁ (supp g has length T₁). k ≥ 0 everywhere since g ≥ 0 (not needed).
- K(y) = (1/2π)∫k(t)e^{ity}dt = |ĝ(y)|²/(2π‖g‖²) ≥ 0, even (g real), entire (k compactly supported) hence C²,
  ĝ(y) = Σ_j b_j(e^{iy(a_j+h)} − e^{iy a_j})/(iy) = O(1/|y|) so K = O(y⁻²) = O((1+y²)⁻¹). All S1 hypotheses hold.
- Independent recheck (L): R(d) from the fmpz_poly product b(x)·x^{M−1}b(1/x) equals C's double sum (asserted at
  d = 0, 1, 2, 500, M−2, M−1 and used for every d); mp_falsify (d): k(γ) by a breakpoint sweep of ∫g(x)g(x−t)dx.
- Rayleigh/eigenvector step (kernel.py) only CHOOSES b; it is not load-bearing: any integer b gives an admissible k,
  and the certificate evaluates the chosen k directly. Algebra checked: Σ_d q(d)R(d) = bᵀQb with Q_jj = q(0),
  Q_jl = q(|j−l|)/2 (j ≠ l), since R(d) for d > 0 counts each unordered pair once.
- The page's "optimum over all admissible kernels, up to the cell width" (Krein / Boas–Kac: a continuous
  positive-definite function on ℝ supported in [−T, T] is g∗g̃ with supp g ⊂ [−T/2, T/2]) is a side remark,
  not used by the bound; not audited further.

## S3. Zeros above γ₁₄₄₀₀ (C l. 138–164)
- Under S1's hypotheses, h_K(y) = 2Σ_{n≤N(T₁)} k(γ_n)Re(c_n e^{iγ_n y}). Split n ≤ 14400 (exact phases) and
  14401 ≤ n ≤ 20311: |2k Re(c e^{iθ})| ≤ 2|k||c| for every θ, so h_K(y₊) ≥ S₊ − P₊ and h_K(y₋) ≤ S₋ + P₋,
  P = 2Σ_{n=14401}^{N(T₁)} |k(γ_n)| |c_n|. Valid for any phases: no information about γ_n y is used.
- Completeness: N(T₁) = 20311 (Arb zeta_nzeros, exact). C asserts last_zero_used = N(T₁) for each kernel; L also
  asserts zero 20312 > T₁ and every zero ≤ 20311 < T₁. Zero 14401 > Weisz's γ₁₄₄₀₀ ball (C and L).
- Zeros 14401–20311 come from acb.zeta_zeros (rigorous isolation on the line, index by Turing's method; FLINT docs,
  acb_dirichlet_turing_method_bound: B good Gram blocks satisfying Rosser's rule give N(g_n) ≤ n+1, N(g_p) ≥ p+1); |c| from
  ζ′ as the order-1 Taylor coefficient on the zero's ball, 1/(|ρ||ζ′|) as an arb ball: rigorous upper bound.

## S4. Code audit of certify.py (C) and the input balls
- Ball rule: C `dec_ulp`/`ball()` (l. 78–83) = arb(mid) ± r ± u, u = one unit in the last printed digit: exactly W's
  data/README.md rule. arb(str) parses rigorously at the working precision. Arrays are strings: gamma <U9860
  (all 9860 chars), gamma_rad <U13 (W: all ≤ 8.6e-9840), c_re/c_im <U66 (60 sig. digits), radii <U12.
  No float underflow of radii (they are strings).
- Phases (l. 101, 110–120): p_hi = max(⌊9627·3.322⌋+400, ⌊9860·3.322⌋+64) = 32818 bits (page says 32,819: off by one,
  immaterial). |γy| ≤ 1.4·10⁴·5·10⁹⁶²⁶ ≈ 7·10⁹⁶³⁰ ≈ 2^31994, so ~824 bits (~248 digits) remain after the point;
  the data term |y|·rad(γ) ≤ 5·10⁹⁶²⁶·8.6·10⁻⁹⁸⁴⁰ ≈ 4·10⁻²¹³ dominates. Enough, with ~210 digits to spare.
  kk = ⌊mid(θ)/mid(2π)⌋ is any exact integer; θ − kk·2π is a rigorous ball for every choice (cos, sin 2π-periodic).
- Rounding to 320 bits by unary + (arb_set_round adds the error to the radius); sums at 320 bits: rigorous.
- kernel_at (l. 89–98): exact integer R, asserts the γ/h ball lies in one cell [d, d+1]; 0 beyond M. Correct.
- JP check (l. 124): f(u) = (1−u)cos πu + sin πu/π, u = γ/γ₁₄₄₀₀: W's kernel.
- Extension (l. 138–164, 96 bits): T₁ exact; N(T₁) exact (assert); zero 14401 on the line and above γ₁₄₄₀₀;
  penalty 2|k||c| with |c| an arb ball from ζ′ on the zero ball. A ball straddling T₁ would trip an assert.
- Final test (l. 166, 175–179): arb `>`/`<` are certain comparisons (whole ball), hk = tot − pen (resp. + pen).
  It certifies exactly h_K₊(y₊) > 2.0739001417 and h_K₋(y₋) < −2.0764178466, i.e. the page's two claims via S1.
- Cosmetic: certify.log's "± 4·10⁻⁴⁰" and "± 3.2·10⁻²⁴" are print rounding of str(40)/str(30); true radii (L):
  known sums 9.9·10⁻⁶⁰ (as W's 1.3·10⁻⁵⁹), final bounds 2.5·10⁻³⁴ (L at 128 bits). Conservative either way.
- Not load-bearing: SHA-256 vs the package's own SHA256SUMS.txt (integrity only); --phases floats (kernel.py input).
- Memory: C holds gamma as numpy <U9860: 14400·9860·4 B = 568 MB, plus zip bytes: the 0.9 GB. L streams rows.

## S5. Falsification and independent recomputation (logs in this directory)
- lean_certify.py (lean_certify.log): S₊, S₋, P₊, P₋ and both bounds reproduced to all printed digits of certify.log;
  Weisz's JP certificates reproduced inside his radii (|Δmid| 3.6·10⁻⁶⁰, 4.6·10⁻⁵⁹). Peak 28 MB (memcap).
- mp_falsify.py (mpmath only; mp_falsify.json):
  (a) 50 extra zeros n = 14401 + 120k (k < 49) and 20311: |γ_mp − γ_Arb| ≤ 4.7·10⁻²¹; every Arb |c| upper bound
      ≥ mpmath 1/|ρζ′(ρ)|, at most 5.7·10⁻¹⁵ relative above (print rounding of the 15-digit file).
  (b) γ₂₀₃₁₁ = 18291.6524 < 18292.5 < γ₂₀₃₁₂ = 18292.8865 (mpmath zetazero): N(T₁) = 20311 confirmed.
  (c) θ/π mod 2 at 33,300 bits for zeros 1, 2, 3, 72, 1000, 5000, 10000, 14399, 14400 at y±: agree with Arb to 4.9·10⁻⁵⁰
      (the dump's 50 digits).
  (d) k(γ) by a direct breakpoint sweep of ∫g(x)g(x−γ)dx (no R, no interpolation formula) at the same zeros, both
      kernels: agree to 3.7·10⁻³¹ (dump's 30 digits); the terms 2k Re(c e^{iθ}) agree to 3.4·10⁻³¹.
  (e) c = 1/(ρζ′(ρ)) from Weisz's γ for zeros 1, 2, 72, 14400 at 75 digits vs his c strings: ≤ 4.8·10⁻⁶².
- mp_fullsum.py (mpmath, all 14,400 zeros, phases at 33,300 bits, 80-digit sums; 25 MB): S₊, S₋, JP₊, JP₋ agree with
  the Arb sums to ≤ 4.4·10⁻⁴⁵ (print precision of the Arb side). Library-independent check of the whole known part.
- fullprec_newton.py (own code, Arb at 33,100 bits; ~65 s and ~190 MB RSS per zero): fresh random zeros (seed 20261120,
  range 73–14400): zeros 7338, 8578, 9656: N(B) ⊂ B, Im N(B) in his γ ball, c(N(B)) in his c balls: 3/3
  (|ζ(s₀)| ≈ 10⁻⁹⁸⁴⁸ at his midpoints; memcap peak 189 MB). B = 1/2 + iγ^W ± 2(rad + ulp) in both coordinates.

## S6. The condition, narrowed (lowprec_all.py, lowprec_all.json; peak 24 MB, 326 s)
- All 14,400 Weisz zeros against Arb's acb.zeta_zeros(n, ·) at 96 bits: |γ_n^W − γ_n^Arb| ≤ 5.2·10⁻²⁵ for every n
  (index n and ~28 digits confirmed; 0 mismatches); δ_n := upper|c_n^W(mid) − c_n^Arb| ≤ 6.8·10⁻²⁵, and every Weisz c ball
  meets our rigorous c ball (0 disagreements).
- c part removed: h_K(y) computed with Weisz's c midpoints differs from the value with the true c_n by at most
  E = 2Σ_n |k(γ_n)| δ_n = 7.9·10⁻²³ (plus), 7.7·10⁻²³ (minus) ≪ margins 5.98·10⁻¹¹ (lim sup), 1.09·10⁻¹¹ (lim inf).
- γ part, tolerance: if every true ordinate satisfies |γ_n − γ_n^W(mid)| ≤ ε, each phase moves by ≤ |y|ε, so h moves by
  ≤ |y|ε·W + (k-variation: |k′| ≤ 1/h, negligible), W = 2Σ|k||c| = 2.7408 (plus), 2.7363 (minus). The bounds survive for
  ε < margin/(W|y|): ε₊ = 1.7·10⁻⁹⁶³⁶ (|y₊| = 1.27·10⁹⁶²⁵), ε₋ = 8.0·10⁻⁹⁶³⁹ (|y₋| = 4.96·10⁹⁶²⁶).
- So the record needs only: the published ordinates γ₁ … γ₁₄₄₀₀ are correct to 8·10⁻⁹⁶³⁹ (≈ 9,643 significant digits);
  Weisz states radii ≤ 8.6·10⁻⁹⁸⁴⁰. Spot checks of this at full precision: 50 zeros (spotcheck2.log, earlier run) + 2
  (spotcheck.log); 4 more here only at 75 digits through c. Not dischargeable cheaply: ~25–50 s per zero at 33k bits.
