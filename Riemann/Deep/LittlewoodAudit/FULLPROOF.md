# LittlewoodAudit: the chain of littlewood-half-plane, checked step by step

Source: `Riemann/Deep/Littlewood78/PROOF.md` (cited as P:<line>), `cert.py` (C:<line>), `cert_output.txt` (O:<line>).
Dusart = arXiv:1002.0442, only version v1 (2 Feb 2010), TeX `Estimates2.tex` (D:<line>). Verdict per step: OK / SLIP (wording, no effect) / GAP.

## 0. Inputs

- H(B), B = 0.87338196661320752 (P:7). Lean statement: L(s, χ) ≠ 0 for Re s > 0.87338196661320752, every modulus and χ
  (Bebel PR 101; `Riemann/Deep/HPAudit/REPORT.md` l.1, l.38: AUDIT PASSED). cert.py takes B as the ball around this exact
  decimal (C:11), so the hypothesis used is exactly the proven one, not a rounding below it. OK.
- Only H(B) for χ itself is needed: if L(ρ, χ) = 0 with β < 1 − B then L(1 − ρ̄, χ) = 0 with real part > B. OK.
- Dusart D:660–670: for x > 1, Σ_{p≤x} 1/p − log log x − B_M ≥ −(1/(10 log²x) + 4/(15 log³x)); for x ≥ 10372 the
  upper bound with the same error. D:692–702: for x ≥ 2974, Σ_{p≤x} log p/p − log x − E ≤ 0.2/log x + 0.2/log²x;
  E ≈ −1.33258 22757 33221 (D:695). PROOF uses exactly these upper bounds (P:153–154) with y ≥ e^10 > 10372 and
  y ≥ e^8 = 2980.96 > 2974. OK.
  SLIP: numbering. Section 6 of the TeX has a numbered Remark (D:368) and a commented-out Corollary (D:529), so the
  two theorems are Dusart Thm 6.10 (Σ 1/p) and Thm 6.11 (Σ log p/p); PROOF.md P:153–154, P:252 and dusart_check.py call
  them 6.11 and 6.12 (6.12 is Mertens' product, D:718). Statements quoted correctly; labels off by one.
  Dusart's proofs use RH verified to height 2.445·10^12 (Gourdon); now rigorous by Platt–Trudgian 2021 (3·10^12).
- E_UP = −1.33258227573 (C:22) ≥ E = −1.3325822757332...: the right direction (E enters W with + and W with −κ'). OK.
- Davenport Ch. 12 (17)–(18) (statement level): for primitive χ mod q, q > 1,
  L'/L(s, χ) = −½ log(q/π) − ½ ψ((s + a)/2) + B(χ) + Σ_ρ (1/(s − ρ) + 1/ρ), Re B(χ) = −Σ_ρ Re 1/ρ. q is the
  conductor because χ is primitive. OK.

## 1. Weight, explicit formula, log L(1, χ) (P:36–69)

- v, V (P:38–42): V(w) = (e^{(1+h)Lw} − e^{Lw})/(hLw²); Taylor: V = 1/w + (1 + h/2)L + O(w): residue 1. OK.
- Lemma 1 (P:44): (1/2πi)∫_{(c)} e^{wy}w^{-2}dw = max(y, 0); difference of y = (1+h)L − log n and y = L − log n over hL
  gives 1 for log n ≤ L, ((1+h)L − log n)/(hL) between, 0 beyond = v(log n/L). OK.
- Lemma 2 (P:49–60): integrand −L'/L(σ + w)V(w); residues: w = 0 gives −L'/L(σ) (V residue 1; L(σ) ≠ 0 for σ ≥ 1);
  w = ρ − σ gives −m_ρ V(ρ − σ) (−L'/L has residue −m_ρ at a zero). Left contour → 0 since |V(w)| ≤ 2e^{L Re w}/(hL|w|²)
  for Re w ≤ 0. Hence Σ_n = −L'/L(σ) − Σ_ρ V(ρ − σ): the sign in P:51 is right. Trivial zeros included (s = 0 is a zero
  for even χ). OK (qualitative inputs only, standard).
- Lemma 3 (P:63–69): log L(σ) → 0 as σ → ∞, so log L(1) = −∫_1^∞ L'/L; L(σ, χ) ≠ 0 on σ ≥ 1 holds unconditionally.
  ∫_1^∞ n^{−σ}dσ = 1/(n log n). Interchange: Σ_ρ ∫|V(ρ − σ)|dσ < ∞ by Lemma 5 (nontrivial) and Lemma 7 (trivial). OK.

## 2. Zeros (P:71–113)

- Lemma 4 (P:73–77): real part of Davenport (17) at s = 1 plus (18): Re L'/L(1) = −½log(q/π) − ½ψ((1+a)/2) + Σ_ρ Re 1/(1−ρ).
  Splitting Σ Re(1/(1−ρ) + 1/ρ) is legal: Σ(1−β)/|1−ρ|² and Σβ/|ρ|² converge absolutely (order 1). Sign and ½ψ term OK.
- Involution (P:78–80): Λ(s,χ) = εΛ(1−s, χ̄) sends ρ to 1−ρ, a zero of L(·, χ̄); conjugation L(s̄, χ̄) = conj L(s, χ) sends it
  to 1−ρ̄, a zero of L(·, χ) at the same height γ, same multiplicity. So the pairing stays inside the zeros of L(·, χ),
  also for complex χ. With H(B): 1−B ≤ β ≤ B. OK.
- Lemma 5 (P:82–91): |ρ−σ|² = (σ−β)² + γ² ≥ |1−ρ|² for σ ≥ 1; |e^{(1+h)Lw} − e^{Lw}| ≤ e^{(1+h)L Re w} + e^{L Re w};
  ∫_1^∞ e^{−cL(σ−β)}dσ = e^{−cL(1−β)}/(cL) gives g_0 exactly as written; g_1 is the σ = 1 value. Both increasing in β. OK.
- Lemma 6 (P:93–107), redone: β ≥ ½, α = 1−β, c = β, u = γ²; target g(β)/(α²+u) + g(1−β)/(c²+u) ≤ K(α/(α²+u) + c/(c²+u)).
  Times (α²+u): g(β) + g(1−β)r ≤ K(α + cr), r ∈ [α²/c², 1], linear in r. r = α²/c²: K ≥ g(β)c/α + g(1−β)α/c (uses α+c = 1).
  r = 1: K ≥ g(β) + g(1−β); K(β) − g(β) − g(1−β) = (c−α)(g(β)/α − g(1−β)/c) ≥ 0. Fixed points β = ½ need K ≥ 2g(½):
  covered (K(½) = 2g(½)). Bounds g(β)β/(1−β) ≤ g(B)B/b, g(1−β)(1−β)/β ≤ g(½). OK.
- Lemma 7 (P:109–113): trivial zeros −m, m ≡ a (2), simple; V(−y) = −e^{−Ly}(1 − e^{−hLy})/(hLy²) ∈ [−e^{−Ly}/y, 0);
  T_1 = Σ_j e^{−L_1(2j+1)} (a = 0 worst), T_0 = Σ_j e^{−L(2j+1)}/L via ∫_1^∞ e^{−L(m+σ)}/(m+σ) ≤ e^{−L(m+1)}/(L(m+1)). OK.

## 3. Loop (P:115–127)

- Lemma 2 at σ = 1, weight (L_1, h_1): Re L'/L(1) = −Re P_1 − Σ_nontriv Re V_1(ρ−1) − Σ_triv V_1 ≤ −Re P_1 + κ_1 Z + T_1.
  Lemma 4 and ψ((1+a)/2) ≤ ψ(1) = −γ_E (ψ increasing): (1−κ_1)Z ≤ −Re P_1 + ½log(q/π) − γ_E/2 + T_1. Division needs
  κ_1 < 1: hypothesis of Prop. 10; enforced in code (C:49 returns None unless the arb ball k1 < 1) and in the tail by
  Lemma T(i) (κ_1 decreasing in L at fixed h; certified κ_1(40) ≤ 9.6·10^{−18}, O:20). OK.
- (★): log|L(1)| = Re P_0 + Σ_nontriv Re G + Σ_triv G ≥ Re P_0 − κ_0 Z − T_0; insert the Z-bound (κ_0 ≥ 0). OK.

## 4. Prime sums (P:129–158)

- Lemma 8(i): p < x: v_k = 1 for p^k < x (x = e^L is transcendental, no prime power equals it); Re −log(1 − z/p) ≥
  −log(1 + 1/p); corrections have k ≥ 2, size ≤ ½Σ_{k≥2,p^k>x}p^{−k}: −τ(x)/2. p ≥ x: k = 1 term ≥ −v/p; k ≥ 2 terms
  ≤ 1/(2p(p−1)), Σ_{p≥x} ≤ 1/(x−1)·½. OK.
- Lemma 8(ii): Re w/(1−w) ≥ −r/(1+r) on |w| ≤ r = 1/p (Möbius image of the disc), so ≥ −log p/(p+1); corrections
  ≤ L_1τ(x_1) (log p < L_1); p ≥ x_1, k ≥ 2 terms vanish unless log p < (1+h_1)L_1/2. OK.
- Lemma 9(a),(b): v(log p/L) = (1/h)meas{t ∈ [1,1+h]: p < e^{Lt}}; average of Dusart's bound; (1/h)∫_1^{1+h} log t dt = E_h,
  (1/h)∫ t dt = 1 + h/2; ε_D and 0.2/y-terms decrease. OK.
- Lemma 9(c): 1/p − log(1+1/p) = [1/p + log(1−1/p)] − log(1−1/p²) > 0, total (M − γ_E) + log ζ(2) = 0.181982…;
  tail ≤ Σ_{p≥x}1/(2p²) ≤ 1/(2(x−1)). OK.
- Lemma 9(d): Σ_{m≥M}m^{−k} ≤ M^{−k}(1 + M/(k−1)), M_k = ⌊y^{1/k}⌋+1 > y^{1/k}; k = 2 gives y^{−1/2} + 2/y; 3 ≤ k ≤ K gives
  ≤ (K−2)(y^{1/3}+3)/(2y); k > K: Σ_{m≥2}m^{−k} ≤ 3·2^{−k} (Σ_{m≥3}(2/m)^k ≤ 4(ζ(2) − 5/4) = 1.58), Σ_{k>K} ≤ 6/y.
  Total ≤ τ̄. Monotonicity: each of ℓe^{−ℓ/2}, ℓe^{−2ℓ/3}, ℓe^{−ℓ}, ℓ²e^{−2ℓ/3}, ℓ²e^{−ℓ} has the stated range. OK.


## 5. Master inequalities (P:160–190)

- Prop. 10 (P:162–171): Lemma 8(i) + 9(a),(c),(d): Re P_0 ≥ C_0 − log L − E_h − ε_D(L) − 1/(x−1) − τ̄(x)/2, since
  −(M + …) + (M − γ_E + log ζ(2)) − 1/(2(x−1)) − 1/(2(x−1)); M cancels (only M's existence is used). Lemma 8(ii) + 9(b),(d):
  −Re P_1 ≤ W − (½log(q/π) − γ_E/2 + T_1). Into (★) (κ' ≥ 0): F as written. E enters W with +: E_UP is the right side. OK.
- Code = formula: C:42–60 matches Prop. 10 term by term (main_err = ε_D + 1/(2(x−1)) + τ̄/2 + 1/(2(x−1)) + T_0). My own
  implementation from PROOF.md alone (lwcheck.py, mpmath iv, 256 bits) reproduces the cert's tightest-cell margins to
  15 digits (lwcheck_B.log (1)). OK.
- F is decreasing in ell = log q: ell enters only through −κ'·ell/2, κ' > 0. So the cell rule "evaluate at the right end,
  compare with c/Λ_left" is sound. OK.
- Prop. 11 (P:175–184): Re Σ_k z^k v_k/(kp^k) ≤ −log(1−1/p) (p < x), ≤ v/p + 1/(2p(p−1)) (p ≥ x); Σ_p(−log(1−1/p) − 1/p)
  = γ_E − M > 0 termwise; −Re P_1 ≤ Σ v_1 log p/p + Σ_p log p/(p(p−1)) = Σ v_1 log p/p − γ_E − E; E cancels against
  Dusart's E. Σ_triv G ≤ 0 (Lemma 7). F_up as written; C:105–117 matches. Fup increases in ell. OK.

## 6. Theorem 1 (P:192–233)

- Λ ≤ 3 (P:196): one evaluation at ell = e^3, F decreasing in ell, so every q ≥ 3 with log q ≤ e^3 (also q = 3, 4,
  where log log q < 0.1) gets |L(1, χ)| ≥ e^{−3.287948} = 0.03733 > 0.11/3. Own optimiser: same value. OK.
- Grids (P:197–201): cells [i/1000, j/1000] contiguous from round(Λ_lo·1000) to round(Λ_hi·1000) (C:144–162); LamB is the
  float 3.13 turned into the integer 3130 by round(), so the scan starts exactly at the rational 3.13; 4.0 → 4000; 3.0 → 3000.
  No double-above-decimal gap. Tails start exactly at the grid ends 40 and 100 (arb(40), arb(100)). OK.
- Lemma T (P:204–228), checked: (i) g_1 is a sum of e^{−cL}/(hL), c > 0, at fixed h, β: κ_1(L) ≤ κ_1(Λ_N/b) < 1.
  (ii) ω(8)|_{h=2} ≤ −1.7237 (recomputed, lwcheck (0)); every term non-increasing for L ≥ 8 (Lτ̄ needs ℓ ≥ 3), ω
  increasing in h. (iii) κ_0L, κ_0e^Λ = κ_0e^{bL} decrease ((1+h)/2 − b > 0, ½ − b > 0), 1/(1−κ_1) decreases; ε_D, 1/(x−1),
  τ̄, T_0 decrease. So Err(Λ) ≤ Err(Λ_N). Requires L = Λ/b ≥ 10, L_1 = L ≥ 8: Λ_N/b ≥ 10 asserted (C:178). OK.
- Constants recomputed: κ_1(40, h = 0.015) ≤ 9.578·10^{−18}, c ≥ 0.115659; (100, 0.006): 8.39·10^{−44}, 0.116424;
  (1000, 0.001): 0.116873; (10^4, 10^{−4}): 0.116933. By hand: κ_1(40) ≈ (B/b)(e^{−40} + e^{−40.6})/(hL) = 9.58·10^{−18};
  Err(40) ≈ κ_0e^{40}/2 = 3.55·10^{−3}; c_40 = 0.92356·0.126618·e^{−0.00746−0.00355} = 0.11566. OK.
- Limit (P:230–233): h = 1/Λ_N gives hL = 1/b; κ_1 ≤ 2Be^{−Λ_N} + 2be^{−Λ_N/(2b)}; κ_0e^{Λ_N} ≤ 2/(b·hL²) = 2b/Λ_N;
  κ_0L ≤ 2e^{−Λ_N}; E_h ≤ h/2. Hence liminf ≥ e^{C_0}b = (π²/6)e^{−γ_E}(1 − B) = 0.1169398360… (recomputed). OK.

## 7. Corollary 2 (P:235–241)

- χ_{−d} primitive, real, odd, conductor d; w = 2 for d > 4; h(−d) = w√d L(1)/(2π) = √d L(1)/π. 0.11/π = 0.0350141…
  ≥ 0.0350; 0.115/π = 0.036606 ≥ 0.0366. OK.
- Thresholds: cert cells in ℓ = log d with F at log q = ℓ_{i+1} and √d ≥ e^{ℓ_i/2} (C:254–262); printed decimal threshold
  rounded up (C:272). Independent cells from floor(10^4 log d_0)/10^4 ≤ log d_0 up to the point where (A) takes over
  (log ℓ ≥ 3 and 0.11e^{ℓ/2}/(π log ℓ) > H, increasing since ℓ log ℓ > 2): all five thresholds PASS (lwcheck_B.log (5));
  7/8 thresholds 3.047·10^10, 8.137·10^13 PASS (lwcheck_78.log (5)). OK.

## 8. Side result (P:24–28, 186–190)

- T' (upper tail) as Lemma T with Err_up; bracket −γ_E + 0.2/L + 0.2/L² − ½log π − γ_E/2 + T_1 ≤ −1.40973 at L = 8
  (recomputed). Grid [3, 40] for 15 Λ: independent adaptive cover PASS (worst margin 0.0192 at [3, 3.001]); tails
  14.2223 (40), 14.1288 (100), 14.0746 (1000), 14.0673 (10^4); limit e^{γ_E}/b = 14.06650. OK.
