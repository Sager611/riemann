# FULLPROOF — unit-group-cover, the chain with source line numbers (unitaudit, 19 Nov 2026)

D:n = Riemann/Deep/UnitCover/DERIVATION.md line n; C:n = UnitCover/cert.py line n; E:n = Riemann/Deep/ExplProg/cert.py
line n; page = discoveries/unit-group-cover.md. ✔ = checked by me line by line; ✔c = also recomputed by my own code.

## 0. Inputs (D:8–22)
- (H) half-plane: HPAudit item (B) (Bebel, every modulus and χ, Re s > 0.87338196661320752, includes q = 1 ⇒ ζ; ζ also from
  (A), Re s > 7/8). B = 0.8733819666132076 ≥ 0.87338196661320752 ✔; B − B₂ = 8.0·10⁻¹⁷ > 10⁻¹⁷ (C log line 1) ✔.
  Functional equation: zeros of primitive L(s,χ) symmetric under ρ ↦ 1 − ρ̄ (same L, same γ) ⇒ b ≤ β ≤ B ✔.
- ExplProg kernel, Lemma Z, Zζ, K1, (1.1): audited PASSED (ExplProgAudit/REPORT.md l.1), imported unchanged (C:17–19) ✔.
  Same usage as E:229–237 (omega_parts: Ω ≤ (x^B Z + x^{1/2}Z_L)/2, Ω₁ ≤ x^B z₁ + x^{1/2}z₁^L) ✔.
- Σ_ρ 1/γ² ≤ 2(1 + B²/196)Σ Re 1/ρ (E:195–196): Re 1/ρ = β/(β²+γ²) ≥ (β/γ²)/(1 + B²/196) for |γ| ≥ 14, β ≤ B; pair
  β + (1 − β) = 1 ⇒ Σ_ρ 1/γ² ≤ 2(1 + B²/196)·0.0230957 = 0.046371186 ✔ (needs no ζ zero with |γ| < 14: classical).
- RS62 Thm 9 θ(x) < 1.01624x ✔ (C:27). Elementary sums ✔ (see §3.4 below).

## 1. Construction (D:24–34)
- g_j = 1 (no-op) if p_j | q or p_j repeated; S_j = S_{j−1} ∪ g_jS_{j−1} = {∏_{i∈F} g_i}; Lemma 1.1 ✔ (dropping no-ops leaves
  distinct primes coprime to q). q = 2 trivial ✔. One draw serves all q ≤ Q; union bound < 1 ⇒ a good outcome exists ✔.

## 2. Lemma D (D:36–54; C:58–77)
- Conditional law: p_j ⟂ F_{j−1}, so E[ψ(g_j)|F_{j−1}] = Σ_{p∉A_j}ν(p)ψ*(p) + ν(A_j) (ψ(p mod q) = ψ*(p) for p ∤ q) ✔;
  Re E ≤ |Σ_p νψ*| + 2ν(A_j) ✔; ν(A_j) ≤ ‖w‖_∞((j−1) log y + log q)/M (ν(p) ≤ ‖w‖_∞ log p/M, Σ_{p|q} log p ≤ log q) ✔.
- Prime powers: |Σ_{p^k≤y,k≥2} log p ψ*w| ≤ ‖w‖_∞(ψ(y) − θ(y)) ≤ ‖w‖_∞·1.01624(√y + y^{1/3} log y/log 2) (C:54–56) ✔
  (k = 3..⌊log y/log 2⌋ has < log y/log 2 terms, each ≤ 1.01624 y^{1/3}).
- EF for primitive ψ* (q* ≥ 3; no pole at s = 1; W entire since supp w ⊂ [e⁻⁷, 1]): trivial zeros −m ∈ {0,−2,…} or
  {−1,−3,…}: |y^{−m}W(−m)| ≤ e^{−Ym}W(0) (w ≥ 0, ty ≥ e^Y) ⇒ ≤ W(0)/(1 − e^{−2Y}) ✔. Even χ: s = 0 counted there ✔.
- Zero sum: pair {ρ, 1 − ρ̄}, β ≥ 1/2 member ≤ y^B G(γ), other ≤ y^{1/2}G_L(γ), G, G_L = Tri.ks(Y) envelopes of
  x^{β−B}|W| (β ∈ [1/2, B]) and x^{β−1/2}|W| (β ∈ [b, 1/2]) (E:220–227; 1/(eY) = sup_c c e^{−cY}) ⇒ (y^B Z + y^{1/2}Z_L)/2 ✔.
  Z(L) bounds Σ_ρ G over zeros of any primitive character of conductor ≤ e^L, L ≥ log 9 (ExplProg Lemma Z, audited) ✔.
- ζ (denominator): Σ Λ(n)w(n/y) = yW(1) − Σ_ρ y^ρW(ρ) − Σ_{m≥1} y^{−2m}W(−2m) (−ζ'/ζ regular at 0, W entire) ⇒
  ≥ yW(1) − (y^B z₁ + y^{1/2}z₁^L) − W(0)/(1 − e^{−2Y}); M ≥ that − PP ✔. δ = num/den, den > 0 certified (C:75) ✔.
- Code vs formula: C:70–77 is exactly D:40–44 (num, den, triv, PP, z₁ = Σ1/γ²·k₂/2) ✔.

## 2b. Lemma D, uniform version (D:56–74; C:79–118)
- Monotonicity (D:56–64): with the envelope fixed at the cell corner Y₀ (valid for all Y ≥ Y₀: k's = n + n₀/(eY) etc.
  decrease), n = num/y^B is non-increasing and d = den/y^B non-decreasing in y (y^{1/2−B}; y^{1/3−B} log y for log y >
  1/(B − 1/3) = 1.85; log y/y^B for log y > 1/B; checked X > 2, X > 1/B at C:83) ✔. No L-monotonicity needed: Z at the
  cell's upper end L₁ (conductors ≤ e^{L₁}, Σ_{p|q} log p ≤ L₁), y₀ at L₀, N − 1 = 4L₁ + 200 ≥ 4L + 200 (C:101) ✔.
- Cells (C:97–106): start exactly at log 9, contiguous (next L₀ := this L₁, same ball), stop at the first L₀ ≥ 60;
  the tail starts at that same Le = 60.2753 (C:109–110) ⇒ no gap ✔. 664 cells ✔c.
- Tail (D:69–74, C:79–95): num + δ₀(den-subtractions) ≤ δ₀yW(1), divided by y^B L with y^b = c'L; every term
  non-increasing in L (k's at Y_m, Z/L by the audited majorant E:259–266) ✔; algebra of each term vs code ✔.
- Exponent 7.8978 > E (page claim): y' = (820L)^{7.8978} ≥ (820L)^E (820L ≥ 1), so δ(y') ≤ δ(y) by monotonicity in y ✔c.
- ✔c (indep_delta.py, ExplProgAudit/indep.py quadrature, 256 bits): worst cell 0.00098366393 at L₀ = 3.1308685
  (claim 0.00098366), c' = 820000: 9.825817·10⁻⁷ (claim 9.8258·10⁻⁷); tails 0.344339 / 0.345212 (claim 0.3443 / 0.3452);
  10 cells per c' incl. first, last, 6 random: interior value (true L, true Y) ≤ cell bound; direct δ at L = 61…10⁴ < δ₀;
  table c' (750, 670, 600, 530, 490, 460; 600000, 490000, 460000) with the row's largest N: all δ ≤ δ₀ (tightest
  Q = 10³⁰: 0.000999961) ✔c. logs/indep_delta_*.txt.

## 3. Lemma C (D:76–126; C:121–261)
- 3.1 (D:82–86): μ̂_j = ∏(1 + χ(g_i))/2, supp μ_j = S_j, Parseval C_r = H‖μ_r‖² − 1, C–S ⇒ C_r ≥ h/(H − h) ≥ t ✔.
- 3.2 (D:87–93): expansion weights 2^{−k−|e|} sum to 1 ✔; e ∉ R(a) ⇒ ψ_e nonprincipal, Re E ≤ δ (Lemma D) ✔; tower
  with factors in [0,1] and π + δ(1 − π) ≥ 0 ✔; Markov ✔. Same as Hughes's Lemma for k = 1 (Hughes TeX l.132–177:
  ρ = (1 + δ)/2, β = δ + (1 − δ)(B+1)α, here with no bad characters B = 0) ✔ (TeX read).
- 3.3 (D:94–115): R(a) closed ✔; kernel = Hom(ℤ^k/span R, Ĝ), size H^{k−ρ}∏|G[d_i]| ✔; (3.1) ✔; BFS completeness
  (each closed R reached by adding its elements one by one; intermediate spans ⊆ span R, no unit vector) ✔.
  ✔c (indep_classes.py, own extended-gcd insertion + canonical HNF + own Smith form, BFS over lattices; lattice ↔
  closed set is a bijection): 1, 4, 48, 2340 sets; class tables (π, k − ρ, #d > 1, ∏d, count) identical to cert.py for
  k = 1..4 (57 classes, max D = 16), and the closed sets identical as sets. Brute force on groups ℤ/2×ℤ/8, ℤ/12, (ℤ/2)³,
  (ℤ/3)², ℤ/5 (all 4-tuples), 4 groups × 20000 random 4-tuples, (ℤ/7)² all 3-tuples: every R(a) enumerated and
  #{a : R(a) = R} ≤ H^{k−ρ}∏|G[d_i]| ✔c.
- 3.4 (D:116–121; C:210–215): |G_q[d]| = ∏ gcd(d, φ(p^a)) (odd), gcd(2,d)gcd(2^{a−2},d) (2^a, a ≥ 3) ≤ 2d^{ω(q)} ✔;
  D^{ω(q)} ≤ τ_D(q) ✔; Σ τ_D ≤ Q(1 + log Q)^{D−1} ✔. ✔c (falsify.py A): CRT formula = brute force (q ≤ 400, d ≤ 16);
  for all 24 (k − ρ, d-tuple) types, the exact sums Σ_{3≤q≤Q} φ(q)^m∏|G_q[d_i]| ≤ the stated bound at 17 checkpoints
  Q ≤ 4000 (max ratio 0.64).
- 3.5 (D:122–127; C:217–226): T_{j+1} = T_j ∩ g_{j+1}T_j ✔; E|T ∩ gT| formula (Hughes l.160–167 with B = 0) ✔;
  β-contraction on {|T_r| < αH} ✔; Σ_q φ(q) ≤ Q²/2 ✔; U(Q) (3.2) ✔ = code ✔.
- Reviewer (verdict-only, this chain pasted inline): passes 1–2 flagged only hypotheses my paste omitted (N = r + s,
  common draws; 0 < α < 1), both stated at D:125, D:132; pass 3: all steps CORRECT (.scratch reviewer3).

## 4. Theorem for every Q ≥ 3 (D:128–147; C:263–331)
- r ≥ a₁ log Q + a₀ with the exact a₁ = (k+1)/log(1/ρ_g) ⇒ ρ^r ≤ ρ^{a₁ log Q + a₀} ✔; generic class: Q^{k+1}ρ_g^{a₁ log Q} = 1,
  term ρ_g^{a₀}/(k+1) ✔; other classes ε = a₁ log(1/ρ) − e > 0 (certified, C:290) ✔; sup_{x ≥ log 3}(1+x)^f e^{−εx}:
  maximiser 1 + x = f/ε; C:264–271 returns the value at log 3 only when f/ε − 1 ≤ log 3 is certain, else the global
  max (f/ε)^f e^{ε−f} (≥ the sup on any half-line) ✔. Hole: Q²β^{b₁ log Q} = 1 ✔. N ≤ (a₁+b₁) log Q + a₀ + b₀ + 2 ✔.
- ✔c (indep_union.py, mpmath.iv 256-bit intervals, own class list): a₁ ≤ 2.8895567, 2.1687288, 1.9300679, 1.8131051;
  b₁ ≤ 0.29357643; sup U ≤ 0.618322, 0.46409, 0.726579, 0.546406 (k = 1..4), 0.53876 (k = 4, Q ≥ 10¹⁰⁰, a₀ = 13);
  δ = 10⁻⁶: a₁ + b₁ ≤ 2.0693703 (k = 3, sup 0.652365), 1.949149 (k = 4, sup 0.477007) — all equal to the log. With
  a₀ − 1 every sup exceeds 0.9 (the code's selection threshold, C:313); it stays < 1 for k = 2 (0.927) and k = 4,
  δ = 10⁻⁶ (0.954), so one prime could be shaved there (harmless; a₀ is a choice).
  Statement: 2.1067 ≥ 2.1066815, 58 = 56 + 0 + 2; 1/log 2 + 0.664 = 2.106695 ≥ 2.1066815; 1.9492 ≥ 1.949149, 98 =
  96 + 2; 1 − 0.5464 = 0.4536 ≥ 0.45 ✔c. Actual U(Q) at the true ceilings on 15 points log Q ∈ [log 3, 10⁵]: max 0.5298
  (Q = 9) < 0.5464 ✔c. Lemma D horizon: N ≤ 2.1067 log Q + 58 ≤ 4L + 201 ✔.

## 5. Table (D:149–161; C:333–364)
- ✔c: every printed (Q, k, r, s, α) gives interval-certified U < 1 (e.g. Q = 10³, k = 2, r = 21, s = 6, α = 1/10:
  U ≤ 0.7782), and each is locally least ((r − 1, s), (r, s − 1) give U ≥ 1). Lemma D at the row's c' and largest N ✔c
  (§2b). Prime bounds: log₁₀ y = 29.3355, 31.326, 33.325, 36.0424, 39.9028, 47.5839 (claims 29.34, 31.33, 33.33,
  36.05, 39.91, 47.59) ✔c. So 10³: 27, 10⁶: 44, 10¹²: 77, 10³⁰: 164, 10¹⁰⁰: 498, 10¹⁰⁰⁰: 4858 are certified.

## 6. Small Q (D:163–170; C:366–424)
- ✔c (indep_small.py, own numpy closure S ← S ∪ S[p⁻¹x]): first 17 primes (≤ 59) cover all 2 ≤ q ≤ 2000 and 16 do not
  (hardest 1547, 1705, 1785, 1885, 1955, 1995); the 16 least primes > 2000 (≤ 2113) cover, 15 do not (1999); Q = 100,
  300, 1000: first 11 (≤ 31), 13 (≤ 41), 15 (≤ 47); least primes > Q: 9 (≤ 139), 12 (≤ 373), 14 (≤ 1091). Counting
  bound 2^{|P|} ≥ φ(q) ✔ (log₂ 1998 = 10.964).

## 7. Code audit (cert.py)
- Every proof inequality is pos() on an Arb ball (C:75, 92, 103, 110, 260, 290, 324, 344, 359) ✔; floats only pick r, s,
  α, a₀, b₀, c', cell ends, and each pick is re-certified ✔; ceil_str rounds the ball's upper end up (+10⁻⁹) ✔.
  Rationals enter exactly (A = Fraction → Arb); B, 1.01624, BMOR constants as stated ✔.
- Coverage: L-cells from log 9 (the least L) with no gaps and the tail from the last cell end; the linear theorem's sup
  from log Q = log 3 (row "Q ≥ 10¹⁰⁰" from 100 log 10) ✔; small-Q sweep from q = 2 ✔.
- Note: Arb comparisons in Python return False when undecided; every such branch (C:102 loop, C:266–270) falls back to
  a valid (larger or contiguous) choice ✔.

## 8. Falsification (falsify.py, logs/falsify.txt)
- Tight law ν = (1 − δ)·uniform + δ·1 (Re ν̂(ψ) = δ exactly) on (ℤ/qℤ)^×, q ∈ {5, 7, 8, 9, 12, 15, 16}, k = 2, 3, 4:
  exact E[C_r^k] ≤ (3.1) bound, and the exact Pr(S_N ≠ G) (DP over subsets) ≤ the Lemma C bound: no violation.
- Smallest moduli q = 3, 4, 6 (H = 2, one real χ) are in the union (all-equal real class, π = 1/2) ✔; Q ∈ [3, 9) uses
  L = log 9 ✔; primes of P dividing q are no-ops, not exclusions (KST l.126 must exclude such moduli) ✔.
