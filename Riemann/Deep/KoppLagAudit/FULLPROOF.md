# FULLPROOF: chain of unit-generated-orders (Kopp–Lagarias Problem 1.7), refereed line by line

Auditor ugoaudit, 19 Nov 2026. Sources: KL = Kopp–Lagarias arXiv 2512.11311 v3 (TeX `CLNO-Unit-generated-arxiv-apr-21-2026.tex`, 1767 lines,
latest version per the abs page); D = Riemann/Deep/KoppLag/DERIVATION.md (line numbers "D l.N"); L = discoveries/littlewood-half-plane.md,
audited in Riemann/Deep/LittlewoodAudit/REPORT.md (AUDIT PASSED, item (a)).

## 0. Statement versus source

- Problem 1.7 = KL \begin{prob}\label{prob:2}, l.329–335 (page: l.329–336): "Determine all unit-generated real quadratic orders O_Δ having one class per genus,
  that is, having Cl⁺(O_Δ) = Cl⁺(O_Δ)[2]. More generally, determine all such orders having Cl(O_Δ) = Cl(O_Δ)[2]." Matches the page verbatim.
- KL Prop. 2.3 (prop:23, l.543–556): O = ℤ[ε] unit-generated ⇒ Δ = f²Δ₀ = n² ∓ 4 (N(ε) = ±1); conversely Δ = n² − 4 (n ≥ 3) or n² + 4 (n ≥ 1)
  gives a unit-generated order. Matches the page's parametrisation (D l.7–9).
- "conjecturally complete" = KL l.145 ✓; "seems out of reach at this time, in parallel with the idoneal number problem being unsolved" = KL l.1743 ✓;
  "a total of 86" = KL l.1609 ✓. KL Thm 1.5 (l.303–305) is the finiteness statement, ineffective via Hua (l.307) ✓.
- Tables 2A (l.1629–1668) / 2B (l.1682–1704), parsed by our parse_tables2.py: 2A has 55 entries (53 positive + Δ = −4, −3), 2B has 37;
  90 positive entries, 89 distinct (Δ = 5 in both), 60 with f = 1, 29 with f > 1, max Δ = 68640 (2A l.1661: Cl (ℤ/2)⁴, Cl⁺ (ℤ/2)⁵, f = 2);
  Cl⁺ not 2-torsion exactly for 221, 1517 (ℤ/4, f = 1), 725 (ℤ/4, f = 5), 3965, 7565 (ℤ/2 × ℤ/4): 84 distinct one-class-per-genus orders.
  Every count in the claim matches the tables; the side remark (86 ≠ 89) holds.
- Genus theory: KL l.1525–1537 define μ(Δ) = r (Δ ≡ 1 mod 4 or 4 mod 16), r + 1 (Δ ≡ 8, 12 mod 16 or 16 mod 32), r + 2 (Δ ≡ 0 mod 32),
  r = #odd primes dividing Δ; KL Prop. 3.4 (prop:genusorder, l.1543–1550) = Halter-Koch Thm 6.5.2: genus group order 2^{μ−1};
  KL l.1509–1516: 𝔉_Δ ≅ Cl⁺(O_Δ) (Halter-Koch Thm 6.4.2), genus group = Cl⁺/Cl⁺², Cl⁺ ↠ Cl with kernel of size 1 or 2.
  This agrees with Cox's table (Δ = 4m: m ≡ 1 (4): r; m ≡ 2, 3 (4) or 4 (8): r + 1; m ≡ 0 (8): r + 2). μ − 1 ≤ ω(Δ) (KL l.1552) ✓:
  odd Δ: μ − 1 = ω − 1; even Δ: μ ≤ r + 2 = ω + 1, with equality only if 32 | Δ.
- F2 (class number formula for orders): KL l.836–852 states (O_K/fO_K)^× ≥ φ(f)² and (O_Δ/fO_K)^× ≅ (ℤ/fℤ)^× from √Δ = f√Δ₀ ≡ 0;
  D l.27–31 uses the exact count f² ∏_{p|f}(1 − 1/p)(1 − χ(p)/p) (standard: the p-part is p^{2e}(1−1/p)², p^{2e}(1−1/p²), p^{2e}(1−1/p) for p split, inert, ramified, i.e. p^{2e}(1−1/p)(1−χ(p)/p)) ✓.

## 1. Lemma A (D l.46–68) — holds as written

Statement (D l.46): n ≥ 1, Δ = n² − 4N ≥ 100, α ∈ O_Δ primitive, m = |N(α)| < √Δ/2 ⇒ m = 1. Checked step by step:
- F1 (D l.22–25): x = 2u + vΔ, y = v; Δ ≡ n² ≡ n (mod 2) so x ≡ ny (mod 2); primitive ⇔ gcd(u, v) = 1; [O : αO] = |N(α)| ✓.
- Normalisation (D l.48): (αη)' = α'η', |ε/ε'| = ε² (|εε'| = 1, ε > 1 for n ≥ 1); choose j with 1 ≤ |α/α'|ε^{2j} < ε², sign with α > 0;
  primitivity and |N| are unit-invariant ✓. α² = |α/α'|·m ⇒ √m ≤ α < ε√m, |α'| = m/α ≤ √m ✓.
- y√Δ = α − α' ≥ 0; y = 0 ⇒ α ∈ O ∩ ℚ = ℤ, primitive ⇒ α = 1 ✓. y ≥ 1: x = α + α' ≥ 0 ✓; (A1) y√Δ ≤ α + |α'| < (ε + 1)√m ✓;
  (A2) x = y√Δ + 2α' ≥ y√Δ − 2√m and n + √Δ = 2ε ✓; (A3) z = x − ny even, z(x + ny) = x² − n²y² = 4N(α) − 4Ny², x + ny ≥ 0 ✓.
- (D l.56–58) φ(y) = yε − y² − m − √m > 0 ⇒ x + ny > 2(m + y²) ⇒ |z| < 2 ⇒ z = 0 ⇒ m = y², α = yε ⇒ y = 1 ✓. φ concave, so the endpoints
  y = 1 and y = Y suffice; Y < 1 leaves no y ✓.
- Constants (D l.59–64), re-derived: n² = Δ + 4N ≥ Δ − 4; (√Δ − 4/√Δ)² = Δ − 8 + 16/Δ ≤ Δ − 4 for Δ ≥ 4, so ε ≥ √Δ − 2/√Δ ✓;
  n ≤ √(Δ + 4) < √Δ + 2/√Δ, so ε < √Δ + 1 ✓; √m < Δ^{1/4}/√2 ✓. φ(1) > √Δ/2 − 1 − 2/√Δ − Δ^{1/4}/√2 = 1.5639 at Δ = 100,
  derivative in s = √Δ is ½ + 2/s² − 1/(2√2·√s) > 0 for s ≥ 10 ✓. φ(Y)/√m identity ✓ (Yε = (ε+1)ε√m/√Δ, Y² = (ε+1)²m/Δ);
  bound 10 − 0.4 − 2.2361·2.44 − 1 = 3.1440 at Δ = 100 and increasing ✓.
- Independent numerical check (ours, ugo.c "class"): Lagrange's form of Lemma A (every |a| < √Δ/2 on the principal ρ-cycle is 1) on all
  11998 orders n ≤ 6000 with Δ ≥ 100: 0 failures (logs/class6000.log). Sharpness remark (√Δ/2 not doubleable) not re-checked.
- Reviewer subagent (Lemma A and Corollary B pasted inline, verdict only): CORRECT; it recomputed 3.8 − √5 = 1.5639… and 8.6 − 2.44√5 = 3.1440…

## 2. Corollary B (D l.71–79) — holds as written

- (Δ/p) = 1 ⇒ p ∤ Δ ⇒ p ∤ f ⇒ O_p maximal, pO = 𝔭𝔭̄, 𝔭 ≠ 𝔭̄ invertible of norm p ✓ (𝔭 = pℤ + ((b + √Δ)/2)ℤ, b² ≡ Δ (4p), p ∤ b).
- 𝔭² primitive: 𝔭² ⊆ kO ⇒ k² | p² ⇒ k = p ⇒ 𝔭² ⊆ 𝔭𝔭̄ ⇒ 𝔭 ⊆ 𝔭̄ ⇒ 𝔭 = 𝔭̄ (equal norms), impossible ✓.
- [𝔭]² = 1 in Cl(O) ⇒ 𝔭² = αO with α ∈ 𝔭² ⊆ O primitive, |N(α)| = [O : 𝔭²] = p² ✓; p² < √Δ/2 ⇔ 4p⁴ < Δ ✓; Lemma A ⇒ p² = 1, contradiction ✓.
- 2 never splits: n odd ⇒ Δ ≡ 5 (8); n even ⇒ 4 | Δ ✓. Independent check: no wide-2-torsion order with Δ ≥ 100 and a split p, 4p⁴ < Δ,
  among all 11998 orders n ≤ 6000 (assert in ugo.c, never fired).


## 3. The effective bound D₀ (D l.81–106) — holds as written

- F2 (D l.26–32): h(O) = h(O_K)·f∏_{p|f}(1 − χ(p)/p)/[O_K^× : O^×], [O_K^× : O^×] = R(O)/R(O_K) (both groups ±η^ℤ), h(O_K)R(O_K) = (√Δ₀/2)L(1, χ)
  ⇒ h(O)R(O) = (√Δ/2)L(1, χ)∏_{p|f}(1 − χ(p)/p) ✓ (normalisation: R = log of the fundamental unit, χ = (Δ₀/·) primitive mod Δ₀ ≥ 5, even).
  Our check: h from ρ-cycles against √Δ·L(1, χ)∏/(2R) with our own L(1, χ) = −(2/√q)Σ_{a<q/2} χ(a) log sin(πa/q) and R = log ε_n
  (ε_n is the fundamental unit of ℤ[ε_n] except n = 3, N = +1, where R = log φ): 889 orders (all n ≤ 300, the 82 survivors, every n > 1000 with 47 | n), max |Δh| = 1.7·10⁻¹¹.
- R(O) ≤ log ε_n (ε_n > 1 is a unit of O) ✓; 2 log ε_n < log(Δ + 4): N = +1, ε_n < n = √(Δ + 4); N = −1, ε_n < √Δ ✓ (strict).
- (L) at q = Δ₀ ≥ 5 ≥ 3 (range ✓; L(1, χ) > 0 for real χ so |L| = L ✓), and max(3, log log ·) is nondecreasing, Δ₀ ≤ Δ ✓.
  (L) itself = littlewood-half-plane item (a), audited PASSED (LittlewoodAudit/REPORT.md l.1, table row (a)); natural logarithms throughout.
- ∏_{p|f}(1 − χ(p)/p) ≥ φ(f)/f ≥ 1/P_{ω(f)} ≥ 1/P_{ω(Δ)} (j-th prime of f ≥ p_j, p/(p−1) decreasing, P_k increasing) ✓. Hence (B1), strict ✓.
- 3.3 (D l.88–90): d log g/dΔ > (1/2 − 1/2.99 − 1/60)/Δ = 0.1489/Δ > 0 on [16, ∞) (log 20 = 2.9957; log log Δ > 3 ⇒ log Δ > e³ > 20) ✓; g continuous at the switch ✓.
- 3.4 (D l.91–101): Cl = Cl[2] ⇒ h = |Cl/Cl²| ≤ |Cl⁺/Cl⁺²| = 2^{μ−1} ≤ 2^{ω(Δ)} (Cl⁺ ↠ Cl, KL l.1513–1516; §0 for μ − 1 ≤ ω) ✓, so g(Δ) < T_k = 2^k P_k
  with Δ ≥ Π_k = p₁⋯p_k ✓. Exclusion test "Π_k ≥ D_k⁺" with g increasing and D_k⁺ ≥ 16 ✓.
- Induction (D l.95–99): log(Π_{k+1} + 4) ≤ log p + log(Π_k + 4) (p ≥ 1), log p ≤ log(Π_k + 1) ≤ X + 1, Y := log(Π_k + 4) ≥ X and
  Y/(Y + X + 1) ≥ X/(2X + 1); log log Π_{k+1} ≤ log(2X + 1); hence g(Π_{k+1})/g(Π_k) ≥ √p·X/(2X+1)·log X/log(2X+1) ✓ (log log Π_k > 3 for k ≥ 15,
  so the max is the log log branch) ✓. Monotone in p and X: d/dX[log X/log(2X+1)] > 0 ⇔ (2X+1)log(2X+1) > 2X log X ✓. T_{k+1}/T_k = 2p/(p−1) ≤ 2·53/52 ✓.


## 4. Computation (D l.108–127) — code audit and independent recomputation

Code audit (read line by line; all decisions exact integers or arb balls):
- cert.py Part 1: g_lower (l.35–42) takes the upper end of max(3, log log D) when the ball straddles 3 (smaller g, safe) ✓; Dplus (l.55–66)
  keeps lower(g(hi)) > upper(T) and returns the bracket end, valid by monotonicity (§3) ✓; KMAX = 60, base K = least K ≥ 15 with every
  k ∈ [K, 60] excluded (l.83–88) ✓; step constant asserted in arb (l.100) ✓; D_SIEVE = max(D₀, 10¹⁸) (l.103) ✓; robustness bisection (l.105–114)
  is valid because D₀(c) is nonincreasing in c ✓.
- cert.py Part 2: nmax (l.148–149) is the least n with Δ_n ≥ 10¹⁸ (N = +1: 10⁹ + 1; N = −1: 10⁹), range [1001, nmax] in contiguous chunks
  (l.150–164) ✓; chunk primes satisfy 4p⁴ < Δ at the chunk's smallest n, hence for every member (l.155) ✓; residue tables by Euler's criterion
  (l.122–129) ✓; every chunk survivor re-tested with all primes 4p⁴ < Δ incl. p = 2 (full_ok l.138–144, guard for a short prime list) ✓.
  No gap at either end: n ≤ 1000 exact, n ≥ 1001 sieved (lesson 277 ✓).
- cert.py Part 3: all n ≤ 1000 from n = 3 (N = +1) and n = 1 (N = −1), 998 + 1000 = 1998 orders, plus the 82 survivors = 2080 (l.185) ✓.
- forms.py: reduced-form test and enumeration (l.24–28, 89–106) use the exact integer equivalents of the real inequalities (√Δ irrational) ✓;
  primitive forms only (gcd test l.102) ✓; ρ = Cohen Def. 5.6.4 (l.30–45) ✓; cycles closed and disjoint (asserts l.151, 154) ✓;
  h = h⁺ or h⁺/2 by J = [(−1, b₀, ·)] on the principal cycle (l.162) ✓; the composition test of C² and the composition-free inverse test
  [f]⁻¹ = [(c, b, a)], J[f]⁻¹ = [(−c, b, −a)] are asserted to agree (l.172–173) ✓; group invariants from |G[2^k]| (l.192–237) ✓.
- cert.py Part 4: its transcription KL_A/KL_B (l.226–241) equals our own parse of the TeX tables (compare.py) ✓.
- sievecheck.c: n ∈ [1001, 10⁹] for both N (stops when Δ ≥ 10¹⁸ + 8), 128-bit bounds, QR tables by squaring ✓. checks.py, PARI scripts:
  cross-checks only (C1–C8, quadclassunit, lfun); read, not load-bearing.

Independent recomputation (our code, ugo.c in C + d0check.py with mpmath intervals; nothing reused from KoppLag/):
- D₀ (d0check.log): c = 0.11 gives D₀ = 20485358601780401 exactly, worst k = 14, possible k = 1..14, base K = 15; step constant
  3.0220584767… > 2·53/52 = 2.03846 (log Π₁₅ = 40.9602194…, Π₁₅ = 614889782588491410). c = 0.037: D₀ = 9.6959·10¹⁷ ≤ 10¹⁸, base K = 16;
  c = 0.036: D₀ = 1.0279·10¹⁸ > 10¹⁸ (so 0.037 is the 3-decimal threshold of this method, as claimed).
- Sieve (sieve.log, 190 s, 33.8 MB): Legendre tables by Euler's criterion, early exit, n ∈ [1001, 10⁹ + 1], both N: 58 + 24 survivors,
  identical to KoppLag's survivors.txt.
- Class groups (class6000.log, 33 s, 39 MB): all 11998 orders with n ≤ 6000 (six times the exact range of cert.py; covers every survivor),
  reduced primitive forms, ρ-permutation (in-degree 1 asserted), h⁺ = #cycles, h from J; 2-torsion decided WITHOUT composition by genus
  theory: Cl⁺ = Cl⁺[2] ⇔ h⁺ = 2^{μ−1}; Cl = Cl[2] ⇔ h = |Cl/Cl²| = 2^{μ−1} or 2^{μ−2} according as J ∈ Cl⁺² or not, and J ∈ Cl⁺² ⇔ every
  assigned character is +1 at −1 (odd p | Δ: needs p ≡ 1 (4); 2-adic: none for Δ odd or Δ ≡ 4 (16); ε, which is +1 at −1, for Δ ≡ 8 (32);
  in every other even case δ or δε occurs, which is −1 at −1; Cox Thm 3.15 translated to Δ = 4m).
  Asserted on every order: 2^{μ−1} | h⁺; J = 1 ⇒ J ∈ Cl⁺²; μ equals KL's mod-16/32 table. Result: exactly 90 entries (89 orders) are
  wide-2-torsion, 85 entries (84 orders) narrow; with Cl ≅ (ℤ/2)^{log₂h}, Cl⁺ ≅ (ℤ/2)^{μ−1} or (ℤ/2)^{μ−2} × ℤ/4 (then h⁺ = 2^μ, asserted)
  and f, all 90 match the TeX tables (compare.py). n = 262: f = 2, Δ₀ = 17160, h = 16, h⁺ = 32, μ = 6. All survivors: same h, h⁺ as
  survivors.txt (h from 21 to 192), none 2-torsion. No 2-torsion order outside the list for 262 < n ≤ 6000.
- Analytic class number formula (own L(1, χ) by the log-sine sum, double precision): 889 orders (n ≤ 300, all survivors, 47 | n > 1000),
  max |h_an − h| = 1.7·10⁻¹⁰ (final run `ugo class 6000 300 47`, 214 s, 37 MB; same TORS/SURVCLASS lines as the first run).
- Lemma A at large n (pcyc.log): principal cycle of 4008 orders (2000 pseudo-random n ∈ [10⁶, 10⁹] per family, and n = 10⁹, 10⁹ − 1,
  1001, 1000): always the 2-cycle, no |a| ∉ {1} below √Δ/2.
- cert.py rerun on a copy (.scratch/ugoaudit/certcopy, 137 s, memcap peak 149 MB): cert_output.txt identical after masking timings
  (logs/cert_rerun_output.txt), survivors.txt identical.


## 5. Assembly (D l.129–134) and dependents

- Let Cl(O_Δ) = Cl(O_Δ)[2]. §3 ⇒ Δ < D₀ = 20485358601780401 < 10¹⁸ (needs (L) only here). n ≤ 1000: decided by the exact class groups.
  n > 1000: Δ > 10⁶ ≥ 100, so Corollary B applies and n is a sieve survivor (n ≤ 10⁹ covers Δ < 10¹⁸ in both families: N = +1 needs n ≤ 10⁹,
  N = −1 needs n ≤ 10⁹ − 1); no survivor is 2-torsion. Hence the list is the n ≤ 1000 list = Tables 2A/2B (all with n ≤ 262) ✓.
- (b): Cl⁺ ↠ Cl, so Cl⁺ = Cl⁺[2] ⇒ Cl = Cl[2]; Cl⁺ is computed exactly on the list ✓. (c): §3 unused below 10¹⁸ ✓. (d): §3.5, d0check ✓.
- No unit-generated order is counted twice except Δ = 5: n₁² − 4 = n₂² + 4 ⇔ (n₁ − n₂)(n₁ + n₂) = 8 ⇔ (n₁, n₂) = (3, 1) ✓.
- Conditional structure: if (L) were false, (c) (complete list for Δ < 10¹⁸, unconditional) survives; any c ≥ 0.037 restores the full claim.
- Dependents (rg 'unit-generated|kopplag' over discoveries/, TREE.md, COROLLARIES.md, RECORDS.md, DISCOVERIES.md; 'Kopp' over Riemann/Tree, NOGO.md):
  DISCOVERIES.md l.11 (half-plane consequences line), TREE.md l.15 (N11 entry), TREE.md l.450 + Riemann/Tree/R414.lean (round authenticating
  the inputs). No page in COROLLARIES.md, RECORDS.md, NOGO.md or another discoveries/ file builds on it. All stand.

