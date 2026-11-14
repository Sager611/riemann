AUDIT PASSED: every step of the RH-free Breuillard–Varjú chain (Prop 9′, Prop 19′, Theorem 2′ in the narrowed ranges dτ⁴ ≥ C(log d)² / C(log d)^{12}) and every listed descendant (Littlewood, BV Cor 3, Eberhard, FJSS, Varjú–Xu incl. Cor 1.3 in its stated range, Bary-Soroker–Garzoni–Sodin, Blachar–Breuillard–Kozma) checks line by line against the source TeX with Cureton–Thorner's 31/32 as the only zero input; one wording fix: the discovery page's paraphrase of BV Cor 3 ("roots of unity of order ≤ N") is false as written and must read "0 or a root of unity of degree < N".

# BVIrrAudit: referee audit of discoveries/random-polynomials-without-rh.md (node N11), 14 Nov 2026

Agent bvirraudit (Orca), 20:47–21:15 UTC. Chain with source lines: FULLPROOF.md. Code: code/; outputs and memory peaks: logs/.
Sources read in TeX (.scratch, deleted after): BV 1810.13360, Eberhard 2008.01223, FJSS 2106.04049, VX 2511.04240, BGS 2502.17218,
BBK 2609.40123, BSHKP 2308.04878. Inputs as audited before: CT statement (HPAudit D), R409 (statement authentication).

## Short answer
- The claim holds as written on the page, in its narrowed ranges, with the same rates up to absolute constants.
- Where RH enters BV: only Prop 9 (l.650–700, input Grenié–Molteni l.670–675), used only in Prop 19's proof (l.1199) and in §8.1
  (l.2433–2441 part 1; l.2463 part 2; l.2467 w for ζ_Q). Our rg census of RH/Riemann/hypothes/assum and of every \ref to Prop 9/19
  confirms the BVIrr list; Prop 13/20 (mild hypotheses) feed only Thms 5–7.
- Every replacement fits the exact use: same fields (Q(a), Q(a_1..a_m), all number fields, covered by CT), same prime windows
  (p ∈ (e^X/2, e^X], X inside condiX), uniform in K and in m.
- Corrections needed: (1) page wording of BV Cor 3 (false as written, counterexamples below); (2) HPAudit's dependency table says VX's
  Legendre corollary survives without CT — wrong, it uses VX Thm 1.1, which needs CT. Neither affects the main claim.

## (1) DERIVATION.md refereed (all steps correct)
| step | source | check | status |
|---|---|---|---|
| FE, CT ⇒ 1/32 ≤ β ≤ 31/32, |ρ| ≥ 1/32; trivial zeros ≤ n_K each | Hecke; HPAudit D | Γ_R/Γ_C pole orders recomputed | correct |
| EF (BV Thm 10, l.704–743) | BV | G(ρ) = ĝ(−iρ); Λ_K = Σ A(n) log n over n^m; contour to −∞ | correct |
| ZC (LO Lemma 5.4, BV l.892–897) | BV | statement level | correct |
| Prop 9′ Steps 1–8 | D §1 | B-spline C², |sinh z| ≤ cosh Re z, sandwich, ηx main term, zero sum 5.2C_0 x^{31/32}(34+log 9/η)(log Δ+n log 12/η), trivial 1.6n, η = x^{−1/32}/4 | correct; numerics 0 violations |
| (9′) normalisation | D §2 | prime powers ≤ 2X√x, Λ_K(p^e) ≤ n_K X | correct |
| Prop 19′ (a) main sum | BV l.1198–1204, Lemma 21 | n_{K_O} = |O| ≤ d^m, log Δ ≤ 4m d^{m+1/10} ≤ 4d^m X ⇒ E_1 ≤ 9C d^{2m}X²e^{−X/32} | correct |
| Prop 19′ (b) discrepancy primes | BV l.1205–1231 | resultant bound, ≤ 23 m d e^{X/10} primes, E_2 ≤ 46 m d^{m+1} X e^{−0.9X}; RH-free | correct |
| Prop 19′ (c) margin | D §3 | uniform in m; d^{1/10} ≥ 640 log d ⇔ d ≥ 10^{48.55}; in condiX 10^{25.81} | correct (only < 1/3 needed) |
| Thm 2′ (1) | BV l.2411–2451 | X = (√c_0/10)τ√d/log d; window ⇔ dτ⁴ > (10^{10}/c_0)(log d)²; exponent of Cor exp-root-2 = X | correct |
| Thm 2′ (2) | BV l.2453–2470 | X² = c_0τ²d/(100m(log md)²), m = ⌈30(log d)²⌉; window ⇔ dτ⁴ > (10^{10}/c_0)m⁵(log md)² ≍ (log d)^{12}, upper end τ√d ≫ (log d)^7; |w−1| = O(X²e^{−X/32}); B_m ≤ 2^{m²} ≤ e^{X/100} | correct |
| BV's own "d ≥ C/τ⁴" | BV l.288, 2430, 2463 | at dτ⁴ = C the chosen X is below condiX (ratio → 0): RH-independent understatement; the page's narrowing is right | confirmed |
Exponents recomputed: e^{−X/2} → e^{−X/32} (zero sum), d^{2m}X² prefactor beaten once X ≥ 320 m log d; m⁵(log md)² from X > 10⁴m²/τ with
X ∝ τ/(√m log md); Wielandt (m ≥ 6 log d) covers BV's m (their "(30 log d)²" Bochert–Jordan quote is a misquote, as BVIrr noted).

## (2) Consequences: where GRH/ERH enters, and the replacement
| consequence (page) | GRH/ERH enters (TeX) | replacement covers that use? | verdict |
|---|---|---|---|
| Littlewood ±1: P(irr) → 1, rates, Gal ⊇ A_n | BSHKP l.113 (limit "known under GRH"); via BV Cor 3 l.333 + Thm 2(2) | Cor 3 with N = 2, U_2 = {0,±1}; Thm 2′(2), τ = 0.49 | holds as written |
| BV Cor 3, every finite-variance law | BV l.333–338 via Thm 2(1) (proof l.2494–2521) | Thm 2′(1) at fixed τ(μ) after truncation | holds with BV's U_N; page wording false |
| Eberhard Thm 1.3 | l.604–612 (BV Prop 9) in Prop local-global l.693–700; PNT for Q l.735, l.819 | m = 1: X = 32cn/33 (local bound valid for all p > 2H+1, l.712–721); m ≥ 2: X = cn^{1/2}/m^{5/2} unchanged | holds as written (same rates, smaller c) |
| FJSS Thm 1.1 | l.481–489 (BV Prop 9, Eb Prop 3.5), used l.512–521; Riemann error for Q l.514 | X = cn^{1/4} (top of local range); n³X²e^{−X/32} log n → 0 | holds as written |
| VX Thm 1.1 | Prop 2.1 l.361–458 (= BV Prop 9 at l.384–394); §6 l.972–1031, X = 5 log d | X = 81 log d, l > 81/ε; local Prop l.821–830 allows any q > d (its proof, l.879–968, needs only l ≥ ε^{−1} log q/log d) | holds as written |
| VX Cor 1.3, f > 2^{7.9π(d)} (Bebel) / 2^{32π(d)}d^{65} (CT) | GS Prop 9.1 under RH for Dirichlet L, l.1038–1050, used in the proof l.1050–1086; and Thm 1.1 | Λ_{Q(√D)} = Λ(1+χ_D); error x^T(d+log x)²; a > 1/(1−T*) = 7.89777 | holds in the stated range (d ≥ d_0); needs CT too |
| BGS Thms 1, 2, minors (l.685) | (eq:bv) l.936–947 used l.1149, l.1241, l.2043 | x = e^{c′n}; ε < c′/64 allowed by cor:height (any ε) | holds as written (smaller c) |
| BBK general-law remark (l.189) | none in proofs; remark cites BV Thm 2 | Thm 2′(2) + cyclotomic tail + BBK Thm 1.1 (unconditional) | holds (finitely supported μ, H_2 > 12, ω_0 ≠ 0) |
Not audited (outside the brief, as the page says): 2304.12075's G(n, p) transfer.

## (3) Falsification
- Cor 3 paraphrase, two counterexamples: (i) μ uniform on {0,1}, N = 2: P(x | P) = ½ so P(irr) ≤ ½, but the page's formula gives
  1 − √(2/(πd)) + O(1/d); (ii) μ = ±1, N = 3, n ≡ 3 mod 4: the omitted Φ_4 term P(f(i) = 0) ~ 4/(πn) ≫ n^{−3/2} (exact ratio 40 at
  n = 999, 403 at n = 99999; logs/cor3_misquote.txt). BV's U_N (0 and roots of unity of degree < N) is what BVIrr's REPORT/DERIVATION use.
- 31/32: the Lean statement quantifies over every K with an absolute constant; no n_K, no exceptional-zero clause; real zeros in
  (31/32, 1) are excluded, so Deuring–Heilbronn or Siegel issues cannot enter.
- Small x: the only pointwise use of an RH bound at x polynomial in d is VX's X = 5 log d; raising X to 81 log d is allowed by VX's
  local proposition (checked in its proof). Every other use has log x ≥ a power of d or n.
- Corners: smallest d ≥ D_0 (10^{48.55}, or 10^{25.81} inside condiX) is absorbed by C (τ < 1); the window's lower end forces exactly
  the narrowed size hypotheses; at τ = (C/d)^{1/4} BV's X leaves the window (BV's own understatement, RH-independent); with BV's
  constants set to 1 and τ = ½, part 1's margin is positive from d = 10^{15} and part 2's window opens at d ≈ 10^{50}
  (logs/params_rerun.txt): "d ≥ 10^{26}" is sufficient, not necessary.
- Reruns (copies, memcap 290 MB, single thread): params.py identical to params_output.txt; littlewood.gp n = 2..15 counts identical
  (peaks ≤ 7 MB footprint, logs/peaks.txt); kernel_check.py checks Steps 2, 3, 5, 7 numerically (80,000 samples, 0 violations).

## (4) What holds
- Holds as written: the main claim (both parts, narrowed ranges, same rates), Littlewood, Eberhard, FJSS, VX Thm 1.1 and Cor 1.3
  (stated range), BGS, BBK.
- Needs a wording fix: the page's BV Cor 3 bullet → "P(irreducible) = 1 − P(P has a root that is 0 or a root of unity of degree < N)
  + O(d^{−N/2})".
- Fails: nothing.

## Reviewer
One reviewer subagent on the decisive step (Prop 9′ with exact formulas, its insertion into Prop 19 uniformly in m up to
⌈30(log d)²⌉, hidden n_K or exceptional-zero dependence): VERDICT CORRECT (absolute C₂ uniform in K; uniformity in m via
2m log d ≤ 2X log d/d^{1/10} ≤ X/320; no extra term from a real zero under CT).

## For the orchestrator
- discoveries/random-polynomials-without-rh.md, consequences bullet 1: replace "a root in the roots of unity of order ≤ N" by
  "a root that is 0 or a root of unity of degree < N".
- Riemann/Deep/HPAudit/REPORT.md dependency table, row "Random polynomials …": "the Legendre corollary survives with HB" is wrong
  (VX Cor 1.3 rests on VX Thm 1.1, i.e. on CT); if HCT fails the whole row falls.
- Source-paper slips met (harmless): BV Prop 9's error vanishes for K = Q as printed; FJSS l.519, l.521 "exp(cn^{−1/4})" for exp(−cn^{1/4});
  BGS (eq:bv) has n^k where the orbit count gives n^{2k}; BGS (eq:rh) covers only Q(a) though k = 6 needs Q(a_1..a_6).

RH relevance: N11
