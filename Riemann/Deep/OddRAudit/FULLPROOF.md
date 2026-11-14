# FULLPROOF: the odd Robin chain with source line numbers (oddraudit, 14 Nov 2026)

Sources. [P] = Riemann/Deep/OddRobin/PROOF.md (line numbers l.); [W] = Washington–Yang, arXiv 2008.04787v1 TeX `WashYangRobinLag.tex` (the only arXiv version, 11 Aug 2020; journal version IJNT 17 (2021) 843–870 not accessible, not compared); numbering in [W] uses one shared counter per section (TeX l. 31–40), so its main theorem MainThm1 (l. 87) is Theorem 1.2 (LEAN.md calls it 1.1); Lemma CAnum = 3.4 (l. 383), robinProp = Prop. 3.5 (l. 390), Theorem 3.6 (l. 440), AlphaLemma = 2.2 (l. 140), xkLemma1 = 2.5 (l. 168), xkLemma2 = 2.6 (l. 185).
Notation. (OR): σ(n)/n < (e^γ/2) log log n, odd n ≥ 3. N = 3^4·5^3·7^2·11···61, c0 = 67N.

## C0. Definitions ([P] l. 9; [W] l. 118–197)
- Odd CA ([W] Def. 2.1, l. 118–132): strict max of g(k) = σ(k)/k^{1+ε} over odd k < N, ≥ over odd k > N. At a tie ε the smaller maximiser is the CA one.
- Exponents ([W] Lemma 2.2, l. 140–147): α_p(ε) = ⌊log((p^{1+ε}−1)/(p^ε−1))/log p⌋ − 1, ε ≤ log_3 4 − 1. Checked: at ε = F(p,a), p^ε = 1 + 1/(p+…+p^a), so (p^{1+ε}−1)/(p^ε−1) = p^{a+1} exactly and the floor formula gives α_p = a, the LARGER maximiser. [W]'s proof (l. 153) says "always taking the smaller power"; that contradicts its own formula at ties, harmlessly: [P] l. 9 adopts the formula and notes that every odd CA number is N_ε on an open ε-interval, so the set of odd CA numbers is the same. ✓
- F(x,k) ([W] Def. 2.3, l. 157), x_k ([W] Def. 2.4, l. 162), α_p = k for x_{k+1} < p ≤ x_k ([W] Lemma 2.5, l. 168–183, proof checked: it is the floor formula rearranged). Increment identity [P] l. 9: σ(p^a) − pσ(p^{a−1}) = 1 ✓, so the exponent of p rises a−1 → a exactly at ε = F(p,a) ✓.
- f_{2,1}(n) = σ_{2,1}(n)/(n log log n), σ_{2,1} = σ of the odd part ([W] l. 69–71 with φ(2) = 1) ✓ = [P] l. 81.

## C1. The gap N < n < c0 ([P] l. 13–41; gap.py)
- (1a) l. 17–18. ε* = log(68/67)/log 67 = F(67,1) ✓ (F(67,1) = log(1+1/67)/log 67). Increment log g_p(a+1) − log g_p(a) = log(σ(p^{a+1})/σ(p^a)) − (1+ε) log p and σ(p^{a+1})/σ(p^a) = p + 1/σ(p^a) decreases in a ✓ ⇒ log g_p concave in a, c_p(a) ≥ 0, monotone away from α_p ✓. gap.py checks strict maximality at ε* for every odd p ≤ 997, p ≠ 67 (tie at 67 to 1e−60) ✓; p ≥ 73 by (1e) ✓.
- (1b) l. 20–22. log(σ(n)/n) = Σ_p [log g_p(a_p) + ε a_p log p]; subtract the same for N ✓ (exact identity, any ε).
- (1c) l. 24–26. violation ⟺ log(σ(n)/n) ≥ log(e^γ/2) + log log t, t = log n ⟺ L(n) ≤ D(t) ✓; δ0 = N's log-margin.
- (1d) l. 28. D'(t) = ε* − 1/(t log t) < 0 on [t_N, t_c0] since ε* = 0.0035234646 < 1/(t_c0 log t_c0) = 0.0036790869 (recomputed, below) ✓ ⇒ a violator n ∈ (N, c0) has L(n) ≤ D(t) < D(t_N) = δ0 ✓.
- (1e) l. 30. α_p = 0 for p ≥ 73; c_p(a) ≥ c_p(1) = (1+ε*) log p − log(p+1), d/dp = (1+ε*)/p − 1/(p+1) > 0 ✓, c_73(1) = 0.00151163 > δ0 ✓.
- (1f) l. 32–41. Exponents 0..39 for p ≤ 71 suffice: costs increase away from α_p ≤ 4 and c_p(39) ≫ δ0 ✓. gap.py takes every exponent with cost.lower() ≤ δ0.upper() and every combination with total lower bound ≤ δ0.upper() (a superset) ✓; 10 vectors, 3 in the gap ✓ (rerun and independent recount below).

## C2. Interpolation ([P] l. 43–45; [W] Lemma 3.4, l. 383–388 = Robin 1984 §3 Prop. 1)
- Consecutive odd CA N1 < N2: at the event ε2 between them (one or several (p,a) with F(p,a) = ε2) both maximise g over odd m; g(n) ≤ g(N1) = g(N2) gives log(σ(n)/n) ≤ ℓ(t) := log(σ(N1)/N1) + ε2(t − t1), and ℓ(t2) = log(σ(N2)/N2) ✓. Several primes changing at one ε change nothing: the bound uses only g(n) ≤ max g, valid for every odd n ✓ (and no two events tie below x_1 = 2000: interp.py, min relative gap printed in interp.log).
- h(t) = log(e^γ/2) + log log t is concave for t > 1 (h'' = −(1 + log t)/(t log t)² < 0) ✓; h − ℓ concave and > 0 at t1, t2 ⇒ > 0 on [t1, t2] ✓. Endpoints are odd CA numbers themselves ✓. So: an odd violator in [N1, N2] forces N1 or N2 to violate ([W] l. 384) ✓.
- Use in [P] l. 79: every odd n ∈ [c0, N*] lies in some [N1, N2] of consecutive odd CA numbers in [c0, N*] ✓ (the sequence is unbounded, increasing, starts at c0 after N, and N* is a member).

## C3. Direct check, x_1 ≤ 5·10^8 ([P] l. 47–50; oddca.py)
- Events (p,1) in increasing p ✓ (F(p,1) decreasing in p); all (p,a ≥ 2) with F(p,a) ≥ F(XD,1)(1 − 10^−9) collected from the first 20000 odd primes (asserted enough: F(p_20000, 2) < F(XD,1)) ✓; near-ties (relative gap < 10^−9) certified in arb ✓; merge positions certified by strict arb comparisons with both neighbouring (q,1) events, else FAIL ✓ (so no exact tie (q,1) = (p,a) occurs).
- Increments: S += log(1 + 1/(pσ(p^{a−1}))) = log(σ(p^a)/p^a) − log(σ(p^{a−1})/p^{a−1}) ✓, T += log p ✓; violation tested as m > 0, undecided (0 ∈ m) = FAIL ✓.
- Blocks: S and T increase along a run of (p,1) events; every number of a block has S ≤ S_after, T ≥ T_before; S_after − log(e^γ/2) − log log T_before < 0 certifies all ✓; failing blocks halve down to single exact checks ✓; block products are exact integers ✓.
- Enclosures: 128-bit arb; 26.4 M additions; final rad(T) < 10^−20 asserted ✓. The last violator is the number completed by (61,1) = N, the next event is (67,1) = c0 (interp.py, independently in mpmath) ✓.

## C4. Transfer, 5·10^8 ≤ x_1 ≤ X_R ([P] l. 52–71; transfer.py; RobinCH certify3.margins4)
- Lemma O1 (l. 54–55): σ(2^a)/2^a = 2 − 2^{−a} = 2(1 − 2^{−a−1}); taking logs, (OR) at N_ε ⟺ R(M_ε) > C_odd ✓ (signs rechecked).
- O2(i) (l. 58): k = ⌊log_2 x⌋ ≥ 2; x_k > x^{1/k} ≥ 2 ([W] Lemma 2.6(a), proof: t + … + t^k < k t^k and log(1 + u/k) ≥ log(1+u)/k) ⇒ α_2 ≥ k ([W] Lemma 2.5) ⇒ 2^{−a} ≤ 2^{−k} < 2^{1 − log_2 x} = 2/x ✓.
- O2(ii) (l. 59): ε ≤ F(2,a) (α_2 ≥ a) and F(2,a) = log(1 + 1/(2^{a+1} − 2))/log 2 < 1/((2^{a+1} − 2) log 2), ε = log(1 + 1/x)/log x > 1/((x+1) log x) ⇒ 2^a < ((x+1) log x/log 2 + 2)/2 ✓.
- O2(iii) (l. 60): all odd p ≤ x divide N_ε (floor convention) ⇒ log N_ε ≥ θ(x) − log 2; θ(x) > x(1 − 1/log x) for x ≥ 41 (Rosser–Schoenfeld 1962, Thm 4) ✓. O2(iv) (l. 61): −log(1 − u) ≤ u/(1 − u) ≤ 2u for u = 2^{−a−1} ≤ 1/2 ✓; (log log)' = 1/(u log u) decreasing ⇒ MVT bound ✓; A ≥ a log 2 and L_N ≤ log N with L_N log L_N in the denominator ✓.
- Use of m (l. 64–65): margins4('R') = robin_pos − Σ(negative parts). robin_pos = (√x_a log x_a)·(a lower bound of Σ_{x2hi < p ≤ x_a} p^{−2} + L_2 valid for every x in the box, ≥ 0) ≤ (√x log x)·(the same) for x ≥ x_a ✓; every negative part is an arb enclosure over the ball L = box(La, Lb) ✓. So m.lower() ≤ (√x log x) R(M_ε) on the whole box, granting RobinCH's analytic lemmas ✓; transfer.py requires (m − corr) > 0 with corr an upper bound of (√x log x)·C_odd over the box ✓.
- Not re-derived here: RobinCH's analytic lemmas (Lemmas 9–13, 11_k) and inputs (RH to 3·10^12, Chirre–Helfgott ψ bound = arXiv 2512.15709, checked by CHAudit; Büthe θ < 10^19); inherited as stated in [P] l. 86.


## C5. Theorem ([P] l. 73–82; [W] l. 87–96, 390–456)
- (a) l. 75–76: N < n < c0 by C1; c0 ≤ n ≤ N* by C2 on consecutive odd CA numbers, each satisfying (OR) by C3 (x_1 ≤ 5·10^8) or C4 (5·10^8 ≤ x_1 ≤ X_R) ✓. Ranges: log N_last = 500016054.5165 ⇒ e^{500016054} ≈ 10^(2.17·10^8) ✓; log10 N* ≥ (X_R(1 − η) − log 2)/log 10 ⇒ 10^(10^46.4480), 10^(10^46.7931), 10^(10^46.7900) (rerun) ✓; the 46.79 rows are conditional on the half-plane input ✓ as stated.
- (b) l. 78–79: [W] Thm 1.2 (l. 87–96; proof = Prop. 3.5 + Thm 3.6, l. 390–456): RH ⟺ (OR) for all odd n ≥ c0. (⇒) under RH (OR) holds on [c0, ∞) by [W], on (N, c0) by C1. (⇐) (OR) on (N, ∞) contains (OR) on [c0, ∞), so RH by [W] l. 452–456 (Robin: RH false ⇒ infinitely many CA M with f(M) > e^γ; M = 2^k N', f(M) ≤ (σ(2^k)/2^k) f(N') < 2·e^γ/2; checked) ✓. Sharp: δ0 > 0 ✓. [P]'s longer route (Prop. 3.5 for x_1 ≥ 120409 + C3 + C2 + C1) is also valid ✓ (F5).
- [W] Prop. 3.5's analytic part (l. 404–437, x_1 ≥ 120409, RH): σ(N)/N ≤ ∏_{x2<p≤x1}(1 − p^{−2}) ∏_{2<p≤x1}(1 − 1/p)^{−1} (l. 425–427) ✓ (α_p = 1 exactly on (x_2, x_1]); x_2 < √(2x_1) so the product over (√(2x_1), x_1] dominates ✓; constants recomputed: c = 0.0461914, α(120409) = 2.6615995 ≤ 2.6616, first order −0.77027 ≤ −0.7702, threshold e^{7.1476/0.7702} = 10723.4 ✓ (wyconst.log). Its inputs (Robin 1984 Lemmas 5–6, Schoenfeld under RH, [W] Lemma 3.3, Lemma 2.6(c)) are refereed and not re-derived.
- (c) l. 81–82: odd n > N: f_{2,1}(n) = f(n) < e^γ/2 by (b) ✓. Even n = 2^k m (k ≥ 1), n > N: if m > N, f_{2,1}(n) = σ(m)/(n log log n) < f(m) < e^γ/2 ✓; if m ≤ N, σ(m)/m ≤ (σ(N)/N)(m/N)^{ε*} ≤ σ(N)/N (N odd CA at ε*) and log log n > log log N, so f_{2,1}(n) < σ(N)/(2N log log N) = 0.445809 ✓. f_{2,1}(N) = 0.8916175219 > e^γ/2 = 0.8905362091 ✓ (values.log). "Without RH, on the ranges of (a)": even n ≤ X has odd part ≤ X ✓.
- [W] l. 85: the question is whether f_{2,1} exceeds its limsup infinitely often; [W] answers "eventually below" under RH; (c) gives the exact last exception N (under RH) ✓.

## C6. Residual risks ([P] l. 84–89) — confirmed as listed: §4 inherits RobinCH + Chirre–Helfgott (+ half-plane for 46.79); §§1–3 inherit nothing analytic. Float sorting of events: harmless, every near-tie and merge position is arb-certified (C3).

