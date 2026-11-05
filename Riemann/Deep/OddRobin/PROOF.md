# PROOF: the odd Robin inequality on the whole range above N (closing the gap of arXiv 2008.04787)

Notation. For odd n ≥ 3 the *odd Robin inequality* is

  (OR)  σ(n)/n < (e^γ/2) log log n.

N := 3^4·5^3·7^2·11·13···61 = 277093800961598968947975375 and c0 := 67N = 18565284664427130919514350125 (Washington–Yang, arXiv 2008.04787 v1, Prop. 3.5; Int. J. Number Theory 17 (2021) 843–870). They prove (OR) for odd n ≥ c0 under RH and write (TeX l. 400): "We do not know whether the inequality is satisfied for N<n<c0 … the numbers here are large enough to make brute-force computation infeasible."

Odd colossally abundant (odd CA) numbers (WY §2, Definition `MainDef` and Lemma `AlphaLemma` in their TeX): for ε > 0, N_ε = ∏_{p odd} p^{α_p(ε)}, where α_p(ε) maximises g_p(a) := σ(p^a)/p^{a(1+ε)}. As ε decreases, the exponent of p rises from a−1 to a exactly at ε = F(p,a) := log(1 + 1/(p + p² + … + p^a))/log p (since σ(p^a) − p·σ(p^{a−1}) = 1). x_1(ε) is defined by F(x_1,1) = ε. Convention: α_p(ε) is WY's floor formula, which takes the larger maximising exponent at a tie ε = F(p,a); then p | N_ε iff p ≤ x_1(ε). Every odd CA number equals N_ε for all ε in an open interval, so covering every ε with x_1(ε) in a range covers both maximisers at each tie.

All numerics below are arb ball arithmetic (python-flint 0.9.0) unless said otherwise; every inequality is decided with certified signs.

## 1. The gap N < n < c0 (gap.py)

Fix ε* := log(68/67)/log 67 = 0.0035234646… .

(1a) N maximises every g_p at ε*: for each odd p ≤ 997, g_p(α_p+1) < g_p(α_p) and g_p(α_p−1) < g_p(α_p) (α_p = exponent of p in N), except g_67(1) = g_67(0) (a tie: N and 67N = c0 are both odd CA at ε*). For p ≥ 73, α_p = 0 and log g_p(1) − log g_p(0) = log(p+1) − (1+ε*) log p < 0.
   Concavity: log g_p(a+1) − log g_p(a) = log(p + 1/σ(p^a)) − (1+ε) log p decreases in a, so the costs c_p(a) := log g_p(α_p) − log g_p(a) are ≥ 0 and grow as a moves away from α_p.

(1b) Identity. For odd n = ∏ p^{a_p}, put t = log n, t_N = log N, L(n) := Σ_p c_p(a_p) ≥ 0. Then
   log(σ(n)/n) = log(σ(N)/N) + ε*(t − t_N) − L(n).
   (Both sides equal Σ_p log g_p(a_p) + ε* t, using log(σ(N)/N) − ε* t_N = Σ_p log g_p(α_p).)

(1c) Hence n violates (OR) (σ(n)/n ≥ (e^γ/2) log log n) iff L(n) ≤ D(t), where
   D(t) := δ0 + ε*(t − t_N) − (log log t − log log t_N),   δ0 := log(σ(N)/N) − log(e^γ/2) − log log log N = 0.00121349015854… > 0
   (δ0 > 0: N violates; log-margin of c0 is −9.2356e−5 < 0: c0 satisfies).

(1d) D is decreasing on [t_N, t_c0]: by the mean value theorem log log t − log log t_N = (t − t_N)/(ξ log ξ) ≥ (t − t_N)/(t_c0 log t_c0), and ε* = 0.0035234646 < 1/(t_c0 log t_c0) = 0.0036790869. So a violator in the gap has L(n) < D(t_N) = δ0.

(1e) Primes p ≥ 73 cost too much: c_p(a) ≥ c_p(1) = (1+ε*) log p − log(p+1) for a ≥ 1 (concavity, α_p = 0); d/dp of the right side is ((1+ε*)(p+1) − p)/(p(p+1)) > 0, so c_p(1) ≥ c_73(1) = 0.00151163 > δ0. Hence a_p = 0 for p ≥ 73.

(1f) Enumeration for p ≤ 71 (all exponents 0..39; costs beyond grow): exponents with c_p(a) ≤ δ0 are
   3: 4 (0), 5 (0.00111989);  5: 3 (0), 2 (0.00076009);  11: 1 (0), 2 (0.00090169);  67: 0, 1 (both 0);  71: 0 (0), 1 (0.00103316);
   every other p ≤ 71: only a_p = α_p. Combinations with total cost ≤ δ0: 10 vectors. Their n relative to N:
   N/5 (< N), N, 67N = c0, 71N, 67·71N, 3·67N, 11·67N (> c0), and three in the gap:
   | n | n ≈ | log-margin log(σ(n)/n) − log((e^γ/2) log log n) |
   |---|---|---|
   | 3N | 8.31·10^26 | −3.7812·10^−4 |
   | 67N/5 | 3.71·10^27 | −5.0950·10^−4 |
   | 11N | 3.05·10^27 | −5.9604·10^−4 |
   All three satisfy (OR). **So no odd n with N < n < c0 violates (OR).** ∎

## 2. Interpolation (WY Lemma 3.4 = Robin 1984, Prop. 1)

If N1 < N2 are consecutive odd CA numbers, there is ε with both maximal for σ(m)/m^{1+ε} over odd m; so for odd N1 ≤ n ≤ N2, log(σ(n)/n) ≤ log(σ(N1)/N1) + ε(log n − log N1), a linear function of t = log n through both endpoint values, while the threshold log(e^γ/2) + log log t is concave in t. If both endpoints satisfy (OR), so does every odd n in between.

## 3. Direct check of odd CA numbers, c0 ≤ N_ε, x_1(ε) ≤ 5·10^8 (oddca.py)

Events (p,a) are sorted by F(p,a) decreasing; the (p,1) events come in increasing p (F(p,1) decreases in p); the 3,645 events with a ≥ 2 and F(p,a) ≥ F(5·10^8, 1) are merged in at positions certified in arb (F(q,1) > F(p,a) > F(q',1) for the neighbouring primes q < q'); near-ties among them (relative gap < 10^−9) are certified in arb. 26,355,866 odd primes (segmented sieve).
S := log(σ(N_ε)/N_ε) and T := log N_ε are accumulated in arb (128 bits). Up to p = 10^5 every number is checked one by one (13,236 numbers incl. all a ≥ 2 events); the last violator is the number completed by event (61,1), i.e. N (as WY found), and every number from event (67,1) (= c0) on satisfies (OR). Above 10^5, a run of consecutive (p,1) events p_i..p_j is checked as a block: S and T increase along the run, so S_after − log(e^γ/2) − log log T_before < 0 proves (OR) for all numbers in the block; S_after − S_before = log ∏(p+1) − log ∏ p and T_after − T_before = log ∏ p with exact integer products. Failing blocks are halved (down to single numbers, then checked exactly). Result: PASS, 26,346,275 numbers in 79,866 blocks; log N at the last number = 500016054.5165. Closest normalised margin (m·√P·log P, P the largest prime) among exact checks with P > 100: −0.1125 (P = 113); at c0 itself −0.0032.

## 4. Transfer from Robin's inequality at CA numbers, 5·10^8 ≤ x_1 ≤ 1.43·10^47 (transfer.py)

Lemma O1. Let M_ε = ∏_p p^{α_p(ε)} be the colossally abundant number for ε (same α_p for odd p) and a = α_2(ε), so N_ε = M_ε/2^a. Since σ(2^a)/2^a = 2 − 2^{−a},
   (OR) at N_ε  ⟺  R(M_ε) := γ + log log log M_ε − log(σ(M_ε)/M_ε)  >  C_odd := [log log log M_ε − log log log N_ε] − log(1 − 2^{−a−1}).

Lemma O2 (x := x_1(ε) ≥ 41).
 (i) a ≥ ⌊log_2 x⌋, so 2^{−a} < 2/x. [WY lemma `xkLemma2`(a) = Robin 1984 p. 190: x_k > x^{1/k} (from F(t,k) > F(t^k,1): log(1+1/(t+…+t^k)) ≥ log(1+1/(k t^k)) ≥ log(1+t^{−k})/k by concavity); α_2 ≥ k iff x_k ≥ 2 (WY lemma `xkLemma1`: α_p = k for x_{k+1} < p ≤ x_k).]
 (ii) a log 2 ≤ log(((x+1) log x/log 2 + 2)/2). [α_2 = a needs F(2,a) ≥ ε, F(2,a) < 1/((2^{a+1} − 2) log 2), and ε = log(1+1/x)/log x > 1/((x+1) log x).]
 (iii) log N_ε ≥ θ(x) − log 2 ≥ x(1 − 1/log x) − log 2 (with the convention above all odd p ≤ x divide N_ε; Rosser–Schoenfeld θ(x) > x(1 − 1/log x), x ≥ 41). For a tie ε the smaller maximiser is N_ε' for ε' slightly above ε, still inside the boxes (or inside §3 at the left end).
 (iv) −log(1 − 2^{−a−1}) ≤ 2^{−a}, and log log log M − log log log N ≤ a log 2/(log N · log log N) (mean value theorem).
 So C_odd ≤ 2/x + A(x)/(L_N log L_N) with A, L_N the bounds in (ii), (iii).

RobinCH's Level 3 certificate (Riemann/Deep/RobinCH, certify3.margins4, "R") gives for each box of log x a ball m with (√x log x)·R(M_ε) ≥ m for every ε with x_1(ε) in the box (its positive part is normalised at the box's left end, the negative parts are balls over the box). transfer.py reruns it and requires m > (√x log x)·[2/x + A/(L_N log L_N)] (ball upper bound over the box). Valid range of its lemmas: log x ≥ 20 (Lemma 13), x ≥ 10^7 (Lemma 11_k), x_2 ≥ 1423 (Lemma 9); start 5·10^8 (log = 20.03). Result:
 | zero input above 3·10^12 | boxes | min margin after the odd correction | x_1 ≤ X_R | all odd n ≤ |
 |---|---|---|---|---|
 | none (RH to 3·10^12 + Chirre–Helfgott + Büthe only) | 1221 | 1.07·10^−5 | 6.46·10^46 | 10^(10^46.44) |
 | + half-plane Re ρ ≤ 7/8 above 3·10^12 (conditional on it) | 1235 | 1.07·10^−5 | 1.42·10^47 | 10^(10^46.79) |
 | + half-plane Re ρ ≤ 0.87338196661320752 (conditional on it) | 1209 | 1.07·10^−5 | 1.43·10^47 | 10^(10^46.79) |
 The odd correction is ≤ 0.0029 (normalised) at x = 5·10^8 against Robin margins ≥ 0.49 there, and negligible above (≈ 2 log x/√x).
 "all odd n ≤": the last odd CA number checked, N*, has log N* ≥ θ(X_R) − log 2 ≥ X_R(1 − η) − log 2 (η: RobinCH Lemma 12 at X_R).

## 5. Theorem

(a) With no analytic input (§§1–3: arithmetic and arb only): every odd n with N < n ≤ e^{500016054} ≈ 10^(2.17·10^8) satisfies (OR). With RobinCH's inputs (RH to height 3·10^12 (Platt–Trudgian), Chirre–Helfgott's |ψ(x) − x| ≤ πx/(3·10^12) + 113.67√x (preprint arXiv 2512.15709), Büthe's θ bounds below 10^19) and no zero-free half-plane: every odd n with N < n ≤ 10^(10^46.44). Conditional on the half-plane Re ρ ≤ 7/8 (OpenAI 2026, Lean port) or Re ρ ≤ 0.87338196661320752 (Lean-checked per QuasiRiemannTracker, axioms unchecked) for zeros above height 3·10^12: every odd n with N < n ≤ 10^(10^46.79).
Proof: N < n < c0 by §1; c0 ≤ n ≤ N* by §2 applied to consecutive odd CA numbers, all of which satisfy (OR) by §3 (x_1 ≤ 5·10^8) and §4 (5·10^8 ≤ x_1 ≤ X_R). ∎

(b) RH ⟺ (OR) holds for every odd n > N = 277093800961598968947975375. The bound is sharp: N violates (OR).
Proof: (⇒) WY Prop. 3.5 gives (OR) at every odd CA number with x_1 ≥ 120409 under RH; with §3, §2 and §1 every odd n > N satisfies (OR). (⇐) WY Theorem 3.6: if (OR) holds for all large odd n, Robin's inequality holds for all large CA numbers (σ(2^k m)/(2^k m log log 2^k m) < 2·(e^γ/2)), so RH by Robin. ∎

(c) Corollary (WY's question, TeX l. 85, on f_{2,1}(n) = σ(odd part of n)/(n log log n), limsup e^γ/2): under RH, f_{2,1}(n) < e^γ/2 for every integer n > N, and N is the largest n with f_{2,1}(n) ≥ e^γ/2 (f_{2,1}(N) = 0.89162… > e^γ/2 = 0.89054…); without RH, on the ranges of (a).
Proof: odd n: (a)/(b). Even n = 2^k m, k ≥ 1, m odd: if m > N, f_{2,1}(n) < f_{2,1}(m) < e^γ/2; if m ≤ N, σ(m)/m ≤ σ(N)/N (N is odd CA, hence odd superabundant), so f_{2,1}(n) ≤ σ(N)/(2N log log N) = 0.4458… . ∎

## 6. Residual risks

- §4 inherits RobinCH's certificate (reviewed, rerun ALL PASS) and the unrefereed Chirre–Helfgott bound, and its 10^(10^46.79) rows the half-plane hypothesis; §§1–3 use none of it (pure arithmetic + arb).
- Reviewer (verdict-only, 2026-11-05): §§1–2 sound, Lemma O1–O2 and the transfer sign correct, §5(b) and the even-n reduction correct; flagged "unconditionally" for 10^(10^46.79) (half-plane is a hypothesis) and the tie convention — both fixed above.
- Part (b) uses WY's Prop. 3.5 (RH ⟹ (OR) for x_1 ≥ 120409) only for x_1 > 1.43·10^47, where its final bound exp(−0.7702/(√x log x) + 7.1476/(√x log² x)) has a large margin.
- Ordering of odd CA events relies on float sorting plus arb certification of every near-tie (relative gap < 10^−9) and of every merge position; float error in F(p,a) (log1p form, no cancellation) is ~10^−15 relative.
