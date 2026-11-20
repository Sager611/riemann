OPEN: the needed bound is not falsified (the kept Bessel phase costs nothing in exact Kloosterman-fraction sums, which sit at the random-model size N^{−1/2}, far below the needed Q^{−1/4}X0^δ) and is neither equivalent to nor implied by a named open conjecture; the closest known result is the Bombieri–Iwaniec double large sieve, which keeps the phase and in e1rate's own bookkeeping gives ν < min{2θ−1, θ/2} on (3/5, 1), so route (a)'s target is met on [5/8, 0.8] modulo one unproved spacing lemma, while on (3/5, 5/8) the exact gap is a Deshouillers–Iwaniec-type type-II bound with the phase kept losing less than X0^{(5−8θ)/(4(1−θ))} (1/8 at θ = 0.6).

# e1phase — falsification of ceilaudit10's e1rate lead (type II keeps the Bessel phase)
Orca e1phase, 20 Nov 2026, 16:56–18:2x UTC. Inputs read-only: Riemann/Deep/E1Rate/{DERIVATION.md, exponents.py}, Riemann/Deep/CeilAudit10/{REPORT.md §2, code/e1_routes.py, e1_delta.py}.
Exponents: powers of T; relative = divided by the main term ≍ H = T^θ; ε-losses suppressed. Code/logs: Riemann/Deep/E1Phase/code/ (memcap peaks ≤ 155 MB, max RSS ≤ 176 MB by /usr/bin/time -l; memcap lines in logs).

## 1. The needed bound, exactly (e1rate notation)
Defining lines. DERIVATION.md l. 12–13: O = Σ_g (4π/g) Σ_{(p,q)=1} Σ_n β_p γ_q α_n e(nq̄/p) Φ(n/(pq)), Φ(x) = e(−x)𝒲(x), |β_p| ≤ P^{−1}, |γ_q| ≤ Q^{−1},
|α_n| ≤ d(n)n^{O(1/L)}; l. 11: 𝒲(x) = −π∫A(ξ)ξ^{−z/2}B(4π√(ξx))dξ; l. 16: N0 = PQT/H², X0 = 2√(2π)√(NT/(PQ)) (ξ ≍ T/2π); l. 17: 𝒲 = Σ_± e^{±iX}·amplitude
of size H X0^{−1/2}; l. 22–24: Mellin separation costs Λ = √X0; l. 28–29: R_triv = N^{3/4}(PQ)^{1/4}T^{−1/4}, R_sep = R_triv·√X0; l. 38–41: DI Lemma 1 and its
use (v = p modulus, q = ab, ℓ = n, C(a,ℓ) arbitrary, gain G = {…}/(ABNV)^{1/2}); l. 153–154: the escape route ("keep the Bessel phase … hybrid large sieve").
exponents.py l. 4–6 gDI (DI's gain), l. 7–12 bounds (x0 = log_T X0 = (n+1−pv−kq)/2, etr = log R_triv, esep = etr + x0/2), l. 13–14 eII = esep + gDI, l. 17–31 Bdec
(Vaughan cut u, pieces 2, 3); ceilaudit10 e1_routes.py l. 9–10 (eII_kept = etr + gDI), e1_delta.py l. 8–9 (eII_d = etr + δ·max(0,x0) + gDI); REPORT l. 103–104.
(B_δ) For every block (V = P or Q the modulus, A·B = the other variable, N ≤ N0T^ε) and all |C(a,ℓ)| ≤ 1, (ab, v) = 1:
  Σ_{b≍B} Σ_{v≍V} |Σ_{a≍A} Σ_{ℓ≍N} C(a,ℓ) e(ℓ\overline{ab}/v) Φ_N(ℓ/(vab))| ≪ T^ε X0^δ ‖Φ_N‖_∞ (ABNV)^{1/2}{(VB)^{1/2} + (N+A)^{1/4}[VB(N+A)(V+A²) + NA²B²]^{1/4}},
  Φ_N = Φ·(dyadic cutoff at x ≍ N/(PQ)), Bessel phase ±4π√(ξℓ/(vab)) kept (NOT Mellin-separated), δ < 1/2 (δ = 1/2 is DI after separation).
Load-bearing instance (diag.py → diag.log, e1rate's grid, exponents.py imported read-only). At θ = 0.6 (δ = 0 and 1/16) and θ = 0.7 (δ = 0) the binding block is
P = Q = ν, N = N0 with Vaughan cut U = 1 and A = 1 (q = b whole, outside), binding DI term the spectral one, G = Q^{−1/4} (θ = 0.6: −0.0579 vs diagonal −0.1316,
NA²B² −0.1158). At θ = 0.8 the binding piece is type I (Lemma I, unchanged): phase keeping on type II is not load-bearing there. So the input is
(★) Σ_{p≍P} Σ_{q≍Q, (p,q)=1} |Σ_{n≍N} C(n) e(nq̄/p) Φ_N(n/(pq))| ≪ T^ε X0^δ ‖Φ_N‖_∞ PQN·Q^{−1/4}   (P = Q = y = T^ν, N = N0 = y²T^{1−2θ}, X0 ≍ T^{1−θ}),
  for arbitrary |C(n)| ≤ 1: no Vaughan identity, no μ-structure, absolute values over (p, q). Block condition: etr − ν/4 + δx0 < 0 with etr = 2ν + 1/2 − 3θ/2, x0 = 1 − θ,
  i.e. ν < (2/7)(3θ−1) − (4/7)δ(1−θ) (θ = 0.6: 0.22857 − 0.2286δ, = ceilaudit10's e1_delta.log). Sizes at θ = 0.6, ν = 0.22857: P = Q = T^{0.2286}, N0 = T^{0.2571},
  X0 = T^{0.4}; X0/N0 = T^{0.143} (the chirp e^{iX0√(n/N)} sweeps ≫ 1 frequencies); (★) asks a saving T^{0.0571−0.4δ} over E1's trivial bound PQN‖Φ_N‖_∞.
Observation: for A = 1 the plain large sieve for the fractions q̄/p (spacing ≥ 1/(4P²)) already beats DI Lemma 1 phase-free: gain 2(P/(QN))^{1/2} = 2N0^{−1/2}
  = T^{−0.129} vs Q^{−1/4} = T^{−0.057}; after Mellin separation (×√X0) it gives ν < θ − 1/2 at the block (useless); linearizing the chirp on n-intervals of
  length N/√X0 costs X0^{1/4} (heuristic, assuming the shifted points stay 1/P²-spaced) and gives exactly E1's ν < (3θ−1)/4 at the block.

## 2. Falsification attempts
### 2.1 Direct computation of (★) with the phase kept (ls.py → ls.log; exact q̄ = q⁻¹ mod p; peak 155 MB)
Model: p, q ∈ [y, 2y), (p,q) = 1; n ∈ [N, 2N), w(n/N) = sin²(π(n/N−1)); Φ = w·exp(iX0√((n/N)·y²/(pq))) (phase X0 at the centre, as l. 16–17; amplitude dropped);
sizes tied as at θ = 0.6, ν = 0.22857: N = y^{1.125}, X0* = y^{1.75}. Relative l¹ = Σ_{p,q}|Σ_n C(n)Φ e(nq̄/p)| / (R·Σw). C: d(n)/max (structured), random
unimodular (mean of 3), adversarial = best of 6 Boyd ∞→1 ascents (starts: d(n), dechirped d(n), random, conjugated rows). Control "randfr": q̄/p → uniform random.
| y | R | N | X0* | DI(A=1) relative bound | adv at X0 = 0, X0*, 10X0* | adv(X0*)/adv(0) | random C | d(n): 0 → X0* | randfr adv at X0* | Δ/R (0, X0*) |
|---|---|---|---|---|---|---|---|---|---|---|
| 40 | 986 | 63 | 636 | 0.532 | 0.1603, 0.1573, 0.1602 | 0.981 | 0.1378 | 0.0341 → 0.0419 | 0.1586 | 1.10, 1.09 |
| 60 | 2182 | 100 | 1294 | 0.464 | 0.1258, 0.1245, 0.1250 | 0.990 | 0.1081 | 0.0245 → 0.0350 | 0.1257 | 1.14, 1.09 |
| 90 | 4920 | 158 | 2630 | 0.407 | 0.0995, 0.0981, 0.0976 | 0.986 | 0.0876 | 0.0176 → 0.0294 | 0.0984 | 1.16, 1.08 |
| 130 | 10292 | 239 | 5005 | 0.363 | 0.0803, 0.0792, 0.0781 | 0.986 | 0.0699 | 0.0119 → 0.0225 | 0.0792 | 1.17, 1.07 |
- The kept phase costs nothing: adversarial ratio 0.98–0.99 for X0 from 0 to 10·X0* (no X0^{1/2}, no X0^{1/8}); l² large-sieve constant Δ = s_max² ≈ 1.1R in all cases.
- Kloosterman fractions behave like random fractions (randfr column equal within 1%); the adversarial ascent gains only 15% over random C.
- Scaling: adversarial relative size ∝ y^{−0.582} (40 → 130), the random-model N^{−1/2} = y^{−0.5625}; the needed Q^{−1/4}X0^δ = y^{−0.25+1.75δ}: (★) holds with a margin
  growing like y^{0.31} (δ = 0); DI's bound (0.36 at y = 130, constants 1) is 4.6× above the adversarial value. Structured C = d(n): the phase raises it ×1.2–1.9 (the
  mean value of d(n) acquires stationary points of the chirp) but it stays below random C.
- Limits: y ≤ 130 (T = y^{1/ν} ≤ 2·10⁹), local ascent (not a certified ∞→1 norm), amplitude and the ξ-integral of 𝒲 omitted (they only smooth). No counterexample.
### 2.2 Spectral side: the Kuznetsov route pays the separation price (Weyl-law model; spec.py → spec.log)
Cauchy–Schwarz over (p,q), Poisson in q mod p, Kuznetsov in p: the weight carries the phase Z ≍ X0|n−n'|/N in log p at Kuznetsov argument x ≍ (NX0/(PQ))^{1/2}
(= T^{0.1} at the block) ≪ Z (≤ T^{0.4}). There J_{2it}(x) ≈ (x/2)^{2it}/Γ(1+2it), so |h(t)|/cosh πt ≈ |φ̃(2it)||2t|^{−1/2}/√(2π), φ̃ = Mellin transform of the
test function. Computed (φ(u) = w(u)e^{iZ√u}, u ∈ [1,2]): ∫|φ̃(2it)|dt = 0.240√Z (Z = 1000; 0.268√Z at 300), spectral mass at mean |t| = 0.30Z, Weyl-law mass
M(Z) = (1/6)∫|φ̃(2it)||t|^{1/2}dt = 0.0221Z (Z = 1000; 0.0246Z at 300; 0.93 at Z = 0). The spectral large sieve (K² + M) at K ≍ Z thus charges ≍ Z in l², √Z ≤ √X0 in l¹:
δ = 1/2 returns. The Weyl term is sharp (Σ_{t_j≤K}|ρ_j(1)|²/cosh πt_j ≍ K²), so no proof of (★) with δ < 1/2 for arbitrary C(n) can go through Kuznetsov + spectral
large sieve with |h|; it must avoid spectral theory or use the sign changes of h(t_j) ~ (x/2)^{2it_j} (spectral cancellation, structure of α_n = σ_{−α,β}(n)).
### 2.3 Does (★) or (B_δ) imply or contradict a named statement? (checked; none found)
- Selberg's eigenvalue conjecture: the binding block is level 1 (modulus p, e(nq̄/p); SL2(Z) has λ1 = 91.14 > 1/4): neither needed nor implied. DI's level-dependent
  (V + A²) terms are not binding at δ < 1/8 (diag.log).
- Linnik–Selberg: (★) is an average over PQ moduli/twists with saving T^{0.057} ≪ square root; no pointwise Kloosterman-sum bound follows.
- Lindelöf / sub-Weyl: the P = Q = 1 block is the E(T)-type sum Σ d(n)e^{iX0√(n/N0)} (window mean square), but the minimax never uses (★) there (trivial bound
  suffices at small P, Q); amplification with arbitrary b_h gives only |ζ(1/2+it)|² ≪ H T^ε (|A(t0)|² ≍ y, diagonal ≍ Hy): no subconvexity follows.
- Sharpness of the spectral large sieve: not contradicted ((★) is not a spectral statement); it only blocks the Kuznetsov route (2.2).
- Ω-results for Kloosterman sums / E(T) = Ω(T^{1/4}): not contradicted (averaged statement, saving far from square root).

## 3. Closest known result: the Bombieri–Iwaniec double large sieve keeps the phase
Lemma (Bombieri–Iwaniec 1986, double large sieve). x_r, y_s ∈ R^K, |x_r^{(k)}| ≤ X_k, |y_s^{(k)}| ≤ Y_k ⇒ |Σ_rΣ_s a_r b_s e(x_r·y_s)|² ≪_K ∏_k(1 + X_kY_k)·B(a; 1/Y)·B(b; 1/X),
  B(a; η) = Σ_{|x_r^{(k)} − x_{r'}^{(k)}| ≤ η_k ∀k} |a_r a_{r'}| (B(b; η) likewise).
Use. The kept phase of O (every block, arbitrary β_p γ_q α_n, no Vaughan) is bilinear: e(nq̄/p)·e^{±4πi√(ξn/(pq))} = e(x_r·y_n), x_r = (q̄/p, ±2(pq)^{−1/2}) (q̄/p ∈ [0,1)
  as a real number), y_n = (n, √(ξn)). The rest of Φ — e(−x) (x ≤ T^{1−2θ+ε} ≪ 1), the Hankel symbol X^{−1/2}a_±(X), cutoffs — is a non-oscillating function of x = n/(pq)
  on a dyadic range: Mellin separation costs T^ε (not √X0); the ξ-integral goes outside (∫|G| ≪ H, DERIVATION §7 point 4). Half-ranges after translations (absorbed in
  a, b): X1 = 1/2, Y1 = N/2, X2 = (PQ)^{−1/2}/2, Y2 = (√2−1)√(ξN)/2, so (1+X1Y1)(1+X2Y2) ≍ N·max(1, X0).
- B(b) ≤ Σ_{|n−n'|≤2}|b_n b_n'| ≤ 5Σ|b_n|² ≪ NT^ε (from the first coordinate alone; reviewer); n ≠ n' is even excluded once T > 64πNPQ (ν < θ/2).
- B(a) ≤ (PQ)^{−2}·#{close pairs}: |q̄/p − q̄'/p'| ≤ 2/N and |pq − p'q'| ≲ 60PQ/X0. Spacing lemma (SL): #close pairs ≪ T^ε(R + R²/(N·max(1,X0))), R ≍ PQ.
- Then |O_block| ≪ T^ε H X0^{−1/2}(PQ)^{−1} N (X0R + R²/N)^{1/2}: relative R_DLS = R_triv·(X0/(PQ) + 1/N)^{1/2} = N(PQ)^{−1/2} + (N/X0)^{1/2} (X0 → max(1,X0)).
  Both terms increase in N, P, Q: max at P = Q = y, N = N0: yT^{1−2θ} + yT^{−θ/2}. So the J/Y dual branch is o(H) for ν < min{2θ−1, θ/2}; SL's other steps hold at
  ν < 1/2 (e1rate §4). Price check (lesson 329): the method's price is the factor 1 + X2Y2 ≍ X0, kept; the gain is the 2D sparsity of the points (B(a) ≈ R), which
  the α-spacing alone (Farey spacing 1/(4P²)) does not give (then the relative bound is ≥ (X0/N)^{1/2} > 1). θ = 1 check (consistency only: N0 = PQ/T < 1 for
  ν < 1/2, the dual range is empty there): ν < 1/2, the Balasubramanian–Conrey–Heath-Brown range for arbitrary coefficients ✓.
- SL numerics (dls_count.py → dls_count.log): B_a/R for the arithmetic points vs the same pq with uniform random α, block scaling:
| scaling | y = 40 | 80 | 160 | 320 |
|---|---|---|---|---|
| θ = 0.6, ν = 0.22857 (X0 = y^{1.75}, N = y^{1.125}) | 16.64 / 17.19 | 9.73 / 10.15 | 5.85 / 6.04 | 3.68 / 3.75 |
| θ = 0.7, ν = 0.35 (X0 = N = y^{0.857}) | 159.2 / 162.5 | 352.7 / 354.8 | 792.9 / 795.2 | 1547 / 1549 |
  Arithmetic ≤ random in all 8 cases (ratio 0.97–0.999), both = the heuristic 1 + (x2-neighbours)·2/Y1 within 4%; at θ = 0.6 B_a/R → 1 as predicted (R/(NX0) → 0).
  Stress (dls_count2.py → dls_count2.log, θ = 0.6 scaling): y = 640, all p, q: 2.464 / 2.505 (R = 249 276); p, q prime: 1.107 / 1.116 (y = 320), 1.049 / 1.051 (640),
  1.021 / 1.023 (1280, R = 28 056, X0 = 2.7·10⁵): no arithmetic clustering of (q̄/p, pq) at scale (1/N, PQ/X0).
- e1rate's bookkeeping with this one extra bound (book.py → book_dls.log; min(etr, eDLS) per block, grid G = 16, 17 bisections): ν = max{E1, min(2θ−1, θ/2)} exactly.
| θ | 0.55 | 0.6 | 0.62 | 5/8 | 0.65 | 2/3 | 0.7 | 0.75 | 0.8 | 0.85 | 8/9 | 0.9 | 0.95 | 1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DLS (+SL) | 0.1625 | 0.2000 | 0.2400 | 0.2500 | 0.3000 | 0.3333 | 0.3500 | 0.3750 | 0.4000 | 0.4250 | 0.4444 | 0.4500 | 0.4750 | 0.5000 |
| E1 (3θ−1)/4 | 0.1625 | 0.2000 | 0.2150 | 0.2188 | 0.2375 | 0.2500 | 0.2750 | 0.3125 | 0.3500 | 0.3875 | 0.4167 | 0.4250 | 0.4625 | 0.5000 |
| e1rate ν_K | 0.1625 | 0.2000 | 0.2150 | 0.2188 | 0.2375 | 0.2500 | 0.2750 | 0.3125 | 0.3500 | 0.4000 | 0.4444 | 0.4571 | 0.5000 | 0.5000 |
| route (a), δ = 0 (conditional) | 0.1857 | 0.2286 | 0.2457 | 0.2500 | 0.2714 | 0.2857 | 0.3143 | 0.3571 | 0.4000 | — | — | — | — | — |
  Binding block P = Q = ν, N = N0 (θ ≤ 2/3: term yT^{1−2θ}; θ ≥ 0.7: term yT^{−θ/2}). DLS+SL ≥ route (a) on [5/8, 0.8] (equal at both ends), > E1 on (3/5, 1),
  > e1rate's ν_K on (3/5, 8/9); at θ ≤ 3/5 nothing ((2θ−1) ≤ (3θ−1)/4). In route (a)'s currency DLS+SL equals loss δ_DLS(θ) = (5−8θ)/(4(1−θ)) on the A = 1 block
  (θ = 0.6: 1/8 — exactly ceilaudit10's tolerance, so at θ = 0.6 itself DLS gives E1 and nothing more).
  Combined with e1rate's separated DI/Weil blocks (book.py all → book_all.log): 0.42500 (θ = 0.85), 0.44445 (8/9), 0.45714 (0.9) = the larger of the two; no synergy.

## 4. The exact gap and what is missing
- Status of the spacing lemma SL (the only unproved input of §3): open. Necessary-condition form: a close pair forces k·pp' = h·qq' + j with |h| ≤ 8P²/N, |j| ≤ M ≍ PQ/X0
  (from qq̄ = 1 + ps, k = sq' − s'q); the random count of its solutions is ≍ R²/(NX0) (matches the numerics), but proving it needs equidistribution of hqq'/(pp')
  mod 1 at scale M/P² on average over (h, p, p') — per (h, p, p') the Weil/complete-sum error QP exceeds the main term Q²M/P² by P²X0/Q² — i.e. averaged
  Kloosterman-type bounds (DI/BC class, but phase-free: the Bessel phase is gone, only the shell width 1/X0 remains). The shell-wise large sieve (Farey spacing 1/(4P²))
  gives #pairs ≪ R²/(NX0) + RP²/N, which makes the DLS bound useless (relative N^{1/2}): the needed gain is exactly beyond Farey spacing. SL is not implied by, nor
  equivalent to, any conjecture named in the brief (it is a counting statement with saving R/X0 = T^{0.057} over the trivial count at θ = 0.6).
- θ ∈ [5/8, 0.8] (and up to 8/9 against e1rate): route (a)'s target is met or exceeded by DLS + SL; the missing piece is SL alone, not a hybrid Kuznetsov bound.
- θ ∈ (3/5, 5/8): DLS + SL gives 2θ − 1 < (2/7)(3θ−1); beating it with route (a) needs (B_δ) with δ < (5−8θ)/(4(1−θ)) (θ = 0.6: 1/8; 0.62: 0.026; 5/8: 0); beating
  E1 needs δ < (3θ−1)/(16(1−θ)) (θ = 0.6: 1/8). Closest known DI-type results: DI 1982 Thm 12 / 1984 Lemma 1 (smooth weights only: oscillation X0 ⇒ δ = 1/2 via
  separation, §2.2 via Kuznetsov), Bettin–Chandee trilinear (no phase), Gallagher-type hybrid large sieve (continuous t-range: constant N + P²X0, relative (X0/N)^{1/2} > 1,
  useless). Missing: a Kloosterman–archimedean bound exploiting sign changes of h(t_j) at |t_j| ≍ X0, or SL-type sparsity inside DI's proof.
- θ ≤ 3/5: unchanged (DLS ≤ E1; route (a) still needs (B_δ), δ < (3θ−1)/(16(1−θ)), no candidate known).
- Zeros (lb_dls.py → lb_dls.log; SL/E1 detector, kappa_opt deg 6): DLS + SL raises LB = max{0, c, (c+2κ)/3, h} only on θ ∈ (0.6, ≈0.61): 0.143152 → 0.144385
  (θ = 0.602), 0.145503 → 0.147969 (0.604), 0.149093 → 0.151542 (0.606), 0.153868 → 0.155106 (0.608); ceilaudit10's route (a) at θ = 0.6 (0.1408 → 0.1550) stays
  conditional on (B_δ), δ < 1/8 — DLS lands exactly on that threshold there and gives nothing at θ = 0.6.
- Not done (time): certified ∞→1 norms; y > 130 for (★), y > 320 for SL; Arb for κ/LB; reviewer verdict and novelty check on §3 (no discovery claimed).

## 5. Review (own checks + one verdict-only reviewer)
- Own: (★) = e1rate's own binding instance (diag.log reproduces ceilaudit10's ν = 0.22857 − 0.2286δ and the binding corner); book_dls.log reproduces
  ν_DLS = max{E1, min(2θ−1, θ/2)} on the grid at 14 θ (worst block P = Q = ν, N = N0); spec.log reproduces e1rate F1's √X0 law (0.240√Z).
- Reviewer (verdict-only on §3, inputs pasted; timed out at 170 s, recap kept): (1) phase identity and T^ε separation CORRECT; (2) DLS form CORRECT;
  (3) GAP, constant only: diagonality needs T > 64πNPQ, but B(b) ≪ NT^ε follows from |n−n'| ≤ 2 alone (fixed above); (4) bound and ν < min{2θ−1, θ/2}
  CORRECT conditional on SL; (5) θ = 1 check is consistency only (N0 < 1, empty dual range; fixed above); (6) SL is the only substantive unproved input,
  no omitted power (separation, translations, coefficient weights, (p,q) = 1, ξ-integral); full moment reduction and N > N0 tails (E1 Step 4) not re-verified.
  Overall: GAP unconditionally (SL), the conditional deduction works.
- Novelty (light, no claim): arXiv search (approach arxiv) "double large sieve Kloosterman fractions Bessel phase", "twisted second moment zeta short intervals
  mollifier Kloosterman", "double large sieve … short interval zeta": no hit applying the double large sieve to the short-window dual (only Tang 2608.14852,
  2401.01057, θ = 1 or single pairs); our repo: "double large sieve" only in Riemann/Deep/ExploreDeterminant (unrelated). Before any DISCOVERIES entry: prove SL,
  full reviewer, monthly-listing novelty scan (AGENTS.md).

## Suggested NOGO line
- e1phase (Orca, N8 falsification of ceilaudit10's e1rate lead): OPEN, not falsified — its load-bearing instance is (★) Σ_{p,q}|Σ_n C(n)e(nq̄/p)Φ(n/(pq))| with A = 1 (no Vaughan, arbitrary C), which with the Bessel phase kept costs nothing numerically (adversarial ratio 0.98–0.99 up to 10X0*, size N^{−1/2} ≪ needed Q^{−1/4}); Kuznetsov + spectral large sieve with the phase in the test function pays the sharp Weyl term (mass ≍ X0, δ = 1/2 returns), so any proof must be non-spectral or sign-exploiting; the Bombieri–Iwaniec double large sieve keeps the phase bilinearly and gives ν < min{2θ−1, θ/2} (> E1 on (3/5,1), ≥ route (a) on [5/8,0.8]) modulo the open spacing lemma SL for (q̄/p, pq) at scale (1/N, PQ/X0) (numerically true to 1%, beyond Farey spacing); zeros: LB +≤ 0.0025 only on θ ∈ (0.6, 0.6101); route (a) at θ = 0.6 still needs δ < 1/8.

## Files
Riemann/Deep/E1Phase/: REPORT.md; code/diag.py → diag.log (binding block/piece/DI term, route (a), θ = 0.6, 0.7, 0.8); code/ls.py → ls.log ((★) with exact q̄, X0 scan,
adversarial/random/d(n), random-fraction control; peak 155 MB); code/spec.py → spec.log (Weyl-law mass of the Kuznetsov majorant); code/dls_count.py → dls_count.log and
code/dls_count2.py → dls_count2.log (SL counts, all p, q and primes); code/book.py → book_dls.log (bookkeeping with the double large sieve), book_all.log (with e1rate's
DI/Weil blocks, θ = 0.85–0.9); code/lb_dls.py → lb_dls.log (LB effect, crossover 0.6101). .scratch/e1phase/STATE.md. No downloads; nothing written outside these paths.
