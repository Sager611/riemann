# AlpHost — PROBES.md (alphost, 20 Nov 2026)

All runs under bash .scratch/memcap.sh 290 <s> (peaks in logs/*.log), single-threaded, uv-only Python. Code in code/, outputs in logs/.

## (a) Mechanism at small scale: quadratic twists ξ·χ_d, χ_d = (·/d), d prime (code/probe_a.py)
Setup: prime powers n ≤ Z = 10^6; ψ(y, χ) = Σ_{n ≤ y} Λ(n)χ(n) on 241 log-spaced y ∈ [10², 10^5]; S(Z, y) = Σ_{n ≤ Z} Λ(n)χ(n) min(1, y/n) (the archive's QH.S, truncated at Z, y ≤ 10^4); L'/L(1, χ) ≈ −Σ_{n ≤ Z} Λ(n)χ(n)(1 − n/Z)/n (Riesz mean; error ≈ Z^{−1/2} log d under GRH; NOT certified). ξ = 1 (M = 1) and ξ = χ_{−4} (M = 4).

### (a1) WinPair-like class: d ≡ 7 (16), d ≡ −1 (105) (z = 7 mimic), all 194 such primes d ∈ [839, 998759] (logs/probe_a_win.log, peak 21 MB)
| M | max_y\|ψ\|/√(y log y) q0/q.1/q.5/q.9/q1 | max_y\|S\|/√(y log y) (y ≤ 10^4) | L'/L(1) q0/q.5/q1 | share \|ψ\| ≤ y^{1/2+1/104} ∀ grid y | share \|S\| ≤ y^{1/2+1/104}/2 |
|---|---|---|---|---|---|
| 1 | 1.09 / 1.46 / 1.91 / 2.45 / 3.05 | 0.35 / 0.70 / 1.47 / 2.17 / 2.53 | −2.65 / −2.05 / −0.86 | 0/194 | 0/194 |
| 4 | 1.15 / 1.55 / 2.02 / 2.52 / 3.65 | 0.38 / 1.16 / 1.97 / 2.69 / 3.81 | −1.27 / −0.04 / +1.51 | 0/194 | 0/194 |
- Mean over d of |S(10^6, y)|²/(y log y) at y = 10², …, 10^4: M = 1: 1.47, 1.20, 1.15, 1.16, 1.09, 0.87, 0.64, 0.56, 0.59; M = 4: 2.48 … 1.49. Bounded by an absolute constant, as QH.meanSqK (MeanSq C) asserts; the small-y excess is the coherent part from χ_d(p) = +1 for p ≤ 7 (the 'RR' term of meanSq_abs).
- L'/L(1, χ_d) ≈ −2.05 (median, M = 1): the forced χ_d(p) = +1 for p ∈ {2, 3, 5, 7} contributes −Σ_{p ≤ 7} log p/(p − 1) = −1.97 (HOSTS §2 structural claim confirmed); |L'/L(1)| ≤ 0.58·log 10² ≪ 4 log y.
- Every d fails the GRH-size test at this scale: the margin y^{ε/2} = y^{1/104} ∈ [1.045, 1.117] on [10², 10^5] is far below the typical normalized size ≈ 2√(log y) ≈ 4–7. This is NOT a counterexample: PairEx is ∃ Y1 ∀ Y ≥ Y1 (see (d) for the scale).

### (a2) Random sample: every 10th prime d ∈ (1000, 10^6], 7833 d (logs/probe_a_sample.log, peak 41 MB)
| M | max_y\|ψ\|/√(y log y) q0/q.1/q.5/q.9/q1 | max_y\|S\|/√(y log y) | L'/L(1) q0/q.1/q.5/q.9/q1 | GRH-size share |
|---|---|---|---|---|
| 1 | 0.86 / 1.57 / 1.99 / 2.51 / 3.57 | 0.21 / 1.06 / 1.79 / 2.47 / 3.89 | −2.58 / −1.69 / −0.61 / +0.68 / +2.49 | 0/7833 (both tests) |
| 4 | 0.87 / 1.54 / 2.02 / 2.59 / 3.76 | 0.28 / 1.16 / 1.90 / 2.68 / 4.21 | −2.47 / −1.41 / −0.37 / +0.77 / +2.57 | 0/7833 |
- |L'/L(1, ξχ_d)| ≤ 0.56·log 10² for every d: the bound 4 log y of PairGood is loose at this scale (kappa4's 4 log Y has room; it needs only the box, no averaging).
- Mean of |S(10^6, y)|²/(y log y): M = 1: 1.93 → 1.40; M = 4: 2.20 → 1.61 (y = 10² → 10^4). Here Z = 10^6 > d, i.e. outside the archive's regime Z² < q (meanSq_abs hypothesis), so off-diagonal terms do not cancel and the ratio is not the diagonal; see (a3).

### (a3) The archive's regime Z² < d: Z = 1000, all 70435 primes d ∈ (10^6, 2·10^6] (code/probe_a2.py, logs/probe_a2.log, peak 27 MB)
Exact square-pairing diagonal Diag(y) = Σ_{nm = □} c_n c_m (c_n = Λ(n) min(1, y/n), prime powers n ≤ Z: same prime with a + b even, or two even powers), the main term of QH.meanSq_abs (MeanSqAbs.lean:8-14: Σ_p sP² + RR² + cross terms).
| y | 10 | 20 | 50 | 100 | 200 | 500 | 1000 |
|---|---|---|---|---|---|---|---|
| M = 1: mean\|S\|²/(y log y) | 2.551 | 2.453 | 2.310 | 2.146 | 1.921 | 1.478 | 0.913 |
| M = 1: Diag/(y log y) | 2.610 | 2.525 | 2.397 | 2.245 | 2.030 | 1.574 | 0.977 |
| M = 1: empirical/Diag | 0.978 | 0.972 | 0.964 | 0.956 | 0.946 | 0.939 | 0.935 |
| M = 4: empirical/Diag | 1.001 | 1.002 | 1.004 | 1.006 | 1.007 | 1.005 | 1.011 |
- With Z² < d the family mean square equals the diagonal to 1–7 % (the M = 1 deficit: (nm/d) for nm ≈ 10^6 ≍ d is not yet equidistributed over d ≤ 2·10^6; the archive has nm ≤ Z² ≤ (log q)^{O(1)}, the Siegel–Walfisz range, where it is). This is the mechanism of QH.meanSqK: orthogonality of (·/d) over primes d via reciprocity and PNT in APs, no large sieve. max_y|S|/√(y log y) quantiles 0/.5/.9/.99/1: M = 1: 0.10/1.74/2.91/3.57/4.66; M = 4: 0.10/1.59/2.85/3.64/5.03.

### (a4) From which Y can PairGood's bounds hold for a positive share? (see (d) for formulas)
- Pointwise: |ψ| ≤ y^{1/2+ε/2} needs y^{1/104} ≥ r√(log y) with r ≈ 3 (observed normalized maxima): log y ≳ 430 (y ≳ 10^{187}); the archive's box_of_S test |S| ≤ y^{1/2+τ}/2 needs log y ≳ 500.
- Markov step of QH.bad_share ((32C/τ)Y^{−τ} ≤ η/(2R), τ = ε/2 = 1/104, C ≈ 2.5 from (a3), η = 1/2, R = 2φ(8M)): Y1 ≳ 10^{540} for M = 1 (table in (d)). Small-scale failure (0 good d out of 7833) is therefore expected and is NOT a counterexample: PairEx is ∃ Y1 ∀ Y ≥ Y1, and its Y1 is also ineffective (Siegel's theorem enters through ExcBound/PNTAPRelW).


## (b) Zero-free boxes: certified argument principle (code/probe_b.py; python-flint 0.9.0 arb balls, prec 64)
Method: winding number of L(s, χ_D) along the boundary of [σ0, 3] × [−T, T]; each boundary segment's ball enclosure must exclude 0 with angular width < π/2 (else bisect), so each increment is certified in (−π/2, π/2) and the sum is 2π·#zeros. Positive control (logs/probe_b_control.log), σ0 = 0.4, T = 10: D = −4: 2 zeros (γ ≈ ±6.02), D = −3: 2 (γ ≈ ±8.04), D = 5: 4.
Main run, σ0 = 1 − 1/400, T = min(8|D|, 200) (the j = 0 height 8M·e^z mimicked by 8·conductor) (logs/probe_b.log, probe_b2.log; peaks ≤ 15 MB):
| D (conductor) | −3 | −4 | 5 | −7 | ±8 | 12 | 13 | −15 | 17 | −19 | −20 | 21 | −23 | ±24 | −163 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T | 24 | 32 | 40 | 56 | 64 | 96 | 104 | 120 | 136 | 152 | 160 | 168 | 184 | 192 | 40 |
| zeros in Re s > 1 − 1/400, \|Im s\| ≤ T | 0 | 0 | 0 | 0 | 0, 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0, 0 | 0 |
- ZeroBox's form is satisfiable and holds for every small real character tested (as the GRH verifications of small conductors predict). D = −163 at T = 200 was stopped after 9 min (rerun at T = 40: 0 zeros, 547 s); D = 1009 hit the 600 s cap (not certified). Not a test of the family statement: ZeroBox for the archive's characters (conductor ≈ e^{Y^{ε/5}}) is only claimed for the chosen pair.
- j = 1 height vs the classical density count (code/probe_cd.py (c2), logs/probe_cd.log): with Q = T = height(M, ε, Y, 1) ≈ 8M e^{z + z log²Y}, Montgomery's (Q²T)^{2(1−σ)/σ}(log QT)^{14} at σ = 1 − 1/400 is e^{0.0150 log Q} (Y = 10^{600}: log Q = 3.88·10^8, log count = 5.83·10^6), the archive's C Q^{1/4} is e^{0.25 log Q}, the d-family has ≥ e^{0.9 z'} = e^{0.9 log Q (1 − o(1))} members: bad share ≤ e^{−0.65 log Q}. The height exceeding d costs nothing.

## (c) ThinH exponents (code/probe_cd.py (c))
- Classical saving at the left end y = e^{x^{ε/3}}/3, λ = 2x^{ε/5}log²x: c log y/(λ + √log y) ≈ (c/2)x^{2ε/15}/log²x; required x^{ε/8} (archive: c = 1); margin exponent 2ε/15 − ε/8 = ε/120. With c_dlVP ∈ {0.1, 0.02} the inequality (with the (log u + λ)² prefactor) holds from log x ≈ 1.7–1.8·10^5 (ε = 1/52), 3.6–3.8·10^5 (ε = 1/104), 2.0–2.1·10^6 (ε = 1/520): x0 astronomically large but finite; c = 1 works for every c_dlVP > 0 (only x0 changes).
- Exceptional zeros: with a Siegel zero, (1 − β) log y ≥ c(δ)N^{−δ}x^{ε/3} = c(δ)·exp(−2δx^{ε/5}log²x + (ε/3)log x) → 0: no saving at all. So RealBox is NECESSARY for real characters with q | cond or d | cond, and it is supplied by PairGood.boxes (ZeroBox contains real zeros; UpperTheory.Top.realBox, AA/CharacterBoxes/RealBox.lean:99-112), i.e. by the family average. Real host characters with neither q nor d in the conductor have conductor | 8·odd(QS) ≤ 8B (QH.real_cond) and are covered by QH.small_zeros: c is uniform in B, x0 = x0(B, ε).
- (e) Principal character (logs/probe_cd.log): G(y) = Σ_{n ≤ y} Λ(n)(y/n − 1) − y(log y − 1 − γ): G(10³, …, 10^7) = 1.53, 1.75, 3.67, 11.0, 27.9, G/√y = 0.048 → 0.009 (≤ Σ_ρ 1/|ρ|² = 0.046 as the explicit formula predicts under RH): the main-term-subtracted definition of gChi(1, ·) is the right object (reviewer round-1 item 7).

## (d) Uniformity of Y1 in M, ε, A (code/probe_cd.py (d))
- Markov step alone ((32C/τ)Y^{−τ} ≤ η/(2R), τ = ε/2, η = 1/2, R = 2φ(8M), C ∈ {1.5, 2.5} from (a3)): ε = 1/52: Y1 ≥ 10^{541}–10^{564} (M = 1), 10^{604}–10^{627} (M = 4), 10^{739}–10^{762} (M = 100); ε = 1/104: 10^{1145}–10^{1587}. log Y1 ≍ (2/ε)log(φ(M)/ε); A enters through Y0(ε, K) of meanSqK/twistK (Siegel-ineffective). Y1 is NOT uniform in M, ε, A, and need not be: winPairs uses one M = ψ.N and one A = J per (ψ, ε).

