RECORD: for N ≥ 4 distinct primitive Dirichlet characters of one modulus and parity and any nonzero real c_j, at least 6.81·10⁻⁴ (N = 4), 5.75·10⁻⁴ (N = 5), …, 3.36·10⁻⁴ (N = 10) of the zeros of Σ_j c_j ε_j q^{s/2−1/4} L(s, χ_j) lie on the critical line, and at least 3.425·10⁻³/N for every N ≥ 11 (9.69·10⁻³/N for N ≥ 10⁶); was 3.70·10⁻⁹ (N = 4) … 1.47·10⁻⁹ (N = 10) and 2.161·10⁻⁶/(N log N), only for log N ≫ 862 (Dousselin 2025)

## Claim
Setting of Dousselin (arXiv 2311.10285 = Forum Math. 37 (2025) 1871–1912, thesis tel-05317244): q ≥ 3, χ_1, …, χ_N distinct even primitive characters mod q (odd: "similarly",
his remark), c_j ∈ ℝ∖{0}, F(s) = Σ_j c_j q^{−1/4} ε_j q^{s/2} L(s, χ_j), κ_F = liminf (N_0(2T,F) − N_0(T,F))/(N(2T,F) − N(T,F)). Then, unconditionally,
  κ_F ≥ 2π(1/A − N(√C1(A) + √C2)²/A³) for every A > 1/ϰ, θ ∈ (0,1), ϰ < 1/8, with C1(A) ≤ Ka·A + Kb **linear in A** and C2 = 4(C3 − 1)/(θϰ)² (DERIVATION §3, §9–§10).
Certified values (ϰ = 0.1249999; counts zeros of odd order, distinct):

| N | κ_F ≥ (new) | θ | A | Dousselin (paper / thesis) | ratio |
|---|---|---|---|---|---|
| 4 | 6.8158·10⁻⁴ | 0.150 | 5265.25 | 3.68·10⁻⁹ / 3.70·10⁻⁹ | 1.8·10⁵ |
| 5 | 5.7578·10⁻⁴ | 0.150 | 6196.71 | 2.94·10⁻⁹ / 2.95·10⁻⁹ | 2.0·10⁵ |
| 6 | 5.0081·10⁻⁴ | 0.150 | 7090.67 | — | |
| 7 | 4.4455·10⁻⁴ | 0.150 | 7955.86 | — | |
| 8 | 4.0060·10⁻⁴ | 0.150 | 8798.05 | — | |
| 9 | 3.6520·10⁻⁴ | 0.150 | 9621.28 | — | |
| 10 | 3.3601·10⁻⁴ | 0.150 | 10428.5 | 1.46·10⁻⁹ / 1.47·10⁻⁹ | 2.3·10⁵ |
| 100 | 5.1386·10⁻⁵ | 0.100 | | 1.45·10⁻¹⁰ / 1.46·10⁻¹⁰ | 3.5·10⁵ |
| 1000 | 6.8655·10⁻⁶ | 0.050 | | 1.43·10⁻¹¹ / 1.44·10⁻¹¹ | 4.8·10⁵ |

Large N: κ_F ≥ c/N for all N ≥ N1 with c = 3.4254·10⁻³ (N1 = 11), 5.1183·10⁻³ (100), 6.8429·10⁻³ (1000), 9.6960·10⁻³ (10⁶); the method's limit is π/(2Ka) ≥ 9.89·10⁻³ (θ = 0.03).
This is the order c/N; Selberg proved c/N² and announced c/(N log N), which Dousselin made explicit. (N = 1, 2, 3 for completeness: 1.90·10⁻³, 1.15·10⁻³, 8.48·10⁻⁴; for N = 2, 3 our
LinComb3 bounds 34.707%, 2.061% are far better, so the record is N ≥ 4.)

## Before
Dousselin, Theorem 'application' (TeX l.619–640): 7.38, 4.91, 3.68, 2.94, 1.46·10⁻⁹ for N = 2, 3, 4, 5, 10, and κ_F ≥ 2.161·10⁻⁶/(N log N)·(1 − 3ε − 1.14 log log N/log N − 862/log N − …) for
N ≥ 2.9·10⁻¹¹/ε³ (positive only for log N > 862); thesis (2025, Thm I.1.2): 7.41, 4.93, 3.70, 2.95, 1.47·10⁻⁹. Reproduced in Arb to all printed digits (cert.log Part 1).
No other explicit bound for N ≥ 4 (NOVELTY.md). Our LinComb3 method gives nothing for N ≥ 4 (needs a non-concentration lemma).

## Method (each step re-derived in DERIVATION.md; inputs from Dousselin's paper are cited by TeX line)
1. Counting (§1): #sign changes ≥ meas(G)/H instead of Dousselin's pairs of intervals: main term 1/A instead of 1/(2A).
2. Combination (§2): his N = 1 argument (Cauchy–Schwarz on the no-sign-change set) run inside each dominance set S_j^*: N(√C1+√C2)² instead of 4N(C1+C2). Steps 1–2 alone: ×12.3.
3. M-estimate (§3): Minkowski over the rectangle's sides and the first moment ∫Lη² = T + o(T) (same rectangle, pointwise bounds at two heights): C2 = 4(C3−1)/(θϰ)²
   instead of 6(C3+1+2√C3)/(θϰ)².
4. Exact main terms (§4): Δ_2(X) = −2√X e^{−iX} − (1+2iX)ψ(X), ψ = Fresnel-type erf; Dousselin's own Lemma 12/eqintero (and its v → 0 limit) give C5' = 2/√(πϰ(1−θ))
   (his Lemma A contour bound is 2.66× larger) and exact profiles S_X(u), F4_X(u) for the I-integral's small-v range (instead of the contour bounds C6, C4).
5. Weights (§5, Lemma Φ): d | δδ' forces u_δ + u_{δ'} ≥ log d/log ξ, so the δ-sums are bounded by Φ(t_d)/φ(d) with Φ(t) = sup_{u+u'≥t} f(u)g(u'), integrated in t (Lemma W).
6. Small v (§7): exact cancellation of the local factors of φ_{iv}(d) against Π_{qd}(v)² (Euler product P2 ≈ 79.9 → P1 ≈ 11.54), |sin(vΛ_d)|/v averaged over d instead of bounded by (1/2+2ϰ)log T.
7. Large v (§8): integration by parts against (√𝒯 d/(δ1δ2))^{iv} (twice): no A log A term; C1 is linear in A — this is what turns c/(N log N) into c/N;
   log ξ/Λ ≤ ϰ/(1/2 + (t_d − 2)ϰ) is kept inside the d-average.
8. Small primes (§10): at every prime p ≤ 10⁶ the signed local δ-sum Σ_{a+b≥k}τ_{−1/2}(p^a)τ_{−1/2}(p^b)p^{−a−b} = −1/p (k = 1), 0 (k ≥ 2) replaces the absolute 1/φ(p^k);
   the Euler factor of each such p becomes exactly 1, so P1 → ∏_{p>10⁶}(1+f1(p)) ≤ 1.000003001. Ablation (N = 4): 4.54·10⁻⁸ (steps 1–2) → 4.15·10⁻⁵ (+ 3–7 with sup weights,
   C2 = 4(√C3+1)²/(θϰ)², constant log ξ/Λ bound) → 9.73·10⁻⁵ (+ Lemma Φ, first moment, r(t)) → 6.82·10⁻⁴ (+ 8).
9. 2609.15329 (optimised Selberg method, ζ only) was read (intro, §4): it uses Zhuravlev's Littlewood–Jensen detector and an exact four-variable Euler-product evaluation, not the I/J/M
   scheme; none of its steps transfers directly to the dominance sets S_j (Caveat). Our §10 is the analogue of its "exact local factors" observation; as a check, our bound for the
   mollified second moment is 1.13/log ξ against its exact 0.68–0.89/log U.

## Certificate
cd Riemann/Deep/SelbergN4 && uv run --no-project --with python-flint --with mpmath --with numpy python3 selbergn4.py   (300 s, one thread, 95 MB; logs/cert.log, logs/result.json)
Arb (python-flint 0.9.0, 128 bits). Part 1 reproduces Dousselin (asserted within 0.6% of every printed value; P1, P2 with rigorous tails over primes ≤ 10⁶). Part 2a–2c: ablation.
Part 2d: special functions as Arb balls (series with tail bound for E, IE on [0,1]; erf closed form elsewhere; self-check against the defining double integral), sups over 128 X-boxes ×
640 u-boxes, Φ on 200 t-cells, all floats rounded upward (×(1+10⁻¹²) products, ×(1+10⁻¹⁰) sums), κ_F = lower end of the ball; free parameters (θ, A) from a float search, only
admissibility matters. Asserts κ_F(4), κ_F(10) > 10⁵ × Dousselin and Kb ≤ 0.

## Checks
- Reproduction: N = 1…5, 10 match to the printed digits (3.682008e-09 for N = 4). Paper inconsistencies found: K2's (log ϰ − ϰ) should be (log ϰ − 1) (his choice is conservative);
  his Σ2' bound (l.1487) drops a factor 2 of his own displayed steps; with it restored his table changes only at N = 1 (5.45 → 5.43·10⁻⁸). Neither affects this result.
- Normalisations of the exact main terms checked numerically (checks/mtcheck.py, q = 1, ξ = 10⁶, θ = 0.1; logs/mtcheck.log): sum/main term = 1.04 (B, and S at X = 0.25, 0.5),
  1.03 (X = 1), consistent with O(1/log ξ) corrections.
- Reviewers (7 verdict-only runs on pasted claims, logs/review.log): counting/combination SOUND; first moment and t-dependent log ξ/Λ SOUND; integration by parts SOUND (uniformity: answered by the certified box sups);
  local-factor identity: a missing (q/φ(q)) in the text, fixed (code unaffected); Lemma W: additive error made explicit; C3/Σ2': rephrased as absolute bounds; main terms: C5' SOUND,
  derivative route made explicit; small primes (§10): decomposition and limits SOUND, two exact-constant statements are upper bounds (as used), Im part handled per δ.
- Every constant of DERIVATION §11 rg-found in logs/cert.log (65/65).

## Novelty
arXiv listings math.NT/CV/CA Oct–Nov 2026, arXiv API (Dousselin; "linear combination" + "critical line"; "Selberg's method"), HAL (paper and thesis), Zenodo, one GXL evidence query:
no explicit bound for N ≥ 3 except Dousselin's; no c/N anywhere (NOVELTY.md, logs/novelty.log). 2609.15329 (ζ only) is the closest modern optimisation of Selberg's method.

## Caveat
- Rests on Dousselin's analytic inputs as proved in his paper (not re-proved): dominance decomposition (Tsang-type CLT for differences, moment lemmas), Bauer's twisted second moment,
  the reduction of ∫|I|² to K(v) (eqIntfacile), Lemma 12 and eqintero with their error terms, the Gabriel-convexity argument for M. Our new steps (DERIVATION §1–§10) were checked by
  verdict-only reviewers on pasted statements, not by a specialist; §8 (integration by parts) and §10 (small primes) carry the largest gains and deserve a specialist's reading.
- κ_F is a liminf over T for fixed q, χ_j, c_j; no uniformity in q. ϰ = 0.1249999 < 1/8 as required. Odd characters: by Dousselin's remark ("handled similarly", l.640),
  not re-derived here; the certified computation is for even characters (the constants do not depend on the parity).
- The factor N is structural: the per-j bad sets are bounded by global moments. Removing it needs the joint (non-concentration) lemma of LinComb3; a Zhuravlev–Jensen detector
  (2609.15329) inside the S_j, or an exact evaluation of K(v) (no absolute values in the δ-sums at all primes), are the next steps; neither is done here.

RH relevance: N12, N8 (critical-line proportions for functions with a functional equation and no Euler product — the Davenport–Heilbronn family; first explicit bound of the right
order c/N in the number of summands).

## Files
- Riemann/Deep/SelbergN4/REPORT.md (this), DERIVATION.md (every step and constant), NOVELTY.md, selbergn4.py (certificate; \`--search\` = parameter search), checks/mtcheck.py
- logs/cert.log, logs/result.json, logs/review.log, logs/novelty.log, logs/mtcheck.log
