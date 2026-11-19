# FULLPROOF: the chain behind t_Z([0,1]) ≥ 0.4222354990

Sources: TeX = arXiv 2610.00386v3 `integer-chebyshev.tex` (v3, 12 Oct 2026; newest); V = `Riemann/Deep/IntCheb/verify.py`;
L = `lowcert.c`; G = `gcert.c`; J = `CERTIFICATE_LOWER_new.json`. Line numbers refer to those files.

## 1. The theorem (TeX 259-317, Theorem thm:cert)
Hypotheses (TeX 259-264): μ a probability measure on [0,1]; U_μ(z) = ∫ log|z−t| dμ(t); Q_1..Q_J ∈ Z[x] pairwise distinct,
irreducible, degrees d_j, leading coefficients a_j ≥ 1; ν_j = NORMALISED root-counting measure; L_j = log(a_j)/d_j;
c_j ≥ 0; B = Σ c_j d_j ≤ 1. A = U_μ − Σ c_j log|Q_j|; G_i = L_i(1−B+c_i d_i) + ∫U_μ dν_i − Σ_{j≠i} c_j ∫log|Q_j| dν_i.
Claim (TeX 272-277): t_Z ≥ exp(min{inf_[0,1] A, −Σ c_j log a_j, min_i G_i}).

Proof, step by step (TeX 279-317), each re-derived by hand:
1. log‖P‖ ≥ ∫log|P| dμ = log|lc P| + n∫U_μ dν (TeX 280-285). Correct: |P| ≤ ‖P‖ on supp μ ⊂ [0,1]; Fubini.
2. Split P over distinct irreducible factors (TeX 285-289). Correct if the factors are taken primitive with lc > 0
   (Gauss) and the integer content c ≥ 1 is dropped: log|lc P| ≥ Σ e_i log a^{(i)}, inequality in the right direction.
   The TeX leaves the content and sign normalisation implicit; no gap (content only increases ‖P‖).
3. (eq:res, TeX 291-296) |Res(Q^{(i)},Q_j)| ≥ 1 for Q_j ≠ Q^{(i)}: distinct primitive irreducibles with lc > 0 are
   non-associate, hence coprime, so Res ∈ Z∖{0}. log|Res| = d_j log a^{(i)} + d^{(i)}∫log|Q_j| dν^{(i)} → eq:res. Correct.
4. Case Q^{(i)} ∉ library (TeX 304-306): H_i ≥ (1−B)L^{(i)} + ∫A dν^{(i)} ≥ inf_C A (needs B ≤ 1, L^{(i)} ≥ 0). Correct.
5. Case Q^{(i)} = Q_i (TeX 308-311): omit j = i; the bound is exactly G_i. Correct (equality, no further step).
6. inf_C A (TeX 313-316): A is superharmonic on Ĉ∖[0,1], +∞ at zeros of Q_j with c_j > 0, continuous into (−∞,+∞]
   (bounded density ⇒ U_μ continuous). Minimum principle ⇒ inf_C A = min(inf_[0,1] A, liminf_{z→∞} A).
   At ∞: A = (1−B)log|z| − Σ c_j log a_j + o(1). So A(∞) = −Σ c_j log a_j ONLY IF B = 1; if B < 1, A(∞) = +∞.
   The TeX sentence "A(∞) = −Σ c_j log a_j" is literally true only for B = 1; for B < 1 the middle term is a
   harmless superfluous weakening. Theorem correct as stated; also correct WITHOUT the middle term when B < 1.
7. Degree bookkeeping: Σ m_i = 1, so ∫log|P|dμ ≥ n·min(inf_C A, min_i G_i); ‖P‖^{1/n} ≥ e^{min}. Constants P: trivial. Correct.
Remark: for a library factor H_i = L_i + ∫U_μ dν_i = g_i := (1/d_i)∫log|Q_i| dμ exactly, and G_i = g_i − (1/d_i)Σ_{j≠i} c_j log|Res(Q_i,Q_j)|
≤ g_i, so the theorem also holds with g_i in place of G_i (the subtraction in TeX 297-301 is a needless weakening for
library factors). No gain here: for the binding x, x−1, 2x−1 every |Res| with the 17 active Q_j equals 1 (resone.py), so G_i = g_i.

## 2. Resultant form of G_i (G 2-7, NOTE.md) and the mass normalisation M ≠ 1 (V 3-7, 120-121)
Identity (checked by expansion): with Res(Q_i,Q_j) = a_i^{d_j} Π_{Q_i(α)=0} Q_j(α) and ∫log|Q_i|dμ = μ(C) log a_i + d_i∫U_μ dν_i,
  G_i^res := (1/d_i)[∫log|Q_i| dμ − Σ_{j≠i} c_j log|Res(Q_i,Q_j)|] = L_i(M − B + c_i d_i) + ∫U_μ dν_i − Σ_{j≠i} c_j∫log|Q_j| dν_i,
M = μ(C). For M = 1 this is the TeX G_i. For M ≠ 1 put μ' = μ/M, c' = c/M: μ' is a probability measure, B' = B/M,
A' = A/M, G'_i(TeX) = G_i^res/M. Hence, if B < M (so B' < 1),
  t_Z ≥ exp(min(inf_[0,1] A, min_i G_i^res)/M),
the middle term being unnecessary by step 6 (B' < 1 ⇒ A'(∞) = +∞). This is exactly what V 120-121 computes. Correct.
Value of the dropped middle term for the new certificate: −0.689682 (§4); it would not bind anyway.

## 3. Code audit: verify.py (V), lowcert.c (L), gcert.c (G)
Data path. JSON floats → Python binary64 (V 92) → "%.17g" (V 85-88, 103, 114-115) → strtod (L 141-146, G 41-42, 55):
17 significant digits round-trip binary64 exactly, so the exact Fraction checks (V 66-82) and both C programs see the
SAME doubles. Coefficients are JSON ints (max |a| = 4.27e16 > 2^53, so this matters): parsed as Python ints (V 92), passed
as decimal strings to fmpz_set_str (L 163, G 57). Exact.
V 68-77: edges strictly increasing, e_0 = 0, e_K = 1 (support ⊂ [0,1]); masses, multipliers ≥ 0; lc ≥ 1; distinct;
irreducible and primitive (python-flint factor, V 44-55); M, B as exact Fractions; assert B < M. All hypotheses of
thm:cert (after the §2 normalisation) are checked. V 108: [0,1/2] only if the exact symmetry test passed.
V 111: only c_j > 0 enter lowcert (c_j = 0 terms vanish from A). V 99-105: all 206 enter gcert; "MIN %.17g" of a
double floored from the arb lower bound (G 96-102) → exact double. V 120-121: min(V, gmin)/M in Fractions;
exp_lower floors x·10^40 then takes the lower end of mp.iv.exp: rigorous, except that the final mp.mpf() rounds to
nearest at 53 bits (≤ 3e-17 above the true lower end; REPORT Defect 4); irrelevant for the 10-digit claim.

lowcert.c. U_μ(x) = Σ_{k=0}^{K} β_k F(x−e_k), F(u) = u log|u| − u, β_k = ρ_k − ρ_{k−1}, ρ_{−1} = ρ_K = 0 (L 147-153, β at L 151):
correct, from ∫_a^b log|x−y| dy = F(x−a) − F(x−b) and summation by parts.
Roots (L 160-181): arb_fmpz_poly_complex_roots, rigorous, 160 bits; real roots have exactly zero imaginary part, the
non-real ones are stored once per conjugate pair (Im mid > 0) as −c log((x−a)²+b²); constpart = −Σ c log|lc| (L 169).
On I = [lo,hi], T = fl((lo+hi)/2) ∈ I (monotone rounding), Hb = max(hi−T, T−lo) as an arb ball (L 58-61):
 - near edge (L 69-77): F' = log|u| < 0 on 0 < |u| < 1 and u = x−e_k ∈ [−1, 1/2], so min_I β F(x−e_k) is at hi (β > 0),
   lo (β < 0), union of both (β ball ∋ 0). Correct.
 - far edge (L 78-82): Taylor data F(T−e), F' = log|T−e|, F'' = 1/(x−e) over the ball I. Correct.
 - near real root (L 88-90): −c log|x−r| ≥ −c log(|T−r| + Hb) since c ≥ 0. Correct.
 - far real root (L 91-97): f' = −c/(x−r), f'' = c/(x−r)². Complex pair (L 98-108): f' = −2c(x−a)/((x−a)²+b²),
   f'' = −2c(b²−(x−a)²)/((x−a)²+b²)². Both re-derived. Correct.
 - L 110-118: lb = val − |der|·Hb + min(0, inf curv)·Hb²/2 (Taylor, Lagrange remainder, ξ ∈ I), floored to double.
 - Safety: a "far" term whose singularity lies in I gives a non-finite ball → −1e300 → bisect (L 119). The double-precision
   near/far test (factor-2 slack) can only cost time, never soundness.
 - Ffun (L 31-50): ball ∋ 0 → |F(u)| ≤ r|log r| + r, r = sup|u| ≤ 0.36 < 1/e (r|log r| increasing there). Correct.
Coverage (lesson 277). L 187-192: lo_i = X0 + (X1−X0)·i/4096 and hi_i use the same expression with i+1 (identical
doubles, all exact dyadics for X1−X0 = 1/2); i = 0 gives exactly 0; last hi forced to X1 = 1/2. Bisection (L 126-133)
shares the midpoint. No gaps, starts at 0, ends at 1/2. Depth ≥ 45 → FAIL, rc 1 → V 119 assert.
"20,328 intervals" = all visited nodes (4096 roots + internal + leaves), not leaves; harmless wording.
Symmetry (checked, not assumed). Cell [e_k,e_{k+1}] ↦ [1−e_{k+1}, 1−e_k] = [e_{K−1−k}, e_{K−k}] (exact test e_k + e_{K−k} = 1),
mass m_k = m_{K−1−k} (exact), uniform density preserved ⇒ U_μ(1−x) = U_μ(x) exactly. σ(j) = index of ±Q_j(1−x) (V 57-63,
sign-normalised) is injective on distinct lc>0 polynomials, hence a bijection, with c_σ(j) = c_j (V 81) ⇒
Σ c_j log|Q_j(1−x)| = Σ c_j log|Q_j(x)|. So A(1−x) = A(x) exactly and inf_[0,1] A = inf_[0,1/2] A.

gcert.c. U_μ(α) = Σ β_k Re F(α−e_k) with principal log (G 66-81): valid because for Im α ≠ 0 the path u = α−y, y in a
cell, keeps Im u = Im α ≠ 0 and never meets the cut; real α uses real F. + M log a_i (G 84-86), exact integer
resultants with a zero-resultant abort (G 88-93), ÷ d_i (G 95), floor (G 96-97): this is G_i^res of §2. Correct.
Verdict on the code: no soundness defect found. Cosmetic: "intervals" counts nodes; a non-finite lb could in principle hit
the 1e300 sentinel only by bisecting further.

## 4. Independent recomputation (indep.py: mpmath/numpy, cell-form U_μ, Horner log|Q|, no arb, no resultants)
- exact (exact_new.log): support [0,1], edges strictly increasing, all 1832 masses > 0, c ≥ 0; M − 1 = 6.78248e-17,
  B = 0.999999998920136, B − M = −1.0799e-9 (exact rationals); symmetry exact (edges, masses, library ↔ reflected library
  with equal c, reflection by sympy substitution); all 206 polynomials irreducible, primitive, lc ≥ 1 (sympy, independent
  of python-flint); distinct; library = the paper's ancillary library (same 206, same order); −Σ c_j log a_j = −0.689682.
  V (JSON double) = −0.86219206599999997920; exp(V/M) = 0.42223549904149013751 ≥ 0.4222354990 (margin 4.15e-11).
- rand (rand_new.log): 10,000 seeded uniform points of [0,1], dps 30: min A = −0.862191199; none below V;
  |A(x) − A(1−x)| = 0 at 50 points.
- Grid scans are NOT enough (falsification lesson): uniform 2·10^5 grid → min −0.8621913316; edge-adapted 4.5·10^5
  points + 16 golden refinements → −0.8621914229. Both miss the true well found below (curvature ≈ 2.6e6 from density
  jumps β_k/(x − e_k) at edges 1.9e-5 away; the well is ~1e-6 wide).
- lowcert's tight leaves (bin/lowcert_dbg: same code + print of leaves with lb < V + 3e-9; 6 leaves, tight_leaves_new.txt)
  evaluated by indep.py: the leaves at x ≈ 0.3160865 contain a true local minimum (localmin_new.log)
    A = −0.862192064266183 at x = 0.316086571167473 (mirror 0.683913428832527, same value),
  i.e. A_min − V = +1.73e-9. Consistent with the rigorous inf A ≥ V and with NOTE.md ("fine-grid min −0.862192064");
  the certificate is tight: V cannot be raised by more than 1.7e-9 for this measure.
- gscan (gscan_new.log): 405,945 cell-adapted points (48 per cell + log-spaced offsets 1e-13..3e-4 on both sides of every
  edge + uniform 1e5), ALL 1789 sampled local minima refined (vectorised float64 golden section), lowest 8 re-refined in
  mpmath: GLOBAL MIN A = −0.862192064266183 at x = 0.316086571167473 (and mirror); only this pair lies within 1e-7 of V;
  none below V. Independent of lowcert, consistent with it.
- G (G_new.log, cmpG_new.log): all 206 G_i from mpmath roots (residuals ≤ 1.4e-41) and Horner, no resultants:
  mine − gcert ∈ [4.6e-20, 1.14e-16] for all 206 (gcert's floored doubles are ≤ 1 ulp below, never above);
  min G_i = −0.862191163250123 at i = 0 (x − 1; also x, then 2x − 1), so min G_i − V = 9.03e-7: V binds.
- Published certificate (exact_pub.log, gscan_pub.log): K = 947, not symmetric, 203 zero masses (allowed), library
  irreducible (sympy), M − 1 = −2.05e-17, B − M = −1.0e-9; independent global min A = −0.862208395409424 at x = 0.4036897,
  i.e. min A − V = 1.00e-8 (their V is the min minus 1e-8); −Σ c log a = −0.688183611561 (= TeX A(∞)); exp(V/M) =
  0.422228600019153 (M < 1 lowers exp(V) by 7.5e-18; TeX's e^V = 0.422228600019 unaffected at 12 digits). Reproduced.

## 5. Reviewer (delicate step: M-normalisation + dropped middle term + resultant form)
Pasted inline: thm:cert verbatim, deduction D (t_Z ≥ exp(min(inf_[0,1] A, min_i G_i^res)/M), no −Σ c_j log a_j term,
B < M strict), steps (i)-(v), the numbers. Verdict: CORRECT. Points: (iii) the factor is M − B + c_i d_i (with 1 in
place of M it would be wrong); M > 0 so division keeps the inequalities, and V/M > V for V < 0, M > 1; (iv) legitimate
through the proof (A' → +∞ uniformly at ∞ when B' < 1; superharmonic off [0,1]; the finite A(∞) is only for B' = 1);
c_j = 0 terms add no singularities; (v) min over all 206 members is safe. Session: reviewer 412ddd32 (2026-11-19).

## 6. Conclusion of the chain
thm:cert (correct; §1) + scaling (correct; §2, reviewer CORRECT) + sound arb code (§3) + rerun (verify_rerun.log) give
t_Z([0,1]) ≥ exp(V/M) = 0.42223549904149014 ≥ 0.4222354990 (truncated; margin 4.15e-11). Independent code (§4) finds
the true min A = V + 1.73e-9 and all G_i within 1.1e-16 of gcert: no counterexample to any step. Each of the paper's 23
support intervals in [0,1/2] carries exactly twice its cells (cells.log). Dependents that write 0.4222355 round UP
(REPORT.md, Defects 1).
