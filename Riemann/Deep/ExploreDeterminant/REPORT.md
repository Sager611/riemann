NO IDEA: the determinant method has no grip on μ(1/2) ≤ 13/84. The first spacing problem is already optimal. Huxley–Watt's second-spacing bound already equals the true count for the two conditions it uses (ratio 0.015–0.016 for t = 10^16–10^22, measured). The saving still missing (13/84 → 3/20, rate θ(σ) = (26+5σ)/(12(14+3σ))) lives on conics B_δ(4a,q) = 4(a'q' − aq), and each of these holds ≤ 2 counted pairs (box totals: median 1, max 3; Bombieri–Pila allows ≈ 100). So curve-by-curve determinant counts save nothing (node N9).

# explore49: the Bombieri–Pila / Heath-Brown determinant method → N9 (Orca exploration, 14 Nov 2026, 18:47–19:35 UTC)

## Short answer
- **No transfer.** Bourgain's 13/84 has two counting inputs:
  - The first spacing problem (Bourgain's A_6) is attained by major arcs, so no count can beat it.
  - The second spacing problem is the only non-optimal input. For ζ, the pairs it counts are integer points on conics (exact identity, §3). The generic conics hold ≤ 2 counted pairs, and HW's bound already has the true size for the conditions it uses. Any further saving needs equidistribution of √(t·aq/π) across conics. Determinant counts done conic by conic (or strip, hyperbola, x4-curve) do not see this.
- **The room is real (numerically).** The true count is at the "ideal" level σ = 2. Inside Bourgain's inequality it would give 3/20. The bottleneck is the method of proof, not the arithmetic.
- **By-products (exact bookkeeping; routine but not in the literature we checked):**
  - Suppose the off-diagonal second-spacing count is ≪ HW·H^{−σ}, uniformly near the optimal N (0 ≤ σ ≤ 2). Then |S| ≪ M^{1/2}T^{θ(σ)+ε} with **θ(σ) = (26 + 5σ)/(12(14 + 3σ))**, the same for all M ≥ T^{(6+σ)/(14+3σ)}.
  - In Huxley's resonance-curve route (x2 tightened by V = N/Q), the exponent is θ_V(s) = (32 + 5s)/(12(17 + 3s)). This gives 53/342 at s = 2/3, and **beating 13/84 there needs s > 3/4**.

## 1. Transfer (tool → node)
Tools: (d) the determinant step, together with (a) and (c), the determinant-method counts of integer points on and near curves (Bombieri–Pila 1989, Heath-Brown 2002, Huxley and Swinnerton-Dyer).

Node: N9, the record μ(1/2) ≤ 13/84 (Bourgain, arXiv 1408.5794 = JAMS 2017; the brief's "1407.2567" is a Witten-index paper). Bourgain §3 (TeX l. 760–950) runs: short intervals I, minor arcs a/q ∈ f''(I)/2 with Q ≤ q < 2Q, Poisson summation, then the sixth-power double large sieve:
- |S|⁶ ≪ M⁶/N³ + (M⁵R⁴N/Q⁷)·B1^{1/2};
- B1 = #{(I, I') : |x_j(I) − x_j(I')| < 1/(2Y_j), j ≤ 4} is the second spacing problem;
- Huxley–Watt 1988 (HW) gives B1 ≪ Δ1Δ2(M/N)²(Q/R)⁴ at N = MT^{−2/7}.

Claim tested: *a determinant-method count of the pairs in B1 beats HW (σ > 0) at the binding parameters Q = R, N = T^{1/5…3/14} (M = √T).* It does not, for any count that works curve by curve (§4).

## 2. What the record leaves open (exact; exponents.py, logs/exponents.txt)
Write all sizes as powers of T, with M = T^α and N = T^n. The three terms of |S|⁶ have exponents:
- trivial: 6α − 3n;
- diagonal (B1 ≥ |𝓘| ≍ M/N): α + 2n + 3/2;
- off-diagonal: 9α/2 − 3n/2 + 1/2 − (σ/4)(3n − 3α + 1).
In the V-route the diagonal exponent gains (3n − 3α + 1)/4.

| off-diagonal saving vs HW | route | exponent at M = √T |
|---|---|---|
| σ = 0 (HW, Bourgain) | V-free | 13/84 = 0.15476 ✓ (n = 3/14, Bourgain (3.13)) |
| σ = 1 (x1 or x4 fully used) | V-free | 31/204 = 0.15196 |
| σ = 2 (all four coordinates, observed true size, §3) | V-free | 3/20 = 0.15 (n = 1/5; trivial = diagonal = off-diagonal) |
| s = 2/3 (Huxley's Δ4^{2/3}) | V = N/Q | 53/342 = 0.15497 ✓ (Bourgain (3.18)) |
| s = 3/4 | V = N/Q | 13/84 exactly: the threshold |
| s = 2 (perfect counting) | V = N/Q | 7/46 = 0.15217 (ceiling for V = H^v: (6+v)/(40+6v)) |

All values were checked by an independent exact LP over n (fractions) and by the reviewer. The table assumes the improved bound holds uniformly at the re-optimized N: HW itself is quoted only at N = MT^{−2/7}, where it equals the diagonal.

## 3. Structure for ζ and exact counts (spacing2.c; logs/run0–run9; t = πK, M = √(t/2π)·t^{α−1/2}, Q = R, q odd)
**The four coordinates.** Gauss sum plus stationary phase give
- x2 = −\bar{4a}/q;
- x1 = −b·\bar{2a}/q = 2b·x2 (mod 1);
- x3 = (2/3)/D and x4 = η/D, with b + η = √(t·aq/π) (b the nearest integer) and D ∝ t^{−1/4}(aq)^{3/4}.
So x3 and x4 depend only on P = aq, and x1 and x4 are the integer part and the fractional part of √(tP/π). This is the system of Bourgain–Watt 1505.04161 §10. My normalization fixes the constants below; it does not affect the exponents.

**Conic identity (exact algebra, reviewer-checked).** Let A = 4a, γ = [[A, c], [q, d]] ∈ SL₂(Z), and δ = γ'γ^{−1} = [[α, β], [g, e]], so (A', q') = δ(A, q) and ‖\bar{A}/q − \bar{A'}/q'‖ = ‖g/(qq')‖. Then
- Aq − A'q' = −B_δ(A, q), with B_δ = αg x² + 2βg xy + βe y² and disc = −4βg;
- so the (x2, x3) pairs are integer points on the conics B_δ = 4(a'q' − aq) with |g| ≲ Q²/H²;
- the identity was verified on every recorded pair (≤ 10⁶ per run), with 0 failures.

**Smallest case**, t = 10¹⁰: H = 2.04, 77 arcs, 12 off-diagonal (x2, x3) pairs on 12 distinct conics. Example: 41/51 ~ 44/43, δ = [[169, −540], [41, −131]], Aq − A'q' = 796 = −B_δ(164, 51).

**At scale** (α = 1/2; full table logs/table_final.md):
| log10 t | ν | H | arcs | (x2,x3) pairs | ÷HW shape | all-4 pairs | ÷(shape·H⁻²) | true/random (x2x3, all4) | off/diag |
|---|---|---|---|---|---|---|---|---|---|
| 16 | 0.3 | 2.51 | 6027 | 85250 | 0.0152 | 842 | 0.00095 | 1.21, 1.32 | 0.140 |
| 18 | 0.3 | 3.16 | 24310 | 555684 | 0.0159 | 3576 | 0.00102 | 1.23, 1.36 | 0.147 |
| 20 | 0.3 | 3.98 | 97028 | 3526664 | 0.0160 | 14604 | 0.00105 | 1.24, 1.39 | 0.150 |
| 22 | 0.3 | 5.02 | 386246 | 22239878 | 0.0161 | 58108 | 0.00106 | 1.24, 1.38 | 0.150 |
| 22 | 2/7 | 14.87 | 187170 | 66416 | 0.0157 | 20 | 0.00104 | 1.21, 1.9 | 1.1e−4 |

Observations (not theorems):
- (x2, x3) ≈ 0.016 × HW shape for every t, ν, V ∈ {1, N/Q} and α ∈ {0.45, 0.5}: HW is sharp for what it uses.
- all-4 ≈ 0.00105 × shape × H^{−2}, i.e. σ = 2. The truth beats Huxley's Δ4^{2/3} shape by a growing factor (0.00066 → 0.00034).
- The ratio to a random model saturates at about 1.24 and 1.38, so there is no growing arithmetic excess.
- At N = T^{1/5} the off-diagonal count is 0.15 of the diagonal, which is the predicted balance.

## 4. Why determinant counts save nothing here (what is covered, what is not)
1. **Fibres hold ≤ 2 points.**
   - Generic conics (disc > 0 and not a square) hold ≤ 2 counted (x2, x3) pairs (α ∈ {0.42, 0.45, 0.5}); a fixed δ holds ≤ 12 pairs (≤ 41 in the V-route).
   - Over the whole box, 3870 sampled conics hold median 1 and max 3 integer points (logs/run7, run8). Bombieri–Pila allows √(8Q) = 90–124.
   - A bound summed over these fibres gives an integer ≥ 1 to each fibre that can carry a point. So it is ≥ #(x2, x3) pairs/2 ≍ HW, which means σ = 0, whatever the degree (one vanishing determinant needs (d+1)(d+2)/2 > 2 points on the curve).
   - The same holds for strips of fixed δ and for hyperbolas a'q' = P' at fixed (a, q).
2. **x4 by divided differences (degree 1).**
   - In P, the curve y = √(tP/π) has g'' ∈ [0.0156H, H] (Q = R). For distinct integers, det[1, P_i, m_i] ≥ g'' − 4ε, so for H ≳ 64 no Huxley or Swinnerton-Dyer determinant vanishes. Shears do not change g''.
   - This does *not* exclude higher degree: y² = (t/π)P is itself a conic. The x4 points lie at distance ≍ 1/H from it (the thick regime), and we did not bound degree-d auxiliary curves there.
3. **Fixed n = P − P'** is the flat chart (curvature ≈ nH/Q²). But the x2-admissible P form an arithmetic set of density 10^{−5}–10^{−7} per n (runs). A rounding count over all P therefore gives up the x2 saving and ends far above HW.
4. **p-adic.** x4 is archimedean. x1 and x2 are congruences to the varying modulus q, so grouping mod a fixed prime re-sorts the same fibres.
5. **Degenerate families are lines** (unchanged in the V-route): rational line pairs (square discriminant) and, only for α < 1/2, translations g = 0 (q' = q, a' = a + jq).
   - A 1-D rounding count along such a line can save x4 (curvature ≈ 0.15j at α = 0.42).
   - They carry 0 (α = 1/2), 2.1% (0.45) and 32% (0.42) of the pairs, with share ∝ T^{4/5−2α}.
6. **Not excluded:**
   - regrouping across fibres, or other coordinates;
   - joint congruence and archimedean arguments;
   - showing that most fibres carry zero full-count pairs;
   - an exponential-sum count of x4 with s > 3/4 in the V-route (Huxley's resonance curves, s = 2/3).

## 5. Adversarial
- The counting problem depends only on the phase t·log m. Splitting into residue classes mod 5 turns the pieces Σχ(n)n^{−it} of the Davenport–Heilbronn function into sums Σ_m e(−(t/2π)log(5m + r)) of the same shape. So every Bombieri–Iwaniec exponent (and Bourgain's exponent pair (13/84, 55/84)) holds for Davenport–Heilbronn on its critical line, and so does this obstruction or any improvement of it. N9 bounds carry no RH information.
- Shuffling x1, x2 and η independently reproduces the counts up to the constant 1.24–1.38.

## 6. Novelty
**Searches.**
- arXiv listings 2026-10/11 (math.NT, math.CA, math.CO), with a title regex for spacing, Bombieri, Iwaniec, determinant method, Pila, Huxley, exponent pair, subconvex, 13/84, points near.
- arXiv API queries: "second spacing problem", "Bombieri-Iwaniec", "determinant method" with zeta, "Bombieri-Pila" with zeta.
- Zenodo: "Bombieri-Iwaniec", "spacing problem", "13/84".

**What they return.**
- Rational points near curves: 2610.02443, 2610.23797, 2610.03231.
- Determinant methods for Galois groups: 2610.35030.
- The circle problem: Cen–Gao–Zhang–Zhang 2610.21730 (209/670, using Watt's second-spacing bound; nothing for ζ).
- Li–Yang 2308.14859 (first spacing).
- Zenodo: nothing relevant.

**Not claimed new.**
- The SL₂ matrix device: Bombieri–Iwaniec 1986 count matrices; their paper was not checked here.
- The 3/20 limit: derived here from Bourgain's inequality; it may appear in Huxley's book (unchecked).

**Not seen elsewhere:** the measured counts, θ(σ), θ_V(s) and the threshold s > 3/4. These are bookkeeping and data (RECORDS/COROLLARIES level), not a discovery.

## 7. Other bridges, measured and rejected
- N7/N10 (large values): the large-value points t_r are real numbers. Discretized as (round(t_r log p/2π))_p, they are integer points near a line, where determinant counts are trivial. Large-value bounds come from Fourier and energy arguments (explore4, explore46).
- N8/N12: no step counts lattice points near curves.

## Suggested NOGO line
- explore49 (Orca, exploration: the Bombieri–Pila / Heath-Brown determinant method, integer points on and near curves → N9) no transfer: 13/84's only non-optimal input, the Bombieri–Iwaniec second-spacing count, is for ζ a count on conics B_δ(4a,q) = 4(a'q' − aq) with ≤ 2 counted pairs each (box totals median 1, max 3, vs Bombieri–Pila ≈ 100); Huxley–Watt already equals the true count for the conditions it uses (ratio 0.016, t ≤ 10^22, measured); the missing saving toward 3/20 (rate θ(σ) = (26+5σ)/(12(14+3σ)); Huxley's V-route needs s > 3/4) is cross-conic equidistribution of √(taq/π), which fibrewise determinant counts cannot see [ExploreDeterminant](Riemann/Deep/ExploreDeterminant/REPORT.md)

## Files
- Code: spacing2.c (exact counts; usage in the header), conic_points.py, exponents.py, table.py, DERIVATION.md.
- logs/:
  - run0: smallest case;
  - run1, run2, run9: scaling;
  - run3, run8: V-route;
  - run4, run6: degenerate families;
  - run5: α < 1/2;
  - run7: points per conic;
  - exponents.txt, table_final.md.
- Every run was under memcap 290 MB; the peak was 70 MB.
- Reviewer: algebra (θ, θ_V, conic identity, x1 = 2b·x2) correct. The "barrier" claim was too strong and was downgraded to this fibrewise obstruction.

RH relevance: N9 (obstruction for subconvexity via Bombieri–Iwaniec; it does not distinguish ζ from Davenport–Heilbronn).
