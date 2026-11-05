# PROOF: the half-plane does not lower the 0.856 threshold (nor any A(σ) entry on [3/4, σ*))

Node N7. Companion files: verify.py (exact rationals, standard library, 0.01 s), explore2.py (numerical
cross-check, about 10 s per σ). Paper arguments, not refereed, no Lean. The result is negative (a no-go for the
channels listed in §2); it does not exclude a genuinely new averaging argument (§8).

## 0. Setting

B := 0.87338196661320752, (H): ζ(s) ≠ 0 for Re s > B (Lean-checked; see DensityB PROOF.md §0).
τ_B(σ) := (2σ − 2 + B)/(B − 1/2), g_B(τ) := 1 − B + (B − 1/2)τ, σ* := (5 − 4B)/(7 − 6B) = 0.856092011887…
(τ_B(σ*) = 3σ* − 1, verify.py first line). Notation of Tao–Trudgian–Yang arXiv 2501.16779 = ANTEDB blueprint
(github.com/teorth/expdb, main, read 2026-11-05): large value pattern (N, T = N^τ, V = N^σ, a, J, W), |W| = N^ρ;
LV, LV_ζ; double zeta sum S(N, W) := Σ_{t,t'∈W} |Σ_{n∈[N,2N]} n^{−i(t−t')}|². "MC at (σ, τ)" means
LV(σ, τ) ≤ 2 − 2σ. By Cor. zero-large-cor3, MC on [0, τ_M] plus the zeta condition (lvoz) give
A(σ) ≤ 3/(τ_M + 1 − σ).

Everything (H) supplies about ζ (DensityB §1 and Lemma 1, plus (Z1)):
- (Z1) short zeta sums: for I ⊂ [N, 2N] and N ≤ |t|, Σ_{n∈I} n^{−it} ≪ N^B |t|^ε ("pair" (0, B), zeta sums only).
  Proof: truncated Perron at height |t|/2, contour moved to Re w = B, where ζ(B + i·) ≪ |t|^ε (DensityB (C1));
  horizontal sides and Perron error ≪ N^{1+ε}/|t| ≤ N^B. ∎
- (Z2) long zeta sums (DensityB Lemma 1): for πN ≤ |t|, Σ_{n∈I} n^{−it} ≪ N^{1−B}|t|^{B−1/2+ε} ("pair"
  (B − 1/2, 1/2)). (Z1) is the smaller bound for N ≤ |t|^{1/2}, (Z2) for N ≥ |t|^{1/2}.
- (Z3) hence LV_ζ(σ, τ) = −∞ for σ > B (τ ≥ 2) and for σ > g_B(τ) (1 ≤ τ ≤ 2).
- (Z4) moments of ζ on (1/2, B) by convexity (DensityB (C3)); they give no new LV_ζ bound (DensityB Prop. 2).
These are sup-norm bounds at the level N^B, far above the mean-square level N^{1/2}; on the critical line they are
weaker than Bourgain's μ(1/2) ≤ 13/84 (T^{B/2} = T^{0.437} versus T^{1/4+13/84} = T^{0.405} for a sum of length
T^{1/2}).

## 1. The records on [3/4, 0.86] and their binding inputs

| σ range | A(σ) ≤ | source (ANTEDB/TTY Table 3) | τ0 in cor2/cor3 | binding input |
|---|---|---|---|---|
| [0.75, 0.76) | 15/(3+5σ) | Guth–Maynard 2024 | (3+5σ)/5 | GM large-value estimate max(2−2σ, 18/5−4σ, τ+12/5−4σ): Poisson, Hölder-type S_2, additive-energy S_3 (Lemmas gm-1, gm-3) |
| [0.76, 0.7604) | 9/(8σ−2) | Ivić, Lemma gzd, m = 6 | τ0_6 | diagonal 6(τ+2−4σ); start c_6 = 47/62 from the double-zeta term X |
| [0.7604, 0.7647) | 15/(13σ−3) | gzd, m = 5 | τ0_5 | same; c_5 = 127/167 |
| [0.7647, 0.7727) | 6/(5σ−1) | gzd, m = 4 | τ0_4 | same; c_4 = 13/17 |
| [0.7727, 0.7735) | 2/(9σ−6) | TTY, improved Bourgain | 9(3σ−2)/2 | Jutila k = 3 + Bourgain's large-value theorem (double zeta sums at length N) |
| [0.7735, 0.7777) | 9/(7σ−1) | gzd, m = 3 | τ0_3 | diagonal 3(τ+2−4σ); c_3 = 41/53 |
| [0.7777, 0.7954) | 9/(8(2σ−1)) | TTY, improved Bourgain | 8(2σ−1)/3 | as two rows up |
| [0.7954, 0.8) | 3/(2σ) | Bourgain 2000 | 2σ | Bourgain's theorem; its R^(1) terms are Jutila's k = 2 terms 2τ+4−8σ, τ+12−16σ |
| [0.8, 0.875) | 3/(2σ) | Ivić, gzd, m = 2 | 2σ | diagonal 2(τ+2−4σ), i.e. the Montgomery cap 3σ − 1; Type II: twelfth moment (slack) |

Unpublished and excluded from the table: Kerr, 3/(2σ) for σ ≥ 23/29.
Lemma gzd (Ivić (11.74)–(11.77), ANTEDB Lemmas a-ivt-1, gzd): LV ≤ max(2−2σ, m(τ+2−4σ), min(X, Y)) with
X = 2τ/3 + 4m(3−4σ)/3 and Y = max(τ + 3m(3−4σ), ((k+ℓ)τ + m(1+2k+2ℓ)(3−4σ))/k) for an exponent pair (k, ℓ);
Type II by Heath-Brown's twelfth moment LV_ζ ≤ 2τ − 12(σ−1/2); then A ≤ 3m/((3m−2)σ+2−m) for
σ ≥ min(c_m, max(a_m, g_m(k,ℓ))), c_m = (6m²−5m+2)/(8m²−7m+2), a_m = (9m²−4m+2)/(12m²−6m+2).
Where ζ enters: (i) the diagonal term m(τ+2−4σ) (mean value of the m-th power of the dual zeta sums), which sets
the cap τ_m = 4σ − 2 + (2−2σ)/m = τ0_m + σ − 1; (ii) X, Heath-Brown's off-diagonal double-zeta term, which sets
the start c_m; (iii) Y, a pointwise exponent-pair bound for zeta sums, used only to push the start below c_m
(Ivić: 3831/4791 for m = 2); (iv) Type II via the twelfth moment. verify.py A1–A3 recompute the five formulas,
the caps and the starts.

## 2. Lemma R: the common mechanism (reflection principle, corrected)

Lemma R (Ivić (11.58)–(11.59); ANTEDB Theorem "Reflection principle"). Let σ ≥ 3/4, 1 < τ ≤ 2, k ≥ 1, and let a
large value pattern have T = N^τ, V = N^σ, |W| = N^ρ. Then ρ ≤ 2 − 2σ, or for some 0 < α ≤ k(τ − 1)
    ρ ≤ k(3 − 4σ)/2 + (s' − α)/2,   S(N', W) = N^{s'+o(1)}, N' = N^{α+o(1)}.
Proof. Halász: |W|²V² ≤ N Σ_{t,t'} |Z(t − t')|, Z(u) = Σ ψ(n/N) n^{iu}. Pairs with |t − t'| ≤ N^{1−ε} give
|W| N^{2−2σ}. For |u| ≥ N^{1−ε} the B-process gives |Z(u)| ≲ N^{1/2} |D(u + v)| with D(u) = Σ_{m≤4T/N}
m^{−1/2+iu} and some |v| ≤ T^{o(1)}; so |W|² ≲ |W| N^{2−2σ} + N^{(3−4σ)/2} Σ_{t≠t'} |D(t − t')|. Hölder:
Σ |D| ≤ |W|^{2−1/k} (Σ_{t,t'} |D^k(t − t')|²)^{1/(2k)}. Expanding, Σ_{t,t'} |Σ_n b_n n^{−i(t−t')}|² =
Σ_{n,m} b_n b̄_m |R_W(m/n)|² with |R_W|² ≥ 0, so it is ≤ max|b_n|² · S(N', W)/N' after a dyadic split
(b_n = d_k-type coefficients times n^{−1/2}, of size N'^{−1/2+o(1)}). Solving for |W| gives the claim. ∎
Remark. The blueprint prints "ρ ≤ min(2 − 2σ, k(3−4σ)/2 + s' − 1)"; with the diagonal s' = ρ + 2α that
inequality does not involve ρ (verify.py B7), so it is a typo; the form above reproduces Jutila's and Heath-Brown's
terms (cases (a), (b) below). The ANTEDB Python code does not implement this theorem.

Corollary R. Insert the bounds for s':
- (a) diagonal: S(N', W) ≥ |W| N'² exactly (the terms t = t'), and Heath-Brown gives s' ≤ ρ + 2α in this case:
  ρ ≤ k(3 − 4σ) + α ≤ k(τ + 2 − 4σ); MC iff τ ≤ τ_k := 4σ − 2 + (2 − 2σ)/k (k = 2: 3σ − 1).
- (b) Heath-Brown's off-diagonal s' ≤ 5ρ/4 + τ/2 + α: ρ ≤ X = 2τ/3 + 4k(3 − 4σ)/3; MC iff
  τ ≤ 3 − 3σ − 2k(3 − 4σ) (k = 1: 5σ − 3; k = 2: 13σ − 9).
- (c) Heath-Brown's s' ≤ 2ρ + α: impossible for σ > 3/4.
- (d) new, from (H): s' ≤ 2ρ + h(α) with h(α) = 2Bα for α ≤ τ/2 by (Z1), h(α) = 2(1−B)α + (2B−1)τ for
  α ≥ τ/2 by (Z2). Then ρ ≤ k(3−4σ)/2 + ρ + (h(α) − α)/2, a contradiction iff (2B−1)·min(α, τ−α) < k(4σ−3).
  The worst dyadic length is N' = T^{1/2}, so the case is excluded for all α iff
  (2B − 1)·min(k(τ−1), τ/2) < k(4σ − 3). For k = 1 (τ ≤ 2): iff τ < τ_B(σ) (verify.py B2).

## 3. Proposition 1 (the binding diagonal does not depend on ζ)

Case (a) uses only S(N', W) ≥ |W| N'², true for every 1-separated W; no bound on zeta sums enters, so the caps
τ_k (k ≥ 2) are the same with or without (H). It is also not an artefact of including t = t': pairs with
1 ≤ |t − t'| ≤ 2 have |D(t − t')| ≍ X^{1/2} (X = T/N), so after Hölder they contribute just as much, and a large
value set may contain such pairs. Hence on [4/5, σ*) the cap 3σ − 1 of every argument built on Lemma R (Jutila
k = 2, Heath-Brown's Lemma 6, Ivić (11.74), the k = 2 terms in Bourgain's R^(1)) is (H)-independent.

## 4. Proposition 2 (every sup-norm channel of (H) yields exactly τ_B)

- (a) Halász–Montgomery with (Z2) or (Z3) (Montgomery's LV theorem, Lemma hl-improv): MC for τ < τ_B (DensityB
  Thm 1).
- (b) Lemma R with (d): for k = 1 exactly τ < τ_B (with Corollary R(b), MC for τ < max(τ_B, 5σ − 3), capped by
  2σ from Corollary R(a)).
  For k ≥ 2, if k(τ − 1) ≤ τ/2 the condition is again τ < τ_B; otherwise it reads (2B − 1)τ < 2k(4σ − 3), but
  then case (a) caps at τ_k ≤ 3σ − 1. So Lemma R with all inputs gives MC exactly for
  τ ≤ max(τ_B, 3σ − 1) on [4/5, σ*) (verify.py B5), the old range.
- (c) Corollary mu-lv (MC for τ ≤ (2σ − 1 − θ)/μ(θ)) with the new hull point μ(1 − B) = B − 1/2 (and its
  convex combinations): on each segment of the convex hull the ratio is linear-fractional, so its maximum is at a
  vertex; a vertex is either an old point (old range, no gain) or (1 − B, B − 1/2), which gives exactly τ_B; on
  [0, 1 − B], μ(θ) = 1/2 − θ and the ratio increases in θ for σ > 3/4 (verify.py B1 compares θ = 1 − B with
  θ = 1/2).
- (d) Lemma hl-improv beyond τ_B: one needs LV_ζ(σ', τ') ≤ 1 − σ' on the sliver 2σ − 1 ≤ σ' ≤ g_B(τ').
  Reflection (Lemma lvz-basic(iv)) and (Z3) give LV_ζ(σ', τ') = (τ' − 1) sup_{σ_r ≤ σ'' ≤ B}
  (LV_ζ(σ'', τ'') + σ'' − σ_r) with τ'' = τ'/(τ' − 1), σ_r = 1/2 + (σ' − 1/2)/(τ' − 1). The best known
  input, LV_ζ ≤ LV(σ'', τ'') ≤ 2(2 − 2σ'') (Montgomery at τ''/2, squared), gives LV_ζ(σ', τ') ≤ 2τ' − 4σ',
  which is ≤ 1 − σ' at σ' = 2σ − 1 iff τ' ≤ 3σ − 1. So hl-improv returns max(τ_B, 3σ − 1): the reflection maps
  the k = 2 cap onto itself. Going beyond needs LV_ζ(σ'', τ'') < 2(2 − 2σ'') for short zeta sums
  (τ'' ≈ 2.5–3.5) at levels σ'' just below B, where (H) allows values up to N^B and gives no count bound.

## 5. Proposition 3 (interpolations are affine)

- (a) Hölder between k = 1 and k = 2 in Lemma R (Σ|D|^p ≤ (Σ|D|²)^θ (Σ|D|⁴)^{1−θ}, p = 4 − 2θ): for fixed
  (ρ, α1, α2) the bound is ρ ≤ (2 − θ)(3 − 4σ)/2 + [θ(s1 − α1) + (1 − θ)(s2 − α2)]/2, affine in θ.
- (b) Heath-Brown-type bound with the twelfth moment, moved off the critical line: Perron at Re w = 1 − s
  (1/2 ≤ s ≤ B, θ = (s − 1/2)/(B − 1/2)), |Z(u)| ≲ N^{1−s} T^{s−1/2} |ζ(s + iu')|, Hölder with
  ‖ζ(s + i·)‖_q ≤ T^{2/q+ε}, 1/q = (1 − θ)/12 (Gabriel interpolation of M(1/2, 12) ≤ 2 with the sup bound at B).
  The Montgomery condition becomes 12[3/2 − 2σ + θ(B − 1/2)(τ − 1)] + (1 − θ)(2τ + 2σ − 2) ≤ 0, affine in θ;
  θ = 0 gives Heath-Brown's τ ≤ 11σ − 8, θ = 1 gives τ ≤ τ_B. With the fourth moment: 6σ − 4 and τ_B.
- (c) Pointwise-times-mean hybrids |D|^{2k} ≤ (sup|D|)^{2k−2j}|D|^{2j}: the exponent is linear in k for fixed j,
  with slope negative exactly when τ < τ_B.
An affine bound is optimised at an end point, so none of these beats max(τ_B, the old range).

## 6. Proposition 4 (exponent-pair channels)

- (a) Lemma gzd with the zeta pair (k, ℓ) = (B − 1/2, 1/2) (we do not check that Ivić's proof of (11.74) accepts a
  zeta-only pair; even granting it): m = 2 starts at g_2 = (23 + 18(B − 1/2))/(29 + 22(B − 1/2)) =
  0.798639… > 1867/2347 = 0.795483 (Bourgain already has 3/(2σ) there); for m = 3…6, c_m ≤ max(a_m, g_m), so
  the starts c_m are unchanged (verify.py A4). (The Lindelöf pair (0, 1/2) would give 23/29, Kerr's value.)
- (b) Ivić's zeta-from-exp (LV_ζ ≤ max(τ − 6(σ − 1/2), ((k+ℓ)τ − 2(1+2k+2ℓ)(σ − 1/2))/k)) with the zeta pair
  changes only Type II, which is not binding (Prop. 6).
- (c) Exponent pairs inside Halász–Montgomery: Prop. 2(a).

## 7. Propositions 5 and 6 (double zeta sums; Type II)

Prop. 5. Heath-Brown: s ≤ max(ρ + 1, 2ρ, 5ρ/4 + τ/2) + 1 for S(N, W) = N^s. For τ ≤ 3/2 and ρ ≤ 1 this equals
ρ + 2, the diagonal, hence optimal; at every binding point of §1 (τ ≤ 3σ* − 1 = 1.5683, ρ = 2 − 2σ ≤ 1/2) also
5ρ/4 + τ/2 ≤ ρ + 1 (verify.py B4). (H) cannot beat a diagonal-optimal bound (pointwise it gives
2ρ + 2 − 2B + (2B − 1)τ ≥ 2ρ + 1). So the inputs S(N, ·) of Bourgain's theorem, of Lemmas gm-1, gm-3 and of
Heath-Brown's energy relations are unchanged. (H) can only improve Heath-Brown's bound where his |W|^{5/4}T^{1/2} term dominates,
i.e. ρ > 4 − 2τ; for ρ ≤ 1/2 that needs τ > 7/4, beyond every binding point.
Prop. 6. The Type II condition (lvoz) holds with room on every row of §1 by the twelfth moment (verify.py C1,
C2; for 3/(2σ) the margin at τ = 8σ/3 is 8σ/3 − 2 ≥ 0.121), so lowering LV_ζ by (Z3) or Prop. 4(b) changes no
entry.

## 8. Conclusion and what would be needed

Theorem (no-go for these channels). Let σ ∈ [3/4, σ*). Adding (Z1)–(Z4) to the large-value and zeta-large-value
estimates of §1, through Halász–Montgomery (pointwise or with LV_ζ), Corollary mu-lv, Lemma R (any k, with
Heath-Brown's or the pointwise double-zeta bounds, per pattern the better one), Hölder or Gabriel interpolation,
Ivić's exponent-pair lemmas, and the double-zeta inputs of Bourgain's and Guth–Maynard's theorems, gives no bound
on A(σ) below the table value. On [4/5, σ*) the Montgomery range stays max(τ_B, 3σ − 1) = 3σ − 1, so σ* is not
lowered. Below 4/5 the half-plane channel 3/(τ_B + 1 − σ) is above every table entry (verify.py A5).
For Guth–Maynard (row 1) we did not re-derive Props. 6.1–11.1 of their paper. At τ = 6/5, σ = 3/4 the binding
inequalities in ANTEDB's proof are their S_3 bound (Prop. 10.1) with the energy bound E_1(W) ≪ |W|N^{4−4σ} + …
(Prop. 11.1); ANTEDB states these as inequalities in |W|, E_1(W), N, T with no zeta-sum exponent, and their
double-zeta input (second Guth–Maynard relation with Heath-Brown) is diagonal-optimal (Prop. 5). So the GM row
rests on this reading and on explore2.py, not on a full proof.

The binding inputs, in a form that names what a new argument must supply:
1. k = 2 (cap 3σ − 1): a bound for Σ_{t≠t'∈W} |D(t − t')| (D the dual zeta sum of length X = T/N) that does not
   pass through Hölder over all pairs; the fourth moment over W − W is ≥ |W| X² from pairs at distance ≤ 2.
2. k = 1 (its diagonal allows the Montgomery range up to 2σ; Type II would then bind): the off-diagonal double
   zeta sum of short sums,
   X^{−1} Σ_{t≠t', |t−t'|≥X^{1+ε}} |Σ_{m∼X} m^{−i(t−t')}|² ≪ N^{1−δ} for |W| = N^{2−2σ}, X = N^{τ−1},
   τ slightly above 3σ − 1. Heath-Brown's |W|^{5/4}T^{1/2} term gives it for τ ≤ 5σ − 3, (H) for τ < τ_B;
   "Lindelöf on average over difference sets" would give it up to τ = 2σ.
3. hl-improv: LV_ζ(σ'', τ'') < 2(2 − 2σ'') for zeta sums of length T^{1/τ''} (τ'' ≈ 2.5–3.5) at levels just
   below B, i.e. better than the squared Montgomery bound for short zeta sums near the half-plane level.

## 9. Numerical cross-check (explore2.py)

Own pipeline (DensityB's explore.py plus: Ivić's zeta-from-exp corrected to a maximum, the fourth moment,
Ivić (11.74) with m = 2…6 and four classical pairs; new mode adds the zeta pair and DensityB's inputs).
Corollary zero-large-cor-0 on a τ-grid of step 0.004 (grid maxima can undershoot, so the absolute values may sit
up to about 0.002 below the table; the old/new difference is the test):

| σ | table | pipeline old | pipeline new | new − old |
|---|---|---|---|---|
| 0.750 | 2.22222 | 2.22222 | 2.22222 | 0 |
| 0.760 | 2.20588 | 2.20588 | 2.20588 | 0 |
| 0.770 | 2.10526 | 2.10379 | 2.10379 | 0 |
| 0.780 | 2.00893 | 2.00803 | 2.00803 | 0 |
| 0.790 | 1.93966 | 1.93895 | 1.93895 | 0 |
| 0.800 | 1.87500 | 1.87500 | 1.87500 | 0 |
| 0.820 | 1.82927 | 1.82927 | 1.82927 | 0 |
| 0.840 | 1.78571 | 1.78571 | 1.78571 | 0 |
| 0.850 | 1.76471 | 1.76471 | 1.76471 | 0 |
| 0.856 | 1.75234 | 1.75234 | 1.75234 | 0 |
| 0.858 | 1.74825 | 1.74825 | 1.74419 | −0.00407 |
| 0.860 | 1.74419 | 1.74419 | 1.73410 | −0.01008 |
| 0.870 | 1.72414 | 1.71821 | 1.69109 | −0.02712 |

Full run (22 values of σ in [0.75, 0.87]): new − old = 0 for every σ ≤ 0.856, negative above σ*, as DensityB.
Side fix: DensityB's repro/explore.py combined the two terms of Ivić's zeta-from-exp with min instead of max
(too optimistic for Type II). Test (.scratch run, DensityB's script with only that line changed): identical
old/new values at σ = 0.80, 0.84, 0.86, 0.87, so DensityB's cross-check stands; at σ = 0.78 the buggy version
returns 2.00000, a spurious value below the table's 2.00893 (corrected: 2.00803, grid error).

## 10. Self-review

Pass 1 (directions, ranges). MC ranges are intervals [0, τ_M]; Lemma R needs σ ≥ 3/4, 1 < τ ≤ 2, which holds
at all binding points (τ ≤ 1.57). The worst dyadic length in (d) is checked over all 0 < α ≤ k(τ − 1). "Per
pattern the better one" (min of Heath-Brown and (H) in s') is allowed because both bound the same quantity; at
ρ just above 2 − 2σ the two exclusions need τ ≤ 5σ − 3 and τ < τ_B respectively, so the minimum adds nothing.
Pass 2 (what is not proved). This is a no-go for the listed channels, not a theorem that (H) cannot lower A
below σ*. The Guth–Maynard row relies on a structural reading (§8). Admissibility of a zeta-only pair in Ivić's
(11.74) is not checked; it would not matter (Prop. 4(a)). The numerical pipeline is a cross-check, not a proof.
