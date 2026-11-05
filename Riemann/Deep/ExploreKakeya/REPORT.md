NO IDEA: every Kakeya/multiscale transfer (sticky vs non-sticky large-value sets, multiscale additive energy, polynomial partitioning of the spikes) acts on Guth–Maynard's additive-energy and Heath-Brown inputs, and inside their cubic-trace chain these cannot move the density exponent 30/13: it is pinned at (σ, T) = (7/10, N^{6/5}). That point lies on Guth's complexity-barrier line σ = 1 − τ/4, the same line below which Bourgain–Kakeya fat-progression polynomials are genuine examples. Even perfect versions of both inputs only lower A(3/4) from 2.222 to 2.180, and the prime-modulus gap (7/3) has no intermediate scales to induct on.

# Exploration explore4: Kakeya (Wang–Zahl) and multiscale harmonic analysis → RH tree (2026-11-05)

Field: Wang–Zahl's Kakeya in R³ (sticky case, convex Wolff axioms, induction on scales); Bourgain–Kakeya; Guth–Maynard (GM) large values; Heath-Brown's (HB) difference-set bound; Chen–Gupta–Li (CGL); Tao–Trudgian–Yang (TTY). Brief: .scratch/orca/explore4-brief.md. Notation: |b_n| ≤ 1, D(t) = Σ_{n∼N} b_n n^{it}, W ⊂ [0,T] 1-separated with |D| ≥ V = N^σ on W, T = N^τ, E(W) the additive energy of W, A(σ) the zero-density exponent.

## Survey (one page)
- **Kakeya.** Wang–Zahl 2210.09581 (with a sequel) prove the sticky case in R³ (Katz–Tao planiness and graininess, plus Furstenberg-set inputs). 2502.17655 proves the full case: δ-tubes that do not cluster in convex sets (Δ_max ≲ 1) have union volume ≳ |T|·|T|. The general case reduces to the sticky case by maximal-density factoring and induction on scales (Guth's outline 2508.05475). Follow-ups: sticky R⁴ under convex Wolff axioms (2609.22035), curved sticky Kakeya (2511.10918), and a non-sticky null Kakeya set (2506.18142).
- **The bridge.** Bourgain ([B5] in Guth's survey): Montgomery's conjecture implies the "Bourgain–Kakeya" conjecture, which implies Kakeya. The Bourgain–Kakeya conjecture says that fat arithmetic progressions I* + t_I (blocks I of length N T^{−1/2}, spacing 2πn_I, thickness T^{1/2}) cannot overlap heavily. Guth's survey 2503.07410 makes three points. §9: "assuming the Bourgain–Kakeya conjecture as a black box does not currently lead to any improvements", and it is unclear what Wang–Zahl, Dvir or Arsovski methods say about it. §8.1: for σ < 1 − τ/4 (7/10 at τ = 6/5) no large-value estimate beats the basic bound, and for random matrices no low-degree or sum-of-squares certificate is believed to exist there.
- **Large values.** GM 2405.20552, Thm 1.1: |W| ≲ N^{2−2σ} + N^{18/5−4σ} + T N^{12/5−4σ}, giving A ≤ 15/(3+5σ) on [7/10, 8/10]; 30/13 is where this meets Ingham's 3/(2−σ) at σ = 7/10. Their dichotomy: low energy is handled by the cubic Schatten trace and the "spikes" of the difference set at rationals p/q with p, q ≲ T/N (Props 8.1, 9.1, 10.1). High energy is handled by HB Thm 1.6:
  Σ_{t1,t2}|Σ a_n n^{i(t1−t2)}|² ≲ |T|²N + |T|N² + |T|^{5/4}T^{1/2}N.
  That gives Prop 11.1: E(W) ≲ |W|N^{4−4σ} + |W|^{21/8}T^{1/4}N^{1−2σ} + |W|³N^{1−2σ}, whose first term "generally sets the limitations".
- **Characters.** CGL 2507.08296: GM extended to all χ mod q gives exponent 7/3, and 30/13 for T-smooth q. The gap is that subdivision in q needs divisors q₁ | q, and no Huxley-type reflection is known for GM (their §2).
- **Others.** TTY 2501.16779 (the exponent database ANTEDB, energy exponents LV*, A*, Heath-Brown relation Thm 63); 2403.13157 (density implies large values); 2610.18923 (trilinear sums via GM/CGL). Paperclip evidence API: not used (quota shared, arXiv sufficed).

## T1 — what a sticky/multiscale energy theorem could buy, at most (`t1_energy_lp.py`, `t1_out.txt`, `t1_table.md`)
Method: exact exponent bookkeeping of GM's chain |W|³N^{6σ−3} ≲ N³ + S2 + S3, using:
- S2 from Prop 6.1, for every k;
- S3 ≲ T²|W|^{3/2} + TN|W|^{1/2}E^{1/2} (Prop 10.1);
- E(W) from a hypothesis H, with E ≥ |W|² enforced;
- subdivision into intervals of length N^{t₀}, t₀ ∈ [1, 4/3];
- the zero detector with powers (GM §13), taking the worst case over N = T^ν, ν ≤ 1/2.

Sanity checks: the MVT alone reproduces Ingham 3/(2−σ); adding Huxley reproduces 3/(3σ−1); H0 reproduces GM Prop 3.1 ((18−20σ)/5 at τ = 6/5, exact) and 15/(3+5σ) to within 0.003.

| input (A(σ) from the chain; MVT and Huxley always available) | 0.70 | 0.72 | 0.75 | 0.78 | 0.80 | sup |
|---|---|---|---|---|---|---|
| H0: GM Prop 11.1 | 2.307 | 2.272 | 2.220 | 2.172 | 2.142 | 2.3067 |
| H1: HB third term dropped from E (R309 idea) | 2.307 | 2.272 | 2.220 | 2.172 | 2.142 | 2.3067 |
| H1 + Conj 7.3 (HB third term dropped in S2 too) | 2.307 | 2.272 | 2.220 | 2.172 | 2.142 | 2.3067 |
| GM bottleneck term cut by N^{−0.1} / N^{−0.2} | 2.307 | 2.262 | 2.204 / 2.197 | 2.156 | 2.142 | 2.3067 |
| H3: minimal possible energy, \|W\|² + \|W\|⁴/T | 2.307 | 2.262 | 2.197 | 2.156 | 2.142 | 2.3067 |
| H3 + Conj 7.3 (both inputs perfect) | 2.307 | 2.255 | 2.180 | 2.156 | 2.142 | 2.3067 |
| H2: "union of short APs" bound \|W\|²N^{2−2σ} alone | 2.307 | 2.318 | 2.250 | 2.184 | 2.142 | 2.3244 |

Verdicts:
1. **R309's suggestion is void.** It proposed saving on HB's (qT)^{1/2}|W|^{5/4}N term. For N ≥ T^{2/3} that term is already dominated (GM footnote 2). Removing it everywhere, which is Lindelöf-strength on average, changes no A(σ) on [0.60, 0.84]. CGL's limiting terms (their Remark 11.2) are the same two as GM's, so the q-version changes nothing either.
2. **30/13 cannot move, and the density hypothesis at σ = 3/4 stays out of reach.** Even with perfect energy and Conj 7.3, sup A = 30/13 and A(3/4) ≥ 2.180 > 2. GM's chain beats the MVT only if both of these hold:
   - τ ≥ 4 − 4σ, from the S2 diagonal term N²|W|², which gives |W| ≤ N^{5−6σ};
   - τ ≤ 6σ − 3, from the cubic random-W term T²|W|^{3/2}, which gives |W| ≤ T^{4/3}N^{2−4σ}.

   The window is nonempty only for σ > 7/10, with corner (7/10, 6/5) = GM's choice T = N^{6/5}. The energy term T^{1/2}N^{3−4σ} also meets the MVT at that corner, so removing it cannot open the window. Neither of the two binding terms depends on energy or on HB's third term. Both are limits of the cubic-trace method. N²|W|² is the n-diagonal |W|²M of HB's bound for the dual sum of length M ≲ T/N, and it is attained for generic W. T²|W|^{3/2} is the triangle inequality over spikes with |Ŵ| ≈ |W|^{1/2}.
3. **The "sticky" structure is already beaten.** Unions of short progressions have E = |W|²K with K ≤ N^{2−2σ} (Lemma 1.7), which is weaker than Prop 11.1 once |W| > N^{2−2σ} (row H2 is worse than H0). The configurations that saturate Prop 11.1's Hölder step are m generic translates of K-term progressions with K·m^{1/2} = N^{2−2σ}. These are non-sticky at the coarse scale.

## T2 — are the GM-extremal non-sticky configurations large-value sets? (`t2_clusters.py`, `t2_out.txt`)
Setup: W = m random blocks of K consecutive integers with K·m^{1/2} = N^{2−2σ}, σ = 3/4, and Gram matrix G = (Σ_n n^{i(t−t')}). For |b_n| ≤ 1 we have Σ_W|D|² ≤ λ₁(G)·N. So bestL2 = λ₁N/(|W|V²) < 1 proves that W is not a large-value set; this step is rigorous.

| (N, τ), T/(2πN) | m = 4 | m = 16 | m = 64 | trace loss (tr G³)^{1/3}/λ₁ at m = 64 |
|---|---|---|---|---|
| (2000, 1.2), 0.7 | 4.94 | 2.11–2.18 | **0.75** | 3.4 |
| (4000, 1.2), 0.8 | 5.14 | 2.26 | **0.92** | 3.6 |
| (1000, 1.6), 10 | 4.83–4.89 | 1.95–2.00 | **0.70** | 2.8 |
| (600, 1.8), 27 | 4.57–4.73 | 1.79–1.84 | **0.63** | 2.6 |

- λ₁(W)/N stays at single-block size (4.9–11, falling as K shrinks) and does not grow with m, because separated clusters are nearly orthogonal. So random-position block configurations obey |W| ≲ 10·N^{2−2σ}, which is Montgomery's scale.
- Unimodular b built from the top eigenvector makes only 3–18% of the points large at m = 64.
- The cubic trace loses about m^{1/3}, because it counts all m cluster eigenvalues. This is the slack in GM's chain in the non-sticky case.
- Turning near-orthogonality into a theorem for adversarial positions means bounding Σ n^{iu} on the difference set. That is the Halász–Montgomery problem, which says nothing for σ ≤ 3/4. Entrywise Halász–Montgomery already overestimates λ₁ by about 2× at m = 64 (HM/N = 11.5 against λ₁/N = 5.7 at τ = 1.6; 10.9 against 4.9 at τ = 1.8). Beating it needs operator-norm cancellation, which is exactly Guth's §8 barrier.
- **Verdict:** the dichotomy holds in these numerics, but its non-sticky half is the original problem. Caveat: N is small; at τ = 1.2 with N ≤ 4000 there are no dual terms (T < 2πN), which is why the τ = 1.6 and 1.8 runs were added.

## T3 — Kakeya structure does occur, but only below the barrier line (`t3_bush.py`, `t3_out.txt`)
Bourgain–Kakeya bush: D = Σ_I ε_I Σ_{n∈I} n^{it} with random signs ε_I and all translations 0. These are genuine Dirichlet polynomials with |b_n| = 1, and |D|² ≈ Σ|D_I|² ≤ J·L² = N²T^{−1/2}. So they reach V = N^σ if and only if τ ≤ 4 − 4σ. This is the same line as T1(2) and as Guth's §8.1 barrier σ = 1 − τ/4; I found no source that states this coincidence.

Numbers at N = 4000, τ = 1.2: |W| = 400, 208, 78, 44, 10, 0 at σ = 0.70, 0.72, 0.74, 0.76, 0.78, 0.80, against N^{2−2σ} = 145, 104, 75, 54, 38, 28. GM's bound is 1668, 888, 479, 262 over the same range.

So these examples sit at Montgomery's level. They are not obstructions to GM, and the claim that "Kakeya-type structure cannot occur" is false: it occurs exactly where the method fails. At (1000, 1.6), Gaussian tails dominate because N is small.

## T4 — characters: the prime-modulus gap has no scales to induct on (analytic, no test)
CGL's 7/3 for general q (against 30/13 for T-smooth q) comes from the missing subdivision in the q-aspect. A Wang–Zahl-style substitute, in which approximate structures (convex sets, Bohr sets) replace exact ones, would need cells on which χχ̄' has small analytic conductor. For prime q, χ ≠ χ' gives conductor q. And no non-principal χ is uniformly close to 1 on [1, N] once N ≥ q^{1/4+ε} (Burgess), while GM needs N ≥ (qT)^{2/3}. So prime-modulus families are single-scale, like finite-field Kakeya, where Dvir's polynomial method rather than induction on scales won. A Stepanov- or Dvir-type mechanism for character-sum large values is a speculative remark only; I found no concrete step.

## Ranking (RH relevance × chance × gain)
| transfer | node | chance | gain if it worked | verdict |
|---|---|---|---|---|
| T1 multiscale energy / HB third term | N7 | high that some saving is provable | A(3/4) −0.04 at most; 30/13 unchanged; density hypothesis not reached | dead |
| T2 sticky/non-sticky dichotomy for clusters | N7 | low: the non-sticky half is Halász–Montgomery or Guth's barrier | would open the window left of 7/10 | dead |
| T3 Bourgain–Kakeya / Wang–Zahl volume bounds | N7 | black box: no gain (Guth §9) | none known | dead |
| T4 Kakeya-type localization for prime q | N7/N11 | none (single scale) | 7/3 → 30/13 for prime q | dead; the polynomial-method remark is untested |

## Files
`Riemann/Deep/ExploreKakeya/`:
- `t1_energy_lp.py`, `t1_summarize.py`, `t1_table.py`
- `t1_out_raw.txt`, `t1_out.txt`, `t1_table.md`
- `t2_clusters.py`, `t2_out.txt`
- `t3_bush.py`, `t3_out.txt`

All runs were single-threaded under a perl alarm (longest 331 s). Scratch: `.scratch/explore4/STATE.md`, `survey_notes.md` (downloads deleted).

Suggested NOGO line: explore4 (Orca: Kakeya / Wang–Zahl / multiscale) NO IDEA. Energy and Heath-Brown inputs cannot move 30/13 (pinned at (7/10, 6/5) on Guth's barrier line σ = 1 − τ/4, where random-sign fat-progression polynomials first reach N^σ). Perfect inputs give A(3/4) ≥ 2.18. Prime q has no scales. [ExploreKakeya](Riemann/Deep/ExploreKakeya/REPORT.md)

RH relevance: N7 (zero density, large values); secondarily N11 (q-aspect) and N9 (Conj 7.3 is Lindelöf on average over difference sets).
