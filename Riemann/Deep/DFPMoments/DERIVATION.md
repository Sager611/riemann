# DERIVATION: optimising de Faveri–Pandey's moment savings (arXiv 2610.01035 v1)

Notation. Exponents of T throughout. s = log_T V, ν = log_T N, ω = log_T W. N(T,V): largest 1-spaced subset of
{t ∈ [T,2T] : |ζ(1/2+it)| ≥ V}. A bound "E(s)" means N(T,T^s) ≺ T^{E(s)}; μ(A) means M_A(T) ≺ T^{μ(A)}; ≺ hides T^ε.
"script" = `dfp_moments.py --all`, log `logs/all.log` (982 PASS, 0 FAIL). TeX = `zeta_paper_writeup.tex` of 2610.01035.

## 0. From large values to moments
D1. μ(A) ≤ max(1, sup_{0≤s≤13/84} (A s + E*(s))), E* = pointwise min of valid bounds.
Proof: a maximal 1-spaced subset S of {t ∈ [T,2T] : |ζ| ≥ V} has every point of the set within 1, so the set has
measure ≤ 2|S| ≤ 2N(T,V). Split [T,2T] by |ζ| ∈ (V/2, V], V dyadic, V ≤ T^{13/84+ε} (Bourgain): ∫_T^{2T}|ζ|^A ≤ T +
Σ_V V^A·2N(T,V/2) ≺ T + sup_V V^A N(T,V/2); sum over dyadic T. Same reduction as TeX proofs of Thms 1.1, 9.1 and TY24 §5.1.
For the true exponents the sup formula is an equality (ANTEDB Lemma mad), so no better moment bound follows from the same
large-value bounds; Hölder between moments is implied.
D2. Continuity: E_k(σ_-(k)) = 1 − 4σ_-(k) and E_k(σ_+(k)) = 2 − 12σ_+(k) (script part 1, all k ≤ 12; identity in §4 for all k),
so DFP's bounds join the global lines 1 − 4s, 2 − 12s and E* is continuous. A s + E*(s) is piecewise linear; its sup is at a
breakpoint: interval ends or pairwise line crossings, all evaluated exactly (script `vertices`, `sup_moment`).

## 1. Inputs (every bound used is published)
I1 4th moment: 1 − 4s. I2 Heath-Brown 12th: 2 − 12s. I3 Heath-Brown auxiliary (ANTEDB hb-12-aux): max(1 − 6s, 5 − 32s).
I4 Ivić Thm 8.2 (= TY24 Lemma 5.1): for an exponent pair (κ,λ), κ > 0: max(1 − 6s, 1 + λ/κ − 2(1+2κ+2λ)s/κ).
   Pairs: Bourgain (13/84,55/84); TY24 Lemma 1.1 (4); TTY25 Thm "new exponent pairs" (4); Huxley–Watt, Huxley–Kolesnik,
   Huxley Ch. 17, Sargos (6); Huxley Table 17.3, m ≤ 11; derivative tests (Robert, Sargos; 6); Heath-Brown 2017
   (p_m,q_m) = (2/((m−1)²(m+2)), 1 − (3m−2)/(m(m−1)(m+2))), 3 ≤ m ≤ 120; (1/2,1/2). Closed under the A, B and Sargos C
   processes to depth 4: 13761 pairs ("+ε" in pairs costs T^ε only).
I5 DFP Thm 8.1 as printed, k = 2..60, on [σ_-(k), σ_+(k)] (k = 2 is Thm 1.3: (24 − 128s)/17 on [7/60, 5/38]).
I6 Bourgain μ(1/2) ≤ 13/84 (range of s).
Known pairs = TY24 Thm reg_thm_1 (hull H, closed under A, B, C) plus TTY25's four. At fixed s the I4 exponent is
1 − 4s + (λ(1−4s) − 2s)/κ, linear-fractional in (κ,λ), so its minimum over a convex polygon of pairs is at a vertex; all
vertices of H are in the list except (p_m,q_m), m > 120, handled by the tail bound in §4.

## 2. DFP's final system, re-derived (TeX §9, proof of Thm 8.1)
After pigeonholing F, N, Z and G = V² (log T)^{-100}, with F = 1, Z = (TN)^{1/4} (every term decreases in F and Z),
R ≺ T^{max(T1,T2,T3,T4)} with, for N = T^ν, W = T^ω (TeX display after "Multiplying by Z/(TN)^{1/4}"):
- T1 = 1/2 + 3ν/2 − 2s − ω  (N²/Z² term of Thm 5.1)
- T2 = 1/(2k) + (2 − 3/(2k))ν − 2s − (k−1)ω/(2k)  (middle term of Thm 5.1, E_2 ≍ R)
- T3 = 3/4 − k/2 + (k/2+1/4)ν + (k−1)ω/2 − 2s + max(0, ν+ω−1)/4  (third term, Υ^{1/4}, Υ = 1 + NW/T)
- T4 = (6 + (16k−12)ν − 16ks − 4(k−1)ω − 12s)/(8k+1)  (Iwaniec's E_2 term T√R/V^{3/2}, Cor 6.2, cleared from R)
At ν = 1 − 4s these are TeX (eq:curyjytduj): 2 − 8s − τ, 2 − 1/k − (10 − 6/k)s − (k−1)τ/(2k), 1 − (2k+4)s + (2k−1)τ/4,
(16k − 6 − (80k−36)s − 4(k−1)τ)/(8k+1).
Constraints (TeX eq:proof_ranges and the checks in the proof): C0 ν ≥ 1/3 (Atkinson); ν ≤ 1 − 4s (N ≪ T/G²); C1 ν ≤ 3/5;
C2 ω ≥ ν; C3 ω ≤ 4ν/3 (Thm 5.1); C4 ν/2 + 2ω ≤ 3/2 (Cor 7.2, phase truncation); C5 ω ≤ (1+ν)/2;
C6 1 + ω − 3ν > 0 strictly (|η_j| ≤ c_k W/N²); C7a 4s ≥ (1−ν)/2, C7b 4s ≥ ω + ν/2 − 1/2 (absorbing the middle E_2 term).
DFP take ω = τ_k(s)·ν/(1−4s). Script part 1 (811 checks) reproduces, for k = 2..12 at both range ends and the midpoint:
T4 = E_k; T3 = E_k (k ≥ 3), T2 = T3 = 11/8 − 29s/4 (k = 2); T2 ≤ E_k; T1 < E_k; every constraint with slack, at N = T^{1/3}
and N = T^{1−4s}; E_k(σ_±) identities; Thm 1.3 (24/17, 128/17); Thm 1.1 (max(1 + 2·7/60, 24/17 − 26·7/(17·60),
2 − 6·5/38) = 37/30); Thm 9.1 (E_6 = (1070 − 6164s)/619 on [451/3688, 21/158]; B(p_10,q_10) = (64/135,122/243);
Ivić line 593/288 − 3587s/288; σ_2 = 58907/445121; 2 − 520/445121).

## 3. Free parameters of the chain: none gains
(i) W per N. max(T1..T4) is jointly convex piecewise linear in (ν,ω) on a convex polygon; for fixed ν its minimum over ω
is at an endpoint, the kink ω = 1 − ν, or a crossing of two pieces (script `g6`, exact). (ii) N. That minimum is convex
in ν, so the worst N is an endpoint; N = T^{1/3} is never worse than N = T^{1−4s}. Result at every curve vertex and
every k there (script part 6 (a), 53 checks): the optimum equals E_k(s). DFP's τ is optimal.
(iii) k. All k: §4 tail. Letting k (and W) depend on (N, Z, F) changes nothing: the single case N = T^{1−4s}, F = 1,
Z = (TN)^{1/4} cannot be excluded, and there the best choice is min_k E_k(s) by (i).
(iv) G. G ≤ V² is forced by the case F = 1; every term decreases in G, so G = V².
(v) Υ-cap. Thm 5.1's middle term carries sup_{1≤D≪Υ} E_2(D)/√D. With E_2(D) ≤ min((2D+1)R_1, R_1(1+X')) this sup is
R_1·√min(1+X', Υ) (maximum at D = 1+X' if ≤ Υ, else at D = Υ); DFP used √(1+X'). Rearranging gives the extra option
T2 + max(0, ν+ω−1)/(4k) for the middle term: R ≺ T^{max(T1, T3, min(max(T2,T4), T2 + max(0,ν+ω−1)/(4k)))}.
No gain: part 6 (b) (53 vertex checks) and (b') (11 k × 25 points, exact, max gain 0).
(vi) Ranges. The chain with optimal W outside [σ_-(k), σ_+(k)], k ≤ 20, never beats E* (part 6 (c): 120 exact points on
[3/26, 13/84], max gain 0).
So the printed E_k are the chain's optimum; the moment gains below come only from combining all of I1–I5 optimally.

## 4. Moments 4 < A ≤ 12: curve, table, optimality
μ*(A) = max(1, sup_s (A s + E*(s))), E* = min(I1–I5). Exact result (script CURVE):

| A | μ*(A) | maximising s | bounds meeting there |
|---|---|---|---|
| [4, 128/17] | 7A/60 + 8/15 | 7/60 | I1 ∩ E_2 |
| [128/17, 1340/157] | 173A/1342 + 296/671 | 173/1342 | E_2 ∩ E_3 |
| [1340/157, 2564/279] | 235A/1793 + 758/1793 | 235/1793 | E_3 ∩ E_4 |
| [2564/279, 4172/433] | 443A/3361 + 1398/3361 | 443/3361 | E_4 ∩ E_5 |
| [4172/433, 6164/619] | 715A/5409 + 2230/5409 | 715/5409 | E_5 ∩ E_6 |
| [6164/619, 3587/288] | 58907A/445121 + 182838/445121 | 58907/445121 | E_6 ∩ I4 at B(p_10,q_10) |

Breakpoints = slopes −E_k' = 128/17, 1340/157, 2564/279, 4172/433, 6164/619, and 3587/288 (I4 slope).
Before = DFP's printed results with Hölder: μ convex in A, μ(4) = 1, μ(6) ≤ 37/30, μ(12) ≤ 2 − 520/445121, lower convex
envelope. Classical (no DFP) = 1 + (A−4)/8. δ = 1 + (A−4)/8 − μ.

| A | μ* | δ* | δ before |
|---|---|---|---|
| 5 | 67/60 | 1/120 | 1/120 |
| 6 | 37/30 | 1/60 | 1/60 |
| 7 | 27/20 | 1/40 | 0.0140836 |
| 8 | 988/671 | 37/1342 = 0.0275708 | 0.0115005 |
| 9 | 2873/1793 | 0.0226576 | 0.0089174 |
| 10 | 771908/445121 | 0.0158468 | 0.0063344 |
| 11 | 830815/445121 | 0.0085075 | 0.0037513 |
| 12 | 889722/445121 = 2 − 520/445121 | 520/445121 | 520/445121 |

A ≤ 6: μ* = before (the first piece is Hölder between 4 and 6). 6 < A < 12: μ* < before strictly (script, 11 values).
before − μ* is concave on [6,12] (before linear, μ* convex), so it peaks at a breakpoint of μ*: at A = 128/17, μ* = 24/17,
before 1.4284604, gain 0.0166957 (script). At every breakpoint A = b_k, μ*(b_k) = a_k where E_k(s) = a_k − b_k s
(24/17, 242/157, 454/279, 730/433, 1070/619; script "breakpoint" lines).
Optimality over the whole input set (script `vertex_certificate`, exact): at each maximising vertex all 13823 bounds are
≥ E*(s). Tails: 2 − 12s − E_k(s) = 8((4k−3) − (30k−22)s)/(16k²+10k−17) ≤ (4 − 30s)/(2k) for s ≤ 3/22, below each vertex gap
2 − 12s − E*(s) for k > 61. For m > 120, (3m−2)/(m(m−1)(m+2)) ≤ 3/m² (⇔ 5m² ≥ 6m) gives κ_m ≥ 1/2 − 3/m², and λ_m ≥ 1/2;
for s > 1/8 the I4 exponent is then ≥ 2 − 12s − 6(8s−1)/(m²−6), deficit ≤ 2.41e-5 below every gap; for s ≤ 1/8 every pair
(λ ≥ 1/2) gives ≥ 1 − 4s ≥ E*. Hence μ* is the best I1–I6 give.

## 5. Above 12 (the range ν)
For 12 ≤ A < 3587/288, μ*(A) = (58907A + 182838)/445121 < classical μ_cl(A); at 3587/288 both equal 593/288, beyond they
coincide. So ν = 131/288 ≈ 0.4549 (paper: "a small ν"). This line is Hölder between μ(12) ≤ 2 − 520/445121 and
μ(3587/288) ≤ 593/288 (I4 at B(p_10,q_10): 2 − 520/445121 + (131/288)·58907/445121 = 593/288, script exact): explicit
range, not a new bound.

## 6. What blocks; what a better input gives
Sixth. μ*(6) = 1 + 2·7/60: E_2 is T4 with W at the C4 cap (ω = 1/2 + s), meeting 1 − 4s at s = 7/60. If Iwaniec's
short fourth moment (TeX Prop 6.1, term T√(J/L)) gained T^{−η}: T4 loses 2η/(8k+1), E_2 = (24 − 128s − 2η)/17, which
meets 1 − 4s at s* = (7 − 2η)/60; E_2^η − (11/8 − 29s/4) = (5 − 38s − 16η)/136, so T2 = T3 ≤ T4 iff s ≤ (5 − 16η)/38,
true at s* iff η ≤ 1/26 (script (d) algebra). Then 6s + E* ≤ 1 + 2s for s ≤ s* (I1); on [s*, (5−16η)/38] it is
6s + E_2^η, slope 6 − 128/17 < 0; on [(5−16η)/38, 5/38] it is ≤ 6s + 11/8 − 29s/4, decreasing (C3 and C6 hold for
3/26 < s ≤ 5/38); for s ≥ 5/38, ≤ 2 − 6s ≤ 23/19 < 16/13. So sixth = 1 + 2s* = 37/30 − η/15 (exact at η = 1/100, 1/50, 1/30).
Floor 16/13 (η = 1/26, s = 3/26): there C4 (ω ≤ 1/2 + s) meets C6 (ω > 2 − 12s), both 8/13, so for s ≤ 3/26 no W fits at
N = T^{1−4s}; also T2 = T3 = 1 − 4s = 7/13. Below 16/13 needs Thm 5.1 for |η| ≫ W/N² or the next Taylor term of
Atkinson's phase in Cor 7.2.
Twelfth. Fixed by E_6 ∩ I4(B(p_10,q_10)); at σ_2 the k = 6 terms are T2 = 0.410736, T3 = T4 = 0.410760. Removing T4
entirely moves it only to 1.9988280 (δ = 2 − 1.9988280, vs 520/445121; script (d), η = 1/100, 1/26, 1: T2 = T3 bind).
Blocked by Thm 5.1 itself (T2: middle term; T3: the diagonal (ND/V²)R^{2k−1} of Prop 4.1 with Υ^{1/4}) and by the
exponent pair: only a stronger Thm 5.1 or exponent pairs beyond Heath-Brown's (p_m,q_m) near (0,1) (equivalently
B-images near (1/2,1/2)) move it.

## 7. Other inputs: no gain
Dirichlet polynomials via the approximate functional equation: a point with |ζ| ≥ V gives a piece of length N = T^ν,
ν ≤ 1/2, with |Σ n^{−it}| ≥ N^{1/2}V = N^{σ'}, σ' = 1/2 + s/ν. Guth–Maynard's term T N^{12/5 − 4σ'} = T^{1 + 2ν/5 − 4s}
and Huxley's T N^{4 − 6σ'} = T^{1 + ν − 6s} exceed T^{1−4s} at ν = 1/2 (s = 1/8: 7/10 and 3/4 vs 1/2; script (e));
the bound must cover that piece. Same as DFP's remark (TeX §1.2). TY24's M(12+δ) ≤ 2 + δ/8 + Cδ^{3/2} is implied by I4.

## Not re-derived
DFP Thms 5.1, 8.1, 1.3 and Prop 4.1 are used as printed (v1, unrefereed); §2 re-derives only their final exponent
bookkeeping from the TeX, confirmed by their printed identities at every checked point.
