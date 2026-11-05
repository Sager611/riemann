SOUND: Re s > B_2 = (777 − √1081)/852 = 0.873381966613207519929… is zero-free for every finite-order Hecke L-function over every number field (so also for ζ and every Dirichlet L-function), conditional on Liu's arXiv 2610.26835v2 (allowed generators and reciprocity, theta and Poisson set-up, inverse moments, product moment M2 = Thm 12.16, pointwise prime bound Lemma 12.7, the §§12.5–12.7 transform identities, detector, comparison formulas of §§13–14, continuation, general-field transfer). No false step. Every finite inequality, both certificate tables and the endpoint were re-derived exactly. The new analytic steps (selected large sieve, two moment inductions, contour extension, principal box) check out step by step, but the inductions are written at sketch level on top of those cited identities.

# Liu, "Selected character moments, zero-free half-planes for Hecke L-functions over number fields, and Quasi-Riemann Hypothesis", Hexagon 2611.00082v1: audit (node N7, 2026-11-05)

**Short answer.** The boundary is right and the proof hangs together. B_2 itself is not new. It is the limit point of a comparison scheme for Dirichlet L-functions, identified on 31 Oct 2026 by Joseph Bebel and Devin Watson on the public QuasiRiemannTracker, a GitHub leaderboard that only accepts Lean-checked proofs. There, Bebel's Lean proof of Re s > B_2 + 2.2·10⁻²⁵⁴ for every Dirichlet L-function, ζ included, passed the tracker's kernel check on 1 Nov. Liu's new contribution is the transfer of that scheme to ideal characters over every number field. That covers Dedekind zeta functions and all finite-order Hecke L-functions.

Source: TeX zip from hexagonmath.org (only v1 exists, 4 Nov 2026; arXiv 2611.00082 is an unrelated paper). "l." = line of number_field_new_v2-1.tex (4817 lines). [NF3] = Liu, arXiv 2610.26835v2 (number_fields_v4.tex, 13294 lines; "NF3 l."). [V14] = Liu, Hexagon 2610.00655v3 (B_tw, audited in ../Hexagon/).

## 1. What is new, and the chain to B_2
| step | where | content | relative to NF3 / V14 |
|---|---|---|---|
| selected large sieve | Thm op:selected l.1066; Lemmas op:inheritance l.518, op:amplification l.537, op:initial l.687, op:reflection l.728, Prop op:squarefree l.816, Lemmas op:phases l.938, op:natural-near l.972 | for sextic characters χ_k restricted to parameters where a fixed product of short prime sums is large (slope 1/2 < K ≤ 3/4): Σ_{k∈E} \|Σ_{Nn≤X} a_n χ_k(n)\|² ≪ Z^ε (X + H^{1−K} X^K) ‖a‖² | new over number fields; Dirichlet/Eisenstein analogue is Bebel's Lean-checked PR93 |
| sixth-power restoration | Thm op:restoration l.1202 | high-amplitude terminal max{(1−K)m + KA, (m+q)/6 + 5A/6 (+ pz/6 for a block)} | new (PR97 analogue) |
| first induction | Prop mom:fixed-k l.1707, Cor mom:three l.1974 | moment for n_1+n_2+6kz ≤ M with k ↓ 1/2, hence n_1+n_2+3z ≤ M | new; NF3's M2 had 6κ ≈ 4.5 |
| Gauss terminal | Lemma mom:gauss-mean l.2031, Prop mom:gauss-terminal l.2176 | completion bound for Gauss polynomials at tiny free width | new |
| second induction | Lemma mom:marked-low l.2429, Prop mom:face-step l.2598, Thm mom:two l.2834 | faces G(1−δ, L_n), L_n ↓ L_* = 6Bk/(1+k) → 2, so n_1+n_2+2z ≤ M | new (PR97/99 analogue) |
| contours | Lemmas profiles-lemma l.2956, contours l.3027 | numerator contour shift θ up to δ/2 (NF3: θ ≤ 1/10) | new proof; replaces V14 Lemmas 5.1–5.2 |
| comparison | §6, l.3129–4265 | count R_d = 1 − δF_d from the new fourth moment; numerator gain γ_2; active margin g = G_N − E(R_d), in which the old coefficient p_1 cancels | the same scalar system as tracker PRs 101–103 (same contact δ* = (61−√1081)/66, x = 1/2, h = 0.8096863494) |
| principal box | Lemma principal-box l.4298 | principal Euler product holomorphic and nonzero for Re s ≥ 0.87 (NF3: 437/500) | new corner, same local polynomial as NF3 l.8796–8852 |
| Mellin, descent | Prop fixed-field-result l.4433, proof of Thm main l.4505 | as NF3 Lemma 13.11 and Prop 19.1 | unchanged |

V14 is off the chain. Its Parts I–III (B_ex, B_tw machinery) are not used. Only its contour lemma is cited, as a template, and that lemma is re-proved here.

## 2. Findings
1. **Endpoint from first principles (A1–A12).** g is affine in b. At x = 1/2 its b-coefficient is −(33δ²−61δ+20)/(24(5−4δ)), which vanishes at δ_2 = (61−√1081)/66. At δ_2, g = 5/12 − δ_2/2 − (3/4+δ_2)ℓ for every b. So g = 0 fixes ℓ_2, and B = (11/3 − ℓ_2)/4 is exactly (777−√1081)/852, the smaller root of 1278B² − 2331B + 1061. The double root in δ forces b_2 = −11/71 + 595√1081/76751. The factorization (sc:contact-factorization) and the formula for ∂_B g are exact. For B < B_2, g(B,b,δ_2,1/2) = (3+4δ_2)(B−B_2) < 0 for every b, so B_2 is the exact infimum of this system. T_103 − B_2 = 2.162·10⁻²⁵⁴ by an exact sign test.
2. **Certificates reproduced identically (C).** We rebuilt each multiplied function as an exact polynomial over Q or Q(√1081): interpolation, confirmed at random points. Then we ran Liu's Bernstein bisection rule. All 25 entries at B_rat = 0.8733821, all 24 at (B_2, b_2) and the 5 contact entries (App. A, l.4621–4669) give the same depth, the same number of rectangles and the same printed k. Example: the tight entry G_N − E(R_d) − 10⁻⁸ at B_rat gives (13, 14, 194222). The denominator bounds Δ ≥ 12/5 and J ≥ 4373305/6747966 hold, and the four rectangles cover [1/8,5/6]×[0,1/2]. All 23 constant entries (l.4678–4703) reproduce to the printed 12 digits at both parameter sets.
3. **Scalar derivations (B, E).** All of these are exact:
   - the detector balance r_d = 2m_d (so F_I = F_2 identically) and R_d = 1/6 + A·t_d;
   - the crossing in the old count R_s, the numerator gains, A_±′, G_N, G_L, G_R and the reference identities;
   - the optimality identities: (sc:free-short-exclusion) with bound 41671/12000000, (sc:free-high-plane), the positive dual uL_0 + vL_K + wH = B_2 identically, and the negative dual with 13/D_− = 0.894938658758…;
   - the unequal-factor neutrality and the saturation model.
4. **Inputs match NF3.** We parsed NF3's TeX numbering: all 29 cited items exist, with matching titles.
   - Thm 12.16 (NF3 l.7575) is the M2 used here. For z = 0 it allows arbitrary lengths, with the principal character excluded. For z > 0 it needs n_1+n_2+6κz ≤ M, 13/18 ≤ κ ≤ 1 and β_* ≤ (1+κ)/2.
   - Lemma 12.7 is the pointwise prime bound.
   - Lemma 13.10 gives the region Re s ≥ a+16e, Re w ≥ 1−a−6e, Re z ≥ 17/50 for 1/2 ≤ a ≤ 11/12. This paper needs a ∈ [61/100, 33/40].
   - NF3's refined endpoint (l.10091) is E(R) with c_κ = 2/p_1 and R_κ = R_s, and its slope (3R+1)/4 + θ/2 − 17/50 is S(R).
   - κ_1 = 0.749774 is at least the 0.74977 that our Kappa audit requires, and B_♭ < 0.87488664704 < (1+κ_1)/2 (Sturm table, App. B; D9).
5. **The 0.874 floor of our previous audit is removed legitimately.** Every use of 437/500 in NF3 concerns the principal Euler product region (l.693, 1248, 8796–8852, 9644–9699, 9896, 10116, 10711, 10832). Liu recomputes the same cleared local polynomial, which is identical to NF3 l.8805–8809, at the corner (0.87, 19/20, 33/200). There the 12 monomial decay exponents are ≥ 1.81 > 1 and the 5 selected-factor error exponents are ≥ 0.87 (D7). The factorization at w = 1, z = 1/6 is exact (D8). Any B > 0.87 is allowed, and B_2 clears it by 0.0034.
6. **Selected large sieve (§3), step by step (F2, F3).**
   - Initial bound X + H²: Gram matrix with the primitive functional equation. The residue at s = 0 vanishes because the L-function vanishes there to order ≥ r_F.
   - Reflection: the scales XD/(RH), Y_0 = H²R/(XD²) and the budget H^{q+2K−1}X^{1−K}D^{1−q−2K}R^{K−1} are exact.
   - Recursion: q ↦ (1−K)(q+2K−1)/K contracts to 1−K, with L = H^{a/K} ≥ H.
   - Amplification uses only multiplicativity, a divisor bound and the selection itself.
   - The 25-cell table is the digit (i−j) mod 6, and the natural-near-range bounds are as printed.
   - Why the reflection closes: for coprime residual parts the Gauss-sum cross factor χ_a(b)·conj(χ_b(a)) is 1 by reciprocity, so the root number splits into a function of a times a function of b (Lemma op:phases; the finite packets come from NF3 Lemmas 3.4 and 12.9).
   - Selection is what breaks the circularity that keeps the unselected cubic large-sieve conjecture open: amplification gains Z^{−Kw} at short source lengths.
7. **Restoration and the first induction (F1, F4, F5, F8, F9).**
   - Both suprema over the sixth-root split (whole product and block) were re-derived. A random search over 4000×401 samples found no excess.
   - The centering gap γ_k = (2k−1)/(2(6k−1)) is positive as soon as 6k−1 > 2 (NF3 needed 6κ−1 ≥ 10/3).
   - After one return the defect is ≤ 6(w_a+ℓ_a), while the slack is Δ′ ≥ w_a+ℓ_a. So a low block of width D′/(6k) costs exactly D′/6 ≤ Δ′, whatever k is. High cells end with the selected large sieve at slope K ∈ [k, κ_1].
   - Amplitude cells are formed only on returned sides, never passed through Poisson summation.
   - Exact checks: (mom:second-M), (mom:full-allowance), the raw and six-fold defect bounds, the pool errors (σ/3)(−2+3s−si+(t−1)e) ≤ 0, and the local tables for all valuations ≤ 60 (first-brace, second-sign, F_2 ≥ 2c̄_2/3).
8. **Second induction (F6–F8).**
   - Marked identities, all exact: the M_* expansion, E_0 = exterior, (mom:marked-defect), the coefficient table, the high-K and high-R planes and the origin payment. The gate L_0 ≤ 6B_0k makes every coefficient ≤ 0.
   - Face-step constants: 50/291, 41/291, 41/582, the high gates 0.485 and 4.0675, and Φ + 2ℓ/3 ≤ 0.
   - The curve: ρ_n = (a/2)(1+k−6Bk/L_n) exactly, L_n ↓ L_* ∈ (2,3), 6Bk−2(1+k) = δ−6δ², and the final width removal.
9. **Counter-instance search on the Gauss-mean lemma (G).** We tested it over Z with Jacobi symbols for all odd moduli (prime powers included, natural zeros kept).
   - The exact facts it uses hold for u, v < 46: F(u,v;0) = 1_{u=v}φ(u) and |F(u,v;j)| ≤ gcd·1_{gcd | j}.
   - The bound Σ_h |D^{−1/2}Σ_u B(u)G(u,h)|² ≪ (H/D)Σ|B|² + (Σ|B|)² holds with ratio ≤ 0.57 for D = 40…120 and H = D/4…D²/4, including coefficients aligned to one frequency.
   - No counter-instance.
10. **Contours and descent.** Lemma contours stays inside NF3 Lemma 13.10's region: θ ≤ δ/2 = a − 1/2. The AFE weights are bounded on Re ω ∈ [1/7, 6/7] by Stirling. The descent through K(ζ_12) is NF3 Prop 19.1, audited sound in ../LiuHecke/.
11. **Weak points, not flaws.**
    - (a) §4 imports from NF3 §§12.5–12.7, without re-deriving them, the first and second transform widths, the allowances (B_C, Λ_C, λ_0), the pool table, the local Gauss tables and the centred exceptional bound (l.2687). We checked all the algebra built on them, not those formulas themselves.
    - (b) Both finite inductions are argued in prose, by rank and margin allocation. Their Dirichlet/Eisenstein analogues are Lean-checked (PR93/97/99); the number-field versions are not.
    - (c) The tracker's Lean faces carry a width cap c < 1/6 ("s ≤ c₀H", Sil's StepGates.lean). Liu's Thm mom:two has no cap (l.2817), and his numerator band uses prime widths up to 0.20. We found no step of his proof that needs a cap, and a cap would not matter anyway. With widths capped at 3/20 on [0.30, 0.35], the margin becomes g + δx(3h/40 − T/4) > g (F10). Falling back to the baseline there instead would fail (grid minimum −4.9·10⁻⁴).
    - (d) The paper is unrefereed, single-author, and credits AI assistance (l.218).

## 3. Inputs and their status
| # | input as used | source | status |
|---|---|---|---|
| I1 | allowed generators, reciprocity, theta, radial functional equation, orbit Poisson, coprimality expansion | NF3 Lemma 3.4, Props 4.1, 5.2, 8.1, Lemma 12.5 | unrefereed; field step audited (LiuHecke) |
| I2 | M2 (Thm 12.16), including z = 0 at arbitrary lengths; pointwise prime bound (Lemma 12.7); inverse moments (Thms 10.7, 11.2) | NF3 | unrefereed; κ < 3/4 audited (Kappa) |
| I3 | §§12.5–12.7 transform identities; Lemmas 12.6, 12.9, 12.10 | NF3 | unrefereed; used in §4 without re-derivation |
| I4 | detector and Sobolev, direct estimate, comparison formulas E, S, N, continuation | NF3 Lemmas 13.7, 13.9–13.11, Props 13.4, 14.3, 14.6, Lemmas 14.7–14.8, (13.13), Thm 13.15 | unrefereed; the same formulas were audited in ../Hexagon/ |
| I5 | seed β_*(F) ≤ B_♭ | NF3 Thm 17.1 | unrefereed; arithmetic re-checked (D9) |
| I6 | general-field transfer and descent | NF3 Thm 18.2, Prop 19.1; Milne CFT | audited (LiuHecke) / standard |

## 4. Novelty and priority
- **B_2 is not Liu's.** It is the contact of the Dirichlet comparison on QuasiRiemannTracker (github.com/GTimG/QuasiRiemannTracker), a leaderboard that merges a bound only after its CI kernel-checks the Lean proof.
  - Bebel's PR101 (opened 31 Oct 2026, 06:03 UTC) and Watson's PR102 (31 Oct, 16:57 UTC) name (777−√1081)/852. They give the contact point ((61−√1081)/66, 1/2) and h = 0.8096863494, exactly Liu's δ_2 and h. Both were merged on 1 Nov.
  - Bebel's PR103 (1 Nov) proves L(s,χ) ≠ 0 for Re s > T_103 = B_2 + 2.162·10⁻²⁵⁴ for every Dirichlet character, ζ included. Its admission record says PASS (2026-11-01 11:03 UTC). Candidate.bound uses only propext, Classical.choice and Quot.sound.
  - Liu cites PRs 93, 97, 99, 101 and 103 and states that B_2 "already occurs in the Dirichlet comparison" (l.179–188, 3147).
- **What is Liu's.** The same boundary for every finite-order Hecke L-function over every number field, Dedekind zeta functions included. Also the open half-plane at B_2 itself; for ζ that gains only 2·10⁻²⁵⁴ over PR103.
- **Elsewhere.** The arXiv math.NT listings for Oct 2026 (854 entries) and 1–5 Nov 2026 (140 entries) contain no B_2. The lowest arXiv claims are Chang–Wang 2611.00924 (0.874462) and Liu 2610.26835 (B_♭). Hexagon has only v1 of 2611.00082. Our repo never mentioned the tracker.

## 5. N7 context
- TREE N7 lists "no zero with Re s > 7/8 (OpenAI; Lean port)" as the machine-checked record. That is out of date. Since 1 Nov the tracker has had a kernel-checked Re s > T_103 = 0.87338196661320751992902677864… for ζ and all Dirichlet L-functions. Merged Dirichlet bounds there were already below V14's B_tw from 27 Oct (0.873985, PR84), then fell through 0.8739… (28 Oct) and 0.87341198 (PR97, 30 Oct) to that bound.
- For Hecke L-functions over number fields, B_2 (this paper, sound given NF3) supersedes V14's B_tw = 0.874058504 and Chang–Wang's 0.874462.
- The method's floor inside this comparison system is B_2 itself (Thm sc:optimality, checked: B1–B5). The unequal-factor model (l.4167–4265) shows that positive interpolation of these moments gives nothing more. Going lower needs a new input: a phase-sensitive mixed estimate, or a better detector or reference.

## 6. Reproduction
Run python3 verify.py (stdlib + mpmath 1.3.0, Python 3.14): 99 checks, 0 failed, about 20 s, 30 MB, one thread. The Bernstein certificates and the algebraic identities use exact rational or Q(√1081) arithmetic, so a PASS there is a proof. The one assumption is that each multiplied function is a polynomial of degree ≤ 6 in each variable. That holds by construction, and the interpolation is confirmed at extra random points. The checks of affine identities at random points, the restoration search and the Gauss-mean test are evidence, not proofs.

Sources: [Hexagon 2611.00082v1](https://hexagonmath.org/abs/2611.00082) · [arXiv 2610.26835v2](https://arxiv.org/abs/2610.26835v2) · [QuasiRiemannTracker PR101](https://github.com/GTimG/QuasiRiemannTracker/pull/101), [PR102](https://github.com/GTimG/QuasiRiemannTracker/pull/102), [PR103](https://github.com/GTimG/QuasiRiemannTracker/pull/103) · [arXiv math.NT 2026-11](https://arxiv.org/list/math.NT/2026-11) · [previous audit](../Hexagon/REPORT.md) · [LiuHecke](../LiuHecke/REPORT.md) · [Kappa](../Kappa/REPORT.md)
