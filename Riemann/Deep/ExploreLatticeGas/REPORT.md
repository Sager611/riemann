IDEA: a screened sector lemma — Chasse's bound behind the N6 record (16T² − 4)/9 treats an off-line zero as isolated, but its real neighbours pull it to the axis (de Bruijn–Newman attraction, S = Σ_k 1/(γ − γ_k)²), which in toys delays the first non-hyperbolic Jensen degree by the predicted factor 2Sδ²/log(1 + 2Sδ²) (+5.5 %, +6.8 %, +9.9 % at zeros 5, 7, 10) and at T = 3·10¹² would raise the record by a factor of about 1.6 in the worst case the zero count allows, and about 6 for a typical zero (node N6).

# Exploration explore15: lattice gases, the hard-core model, the local lemma and zero-free partition functions → RH tree (2026-11-06)

Field: hard-core lattice gases and the independence polynomial; the Lovász local lemma and Shearer's bound (Scott–Sokal); Barvinok interpolation; Weitz's tree recursion and correlation decay; Peters–Regts; spectral independence; contraction proofs of real roots (Heilmann–Lieb, Chudnovsky–Seymour).
Read first and not repeated: explore1 ([ExploreStatMech](../ExploreStatMech/REPORT.md): Lee–Yang lifts, cluster expansions, Basu–Midha), R34, R48, R56, R358, R381 (matching polynomials, Heilmann–Lieb limits).
Seven scripts, each run as one single-threaded process at nice 19 under a hard alarm. Raw numbers are in `outputs.txt`.

## Short answer
1. **The field's own tools give nothing for RH.**
   - **Weitz-type recursion (A, N6): dies.** Jensen polynomials have an exact tree recursion J^{d+1,n} = J^{d,n} + X·J^{d,n+1}. It gives a true one-step lemma: if J^{d,n} and J^{d,n+1} interlace for all n, then all J^{d+1,n} are real-rooted. But interlacing is not inherited: counterexamples exist for d = 1…5. In a toy, J^{248,0} and J^{248,1} are real-rooted and J^{249,0} is not. Coefficient slack would have to beat 10^{−0.60d}, and no exact Heilmann–Lieb-type recursion exists.
   - **Local lemma over primes and the Mertens gas (C, N7, N10, N12): die at once.** The local lemma over primes reaches σ_s(x) = 1.052 … 1.039 → 1. Over F_q(t) it reaches only |u| < 1/q, short of RH's |u| = q^{−1/2}. M(n) is the hard-core gas of divisor chains at fugacity −1, which lies 0.61n times beyond Shearer's point, and the field's topology bounds |M(n)| only by 2n/π².
2. **The brief's adversarial test of the known wedge produced the one idea (B, N6).**
   - Taylor-coefficient criteria (Barvinok interpolation in disguise) see an off-line zero at height γ and offset δ only from order ∝ γ². Jensen sees it at degree ≈ γ²/(4δ²), Chasse's sector bound D_sec; Li at n ≈ (γ²/δ)·log n.
   - The bound is sharp for the lowest zero (first failure at D_sec + 2) but not above it. The excess is 2.2, 3.0, 4.7, 5.5, 6.8, 9.9 % at zeros 2, 3, 4, 5, 7, 10.
   - Off-line walls above do not reduce it (6.7 % at zero 5).
   - Mechanism: to first order the Jensen deformation at degree d is the de Bruijn–Newman flow at height γ for time γ²/(8d). An isolated pair x ± iδ collides at time δ²/2, which reproduces D_sec exactly. Real neighbours shorten the collision time to log(1 + 2Sδ²)/(4S).
   - Predicted excess (g_screen.py, from the actual zeros): 1.4, 2.8, 3.3, 5.5, 5.8, 7.3, 10.8 %. Measured: 0.8, 2.2, 3.0, 4.7, 5.5, 6.8, 9.9 % for zeros 1, 2, 3, 4, 5, 7, 10. The neighbour sum is taken in its y → 0 limit, which slightly overestimates.
3. **Adversarial tests.**
   - The Davenport–Heilbronn function passes every Jensen test we could compute: all γ_DH(n) > 0 for n ≤ 110, J^{d,0} real-rooted for d ≤ 60, J^{9,n} for n ≤ 40. Its first off-line zero acts only from d ≈ 1.93·10⁴.
   - A re-derives Pólya–Schur, B the Chasse and Li T² wedges, C absolute convergence.
   - The screening idea improves a record, not RH. It needs the real zeros below T (verified) and must survive adversarial zeros above T. Walls are tested; a gap allowed by |S(t)| ≤ 7 lowers the heuristic gain to about ×1.6.

## Survey (arXiv API and abstract pages, 2005–2026)
- **Local lemma = zero-freeness.**
  - Scott–Sokal (J. Stat. Phys. 2005): the local lemma holds for (p_v) iff the multivariate independence polynomial has no zero in the polydisc |w_v| ≤ p_v. Equivalently, Shearer's condition holds: Z_{G[S]}(−p) > 0 for all S.
  - Sharp for trees: (Δ−1)^{Δ−1}/Δ^Δ (Shearer 1985).
  - The dominant singularity of log Z(−λ) lies on the negative axis (coefficients alternate in sign; Pringsheim).
  - Hierarchies between the local lemma and Shearer: 2111.08837. Cluster-expansion lemma: 1711.06797. Multi-body Dobrushin–Sokal polydiscs: 2508.12078.
- **Algorithms.**
  - Barvinok interpolation: on a zero-free neighbourhood of [0,1], log Z is fixed by O(log(n/ε)) Taylor coefficients. Patel–Regts 2017 make this polynomial-time.
  - Harvey–Srivastava–Vondrák 1608.02282: Weitz's algorithm works in every root-free polydisc.
  - Below −λ_s approximation is hard (Galanis–Goldberg–Štefankovič 2017).
- **Tree recursion and decay of correlations.**
  - Weitz 2006 (self-avoiding-walk tree, up to the uniqueness threshold λ_c(Δ)).
  - Peters–Regts 2019: zero-free neighbourhood of [0, λ_c(Δ)), by contraction of the tree map. Connective-constant version: 2604.02746 (Apr 2026).
  - Strong spatial mixing from zero-freeness: 2401.09317. Very strong spatial mixing ⇒ zero-free: 2603.17858.
- **Spectral independence.** Anari–Liu–Oveis Gharan 2020; Chen–Liu–Vigoda 2020–22: stability ⇒ spectral independence ⇒ optimal mixing.
- **Contraction proofs of real roots.**
  - Heilmann–Lieb 1972 (the vertex recursion keeps a Pick property).
  - Chudnovsky–Seymour 2007 (claw-free graphs). Bencs 1409.2527 (Christoffel–Darboux proof). 2610.41502 (real-rooted sink polynomials, claw-free).
  - Marcus–Spielman–Srivastava 2013 (interlacing families → Ramanujan graphs). This is the only Riemann-hypothesis-type statement the machinery proves (Ihara zeta), and only for some member of a family.
- **Contacts with number theory (old, none on zero-freeness).**
  - Björner 1101.5704: M(n) is the reduced Euler characteristic of the complex of squarefree integers ≤ n, a wedge of spheres with Σβ_k = 2n/π² + O(n^θ).
  - Knill 1608.06877 (a Morse picture).
  - Jensen side: GORZ 2019; Chasse 2011; O'Sullivan 2007.13582; Farmer 2008.07206 (the Hermite limit is Berry's "cosine attractor"; Jensen polynomials are not a route to RH); Holland 2608.08682; 2610.38870.
- **Gap.** No paper links the local lemma, Shearer, Barvinok or Weitz to ζ, L-functions or Jensen polynomials (0 relevant hits per pairing). We found no quantitative estimate of the first failing Jensen degree beyond Chasse's sector bound, but we have not run a full novelty check.

## Transfers and tests

### A (N6): a Weitz-type tree recursion for Jensen polynomials (`a_pascal.py`, `b_fragility.py`)
- **Recursion.**
  - J^{d,n}(X) = Σ_j C(d,j)γ(n+j)X^j satisfies J^{d+1,n} = J^{d,n} + X·J^{d,n+1}: two "children", as in Weitz's R_v = λΠ_c 1/(1+R_c).
  - For the ratios r_{d,n} = J^{d,n+1}/J^{d,n}: r_{d+1,n} = r_{d,n}(1 + X·r_{d,n+1})/(1 + X·r_{d,n}).
  - Generating function: Σ_n J^{d,n}y^n/n! = (1 + X∂_y)^d F(y). Here 1 + X∂_y is a Borcea–Brändén preserver (its symbol 1 − Xw is stable), but it carries real stability only downward from F ∈ LP.
- **One-step lemma (true).**
  - Write f ≪ g for 0 > f_1 ≥ g_1 ≥ f_2 ≥ … ≥ g_d. If J^{d,n} ≪ J^{d,n+1}, then f and X·g interlace, so J^{d+1,n} = f + Xg is real-rooted (Hermite–Kakeya–Obreschkoff).
  - Write I(d) for "J^{d,n} ≪ J^{d,n+1} for all n" and H(d) for "all J^{d,n} real-rooted". Then I(d) ⇒ H(d+1), asserted on every sample.
- **Test 1: no invariant region.**
  - I(d) does not imply I(d+1), for each d = 1…5 (perturbed real-rooted moment sequences; found in 4–997 tries).
  - Explicit: γ = (1, 1, 1, 1/100) is log-concave; J^{2,0} and J^{2,1} are real-rooted but do not interlace; J^{3,0} has complex zeros.
  - In the toy of B (zero 1, δ = 0.45), J^{248,0} and J^{248,1} are real-rooted (certified) and J^{249,0} is not.
  - Every degree is a new inequality; all of them together are the Laguerre–Pólya condition (Pólya–Schur).
- **Test 2: fragility.**
  - ε(d) = min over local extrema x* of |J(x*)|/Σ_j a_j|x*|^j: the relative coefficient slack before two zeros can collide.
  - For d ≤ 120: ε ≈ 10^{−0.60d} for Ξ (shifts 0 and 40), 10^{−0.24d} for Hermite He_d (the matching polynomial of K_d, GORZ's limit), 10^{−0.48d} for Laguerre.
- **Test 3: no exact Heilmann–Lieb-type recursion.**
  - J^{d+1,0} = (AX+B)J^{d,0} − C·J^{d−1,0} fails for Ξ, with residual 2·10⁻³ (d = 2) to 3·10⁻² (d = 30). It holds exactly for Laguerre.
  - This fits Chihara's classification of orthogonal Brenke polynomials: here Σ_d J^{d,0}t^d/d! = e^t F(Xt).
- **Verdict:** dies. The only exact recursion does not contract.

### B (N6, N3): Barvinok interpolation as Taylor-coefficient criteria; reach and screening (`c_reach.py`, `e_li_reach.py`, `g_screen.py`)
- **Correspondence.**
  - Barvinok: zero-freeness in a disc ⇔ geometric control of the Taylor coefficients of log Z. For ξ the instances are Li (log ξ(1/(1−z)), whose disc maps onto Re s > 1/2) and Jensen (coefficients at s = 1/2).
  - A zero at height γ and offset δ lies a conformal distance ≈ δ/γ² outside the disc, so a coefficient test sees it only at order ∝ γ².
- **Toy (Davenport–Heilbronn-type).**
  - Θ̃ = Θ·(1 − w/w̃)(1 − w/w̄̃)/(1 + w/γ_k²), where 8ξ(1/2+z) = Θ(z²) and w̃ = (δ + iγ_k)²: Ξ with the k-th pair moved to Re s = 1/2 ± δ.
  - It keeps the functional equation and reality on the line, has no Euler product and has one off-line quartet.
  - γ(n) = 64·4ⁿ·n!/(2n)!·∫₀^∞u^{2n}Φ(u)du, by the trapezoidal rule (positive integrand). It agrees with the Arb power series to 804 digits for n ≤ 160.
  - γ_k is taken to full precision. With a double-precision γ_1 the toy acquires spurious complex zeros from about height 47.
  - First failing degree: track the off-line pair of J^{d,0}(X/d) down from the Hurwitz limit (Newton, up to 2607 digits, bisection). Certify the real-rooted side by sign counts where affordable.

| moved zero(s) (δ = 0.45 unless noted) | Chasse D_sec | first failing d | excess | screening prediction (g_screen.py) |
|---|---|---|---|---|
| γ_1 = 14.13 | 247 | 249 (J^{248,0}: 248/248 real, certified) | +0.8 % | 250 (+1.4 %) |
| γ_1, δ = 0.35 / 0.25 | 408 / 799 | 410 (J^{409,0}: 409/409, certified) / 801 | +0.5 / +0.3 % (always +2 degrees) | excess in degrees independent of δ |
| γ_2 = 21.02 | 546 | 558 | +2.2 % | 561 (+2.8 %) |
| γ_3 = 25.01 | 772 | 795 | +3.0 % | 797 (+3.3 %) |
| γ_4 = 30.42 | 1143 | 1197 | +4.7 % | 1206 (+5.5 %) |
| γ_5 = 32.94 | 1339 | 1412 | +5.5 % | 1417 (+5.8 %) |
| γ_7 = 40.92 | 2067 | 2207 | +6.8 % | 2219 (+7.3 %) |
| γ_10 = 49.77 | 3059 | 3362 | +9.9 % | 3389 (+10.8 %) |
| γ_1, γ_2 (two pairs) | 247 | 249 | +0.8 % | — |
| γ_1, γ_2, γ_3 | 247 | 250 | +1.2 % | — |
| γ_5, γ_6 | 1339 | 1422 | +6.2 % | — |
| wall γ_5 … γ_12 | 1339 | 1429 | +6.7 % | — |

- **Shifts.** For zero 1, every shift n = 1, 2, 3, 5, 10 stays real-rooted up to d = 300 (certified): the derivatives lose the off-line pair, as Pólya–Wiman predicts. The failure sits at n = 0 only.
- **Li.** λ_n(toy) first turns negative at n = 4079, 9640, 14449 for zeros 1, 2, 3, about 16–18 times the Jensen degree. The prediction (γ²/δ)·log(main term/2) is within 5 %.
- **Screening law.**
  - **Exact form.** Σ_d J^{d,0}(X)t^d/d! = e^t·Θ(Xt) is a Brenke generating function. Cauchy's formula on |t| = d gives
    J^{d,0}(X/d) = (d!·e^d/(2π·d^d)) ∫_{−π}^{π} e^{d(e^{iθ} − 1 − iθ)} Θ(X·e^{iθ}) dθ.
    This is an average of rotated copies of Θ with a complex weight ≈ e^{−dθ²/2}.
  - At height γ a rotation by θ is an imaginary shift by ≈ γθ/2, so this is the de Bruijn–Newman flow at time γ²/(8d).
  - To first order in 1/d, J^{d,0}(X/d) = Θ(X) − (1/2d)X²Θ''(X) + ….
  - With X = z², z = ix, near height γ this is ∂_t Ξ = −∂_x²Ξ (the de Bruijn–Newman direction) for the time t = γ²/(8d), plus a harmless dilation.
  - An isolated pair x ± iδ becomes real at t = δ²/2, i.e. at d = γ²/(4δ²) = D_sec: this is Chasse's bound to leading order.
  - Real neighbours add y' = −1/y − 2yS, so the pair becomes real at t_c = log(1 + 2Sδ²)/(4S) and d* ≈ D_sec·2Sδ²/log(1 + 2Sδ²).
  - With S computed from the actual zeros, this matches the table (last column). It also explains why the excess in degrees does not depend on δ for small δ: d* − D_sec ≈ γ²S/4.
  - Walls above do not cancel it: wall members at the same offset exert no vertical force on one another, and their mirror images still attract.
- **The idea (N6).**
  - The record's sector lemma uses only "all zeros in the sector". A screened version would also use "the zeros below T are real and known".
  - At T = 3.000175·10¹², with ρ = log(T/2π)/2π = 4.28 and the mean GUE neighbour sum S = ρ²·I_GUE (I_GUE = ∫(1 − (sin πu/πu)²)u⁻²du = 2π²/3 = 6.5797, numerically), the law gives a reach of 5.9× the record if only the zeros below T screen (S = 60, 2Sδ² = 17), and 9.5× with both sides (S = 120): about 0.9–1.5·10²⁶ instead of 1.6·10²⁵.
  - What a proof needs:
    1. A rigorous form of the Jensen ↔ de Bruijn–Newman correspondence beyond first order, or a direct Obreschkoff-type estimate with real neighbours.
    2. A worst case over zeros above T (count fixed by N(t) up to the S(t) bound; offsets ≤ 3/8). Only walls were tested; clusters and large gaps just above T were not. A gap of about 14 mean spacings (|S(t)| ≈ 7 at this height) still leaves a heuristic factor of about 1.6, because the zeros beyond the gap keep pulling (S ≈ 4ρ²/14 ≈ 5).
    3. All shifts n. In the toys the shifts are easier (Pólya–Wiman), but an infinite wall above T need not lose its off-line zeros.
  - The heat-flow reading of the Jensen deformation is known in kind (Pólya; Farmer 2008.07206, Berry's "cosine attractor"). The quantitative screening estimate of the first failing degree, and its use for the record, are ours; novelty not yet checked.
  - Suggested first step for the N6 owners:
    1. Prove the screened inequality for one off-line pair among known real zeros (the toy setting), from the exact rotation-average form.
    2. Then add adversarial zeros above T, using the barrier arguments of the N2 chains ([DBNB2](../DBNB2/REPORT.md)), which bound exactly this de Bruijn–Newman attraction.
- **N3 by contrast.** One worst-case zero above T gives a Li reach of 174 T²; 1.3·10¹³ phase-aligned ones give 92 T². Our record 105.77 T² lies between, so there is no comparable slack.
- **Verdict:** known wedges reproduced; one N6 lead (screening). It improves a record, not RH.

### C (N7, N10, N12): the local lemma over primes; Mertens as a hard-core gas at fugacity −1 (`d_chaingas.py`)
- **Local lemma over primes.**
  - 1/ζ(s) = Z_G(−p^{−s}) for the edgeless graph on primes. In Cartier–Foata terms ζ is the heap generating function, and Λ(n)/log n are the pyramid weights.
  - Shearer's polydisc is absolute convergence, σ > 1 (explore1 D).
  - Truncation n ≤ x makes A_p = {p | n} and A_q exclusive when pq > x, a threshold graph. Its Shearer threshold is σ_s(x) = 1.052, 1.060, 1.055, 1.049, 1.044, 1.039 for x = 10², …, 10⁷, tending to 1.
  - **Sanity test (function field).** Over F_q(t), 1/Z(u) = Π_P(1 − u^{deg P}) reaches only |u| < 1/q, while RH, a theorem there, lives on |u| = q^{−1/2}.
- **Divisor-chain gas.**
  - H_n is the incomparability graph of the squarefree integers in [2, n] under divisibility; its independent sets are the divisor chains.
  - Z_{H_n}(λ) = Σ_{m≤n} μ²(m)F_{ω(m)}(λ) (Fubini polynomials), and Philip Hall gives Z_{H_n}(−1) = M(n). Checked for n = 10², …, 10⁷: 1, 2, −23, −48, 212, 1037.
  - Shearer's point λ_s satisfies λ_s·V = 1.026 → 1.000, with V ≈ 0.61n vertices. So fugacity −1 lies 0.61n times beyond it, where approximation is hard for general graphs.
  - Morse matchings and homotopy types (as for hard squares at −1) bound |Z(−1)| only by the total Betti number 2n/π² (Björner).
  - Z_{H_n} has a complex pair near −1 from n = 10⁵ on (−1.131 ± 0.052i at 10⁷). This is explore1-C's Lee–Yang-edge reformulation again.
  - Other fugacities λ give, heuristically, D_λ(s) ≈ 1/(1 − λ(ζ(s) − 1)). Its singularities are a-points of ζ, which exist in σ > 1 for every a ≠ 0 (Bohr).
- **Verdict:** dies at once.

## Adversarial tests
1. **Davenport–Heilbronn (`f_dh.py`).**
   - Setup: ξ_DH(s) = (5/π)^{(s+1)/2}Γ((s+1)/2)f(s), f = ((1−iκ)/2)L(s,χ) + ((1+iκ)/2)L(s,χ̄), χ mod 5.
   - Results: odd Taylor coefficients vanish to 10⁻¹²⁰³; all γ_DH(n) > 0 for n ≤ 110; J^{d,0}_DH is real-rooted for d ≤ 60 and J^{9,n}_DH for n ≤ 40 (certified).
   - A, B and the screening law use only zero locations and Taylor coefficients. So they hold for DH as well, below the reach of its first off-line zero 0.808517 + 85.699348i: D_sec = 19,290, plus screening.
   - Hence the screening idea can only ever raise a finite-range record; it cannot separate ζ from DH. C needs an Euler product, which DH lacks.
2. **Known wedges.**
   - A only re-derives Pólya–Schur. B reproduces Chasse's wedge and Li's T² range, and shows Chasse is not sharp above the lowest zero. C reproduces σ > 1 (|u| < 1/q).
   - Holland's n³log²(n+2) ≥ K·d⁵ and GORZ's large-n regime lie in the shift direction, where the toys show no failures at all.

## Ranking (RH relevance × chance × gain)
| rank | transfer | node | relevance | chance | gain | result |
|---|---|---|---|---|---|---|
| 1 | B: screened sector lemma | N6 | medium (a record, not RH) | low–medium: heuristic matches the toys; the proof must handle adversarial zeros above T and all shifts | heuristic ×1.6 on the N6 record (worst case allowed by S(t); ×6 for a typical zero), if provable | lead |
| 2 | A: Weitz/Pascal recursion with a contraction region | N6 | high | 0: no invariant region; slack 10^{−0.6d}; no exact recursion | RH | dies |
| 3 | B′: Li reach | N3 | medium | 0: the record is at its reach | — | known |
| 4 | C: local lemma over primes; chain gas at −1 | N7, N10, N12 | high | 0 | — | σ > 1 only; |M(n)| ≤ 2n/π² topologically |

## Suggested NOGO line
- explore15 (Orca, exploration: hard-core lattice gas, Lovász local lemma/Shearer, Barvinok interpolation, Weitz recursion → N6, N3, N7, N10, N12) IDEA (N6, record only): screened sector lemma — Chasse's bound is sharp only for the lowest zero; real neighbours add the de Bruijn–Newman factor 2Sδ²/log(1+2Sδ²) (toys: +9.9 % at zero 10; heuristic ×1.6 at T = 3·10¹² in the worst case allowed by S(t)). Dead: the Weitz-type Jensen recursion J^{d+1,n} = J^{d,n} + X·J^{d,n+1} has no invariant region (slack 10^{−0.6d}); the local lemma over primes stops at σ = 1 (|u| < 1/q over F_q(t)); M(n) is the divisor-chain gas at fugacity −1, 0.61n beyond Shearer, total Betti number 2n/π². [ExploreLatticeGas](Riemann/Deep/ExploreLatticeGas/REPORT.md)

## Files
- `a_pascal.py` (198 s): Pascal recursion, one-step lemma, counterexamples to inheritance.
- `b_fragility.py` (27 s): fragility ε(d), Favard test.
- `c_reach.py`: toys, first failing Jensen degree, shifts, multi-pair and wall configurations. Run as `c_reach.py K [nocheck] [shift] [cfg=k:δ+…]`. Runs: K = 1100 (330 s), 600 shift (261 s), 1700 configurations (≈ 90 s), 3830 for zeros 7 and 10 (1414 s). All < 100 MB.
- `d_chaingas.py` (197 s, < 300 MB): Shearer thresholds over primes; divisor-chain gas.
- `e_li_reach.py` (1 s): Li reach of the toys and worst cases above T.
- `f_dh.py` (21 s): Davenport–Heilbronn coefficients and Jensen tests.
- `g_screen.py` (seconds): screening prediction from the actual zeros; extrapolation to T.
- `outputs.txt`: raw outputs.
- Commands: `uv run --no-project --with python-flint python3 -u <script>` (add `--with mpmath` for a_pascal, e_li_reach, g_screen; `--with numpy --with mpmath` for d_chaingas).

RH relevance: N6 (lead: screened sector lemma for the record), N3; N7, N10, N12 dead ends.
