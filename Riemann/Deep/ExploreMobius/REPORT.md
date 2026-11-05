NO IDEA: Möbius-disjointness bounds for mean-zero observables also hold for g ≡ 1, so they cannot bound M(x); and every family that passes through the constant returns exactly the zero-free exponent fed into its constant part: for the digit twists e(α s₂(n)) the two-constants inequality on the arc is an identity for the frequency-0 part (1+z)^λ (the sup form even loses log₂(1+sin πδ)), for rotations the major arcs need the half-plane ((2+Θ)/3 ≥ Θ), and the size of M(X) sits in correlations at shifts h ≳ X/16 (measured to 2^22), out of reach of the 2026 two-point Chowla theorem, while even uniform square-root Chowla would give only Θ = 3/4.

# Exploration explore6: Möbius disjointness (Sarnak), Chowla and Möbius pseudorandomness → RH tree (2026-11-05)

Field: Sarnak's conjecture, Chowla, Gowers uniformity of μ, digit and automatic sequences, quantitative rates. Brief: .scratch/orca/explore6-brief.md. Notation: M(x) = Σ_{n≤x} μ(n); Θ = exponent of M; B_2 = (777 − √1081)/852 = 0.8733819666 (the tree's Lean-checked half-plane for ζ and every Dirichlet L); c(α) = −log₂|cos πα|. Derivations: [derivations.md](derivations.md).

## Survey (one page)
- **Sarnak's conjecture.** Proved for rotations (Davenport 1937, uniform (log x)^{−A}), nilsequences (Green–Tao 0807.1736, (log N)^{−A}), horocycle flows (Bourgain–Sarnak–Ziegler 2013: Kátai-type criterion, no rate), automatic sequences (Müllner 1602.03042, qualitative), AC0 and the whole Walsh system (Green 2012; Bourgain 1109.2784), skew products (Liu–Sarnak, Wang; Lau–Ma 2602.03399), systems with singular spectrum (el Abdalaoui 2006.07646). Survey: Ferenczi–Kułaga-Przymus–Lemańczyk 1710.04039. el Abdalaoui 2607.29275: the Möbius flow has entropy (6/π²) log 3. Tao: log-Sarnak ⟺ log-Chowla. Algom–Wang 2406.06956: the decay in log-Sarnak can be arbitrarily slow.
- **Chowla.** Matomäki–Radziwiłł–Tao 2015 (averaged), Tao 2016 (log, k = 2), Helfgott–Radziwiłł (1/√log log), Pilatte (JAMS 2026, (log)^{−c}). **OpenAI, "Ordinary two-point correlations of multiplicative functions" (24 Sep 2026, github openai/math, 87 pp.): the ordinary two-point Chowla conjecture, with a (log X)^{−c} saving, c absolute, for each fixed pair of affine forms.** The method uses high traces of a divisor graph. Follow-ups: Wang–Yi 2610.16394, Lü–Wang–Wu 2610.20242 (Σ Λ(n)μ(n+2) = o(x), the parity-problem model case), Huo–Lian–Ouyang 2610.12869, Wang–Yang 2610.34912. Chinis 2105.14653: Siegel zeros imply Chowla and Sarnak along a subsequence.
- **Short intervals, Gowers norms.**
  - Matomäki–Teräväinen 1911.09076: all intervals, θ > 0.55. Matomäki–Shao–Tao–Teräväinen 2204.03754: nilsequences in all intervals of length X^{5/8+ε}, (log)^{−A}. MRSTT 2411.05770: almost all intervals, X^{1/3+ε}. MRTTZ 2007.15644: averaged version.
  - The brief's 2610.09567 is Matomäki–Pandey–Pliego–Teräväinen–Wang. Λ minus its Cramér model is Gowers-uniform in all [X, X+X^{3/5+ε}], and μ cancels in all [X, X+X^{19/35+ε}]. Method: mean values of triple products of Dirichlet polynomials with Guth–Maynard large values.
  - Pliego–Wang 2610.11487: Σ_{x<n≤x+H} μ(n)e(nα) ≪ H(log x)^{−1/3+ε} for H ≥ x^{0.58+ε}. Teräväinen 2610.03707: first quantitative short-interval Gowers bounds, H ≥ X^{5/8+ε}.
  - On [N]: Tao–Teräväinen 2107.02158 ((log log N)^{−c}); Leng 2212.09635 (U³: (log N)^{−A}, ineffective).
- **Rates.** Every result whose family contains the constant, a character or n^{it} has at most log-power (or Vinogradov–Korobov) savings. Power savings exist only for families bounded away from the major arcs:
  - badly approximable rotations, x^{4/5+ε} (Vinogradov–Vaughan); polynomial phases on minor arcs;
  - digit functions, x^{1−σ(α)} (Mauduit–Rivat 2010; Rudin–Shapiro 2015; Bourgain's large-weight Walsh functions); Zeckendorf digits (Drmota–Müllner–Spiegelhofer 2109.04068).
  - Under GRH, sup_α|Σ μ(n)e(nα)| ≪ x^{3/4+ε} (Baker–Harman 1991; refined by W. Zhang 2204.04613). The L¹ norm of Liouville exponential sums is ≫ X^{1/4} (Li–Sun–Wang 2611.03589).
- **Listings.** math.NT and math.DS for Oct–Nov 2026 were scanned for Möbius, Sarnak, Chowla, Liouville, nilsequence and digit titles. Nothing gives a power saving for a family through the constant, and nothing links disjointness to a zero-free region.

## Adversarial tests used
- **Davenport–Heilbronn (DH).** The mechanism must use something DH lacks: multiplicativity of μ, Λ ≥ 0, or a bounded Dirichlet inverse of the smooth sequence 1. The DH inverse coefficients are unbounded, because the DH function has zeros in Re s > 1 (Davenport–Heilbronn 1936; Titchmarsh §10.25).
- **Power-saving barrier.** Does the mechanism give x^{1−δ} against a family containing 1 or the characters without a zero-free strip as input?
- **g ≡ 1 test (added here).** For a uniquely ergodic system and ∫f = 0 we have Σ_{n≤N} f(T^n x) = o(N), at the rate of unique ergodicity (O(1) for the observable e(θ) of a rotation). So every Möbius-disjointness statement about mean-zero observables also holds for g ≡ 1. An argument that combines such statements with properties 1 shares with μ (multiplicativity, |g| ≤ 1, an Euler product, a Type I/II decomposition) would prove Σ_{n≤N} 1 = o(N). The same holds for g = χ and g = n^{it}. A mechanism must use μ's distance from these at power scale, and the only known source of that is a zero-free region.

## Test A — measured exponents to 2^27 = 1.34·10^8 (`t1_sieve_twists.py`, `t1_out.txt`)
Segmented sieve: M(10^8) = 1928 and L(10^8) = −3884 check; 103 s, 497 MB. m(K) = max_{2^K≤y<2^{K+1}} |Σ_{n≤y} μ(n)f(n)|.

| f | fitted exponent (K = 12..26) | m(K)/2^{K/2} at K = 20, 23, 26 | best proven rate |
|---|---|---|---|
| 1 (Mertens) | 0.508 | 0.54, 0.50, 0.56 | log-power classically; x^{B_2+ε} here |
| λ (Liouville) | 0.498 | 1.70, 1.68, 1.53 | same |
| e(nφ), golden | 0.459 | 1.56, 0.63, 0.67 | x^{4/5+ε} |
| e(n√2) | 0.498 | 1.21, 1.23, 1.74 | x^{4/5+ε} |
| e(n²√2) | 0.519 | 1.71, 1.54, 1.40 | power (minor-arc Type I/II) |
| e(⌊n√2⌋n√3), Heisenberg nilsequence | 0.604 | 1.65, 2.04, 1.48 | (log)^{−A} (Green–Tao) |
| Thue–Morse | 0.455 | 1.50, 0.55, 0.65 | power (Mauduit–Rivat; Bourgain) |
| Rudin–Shapiro | 0.505 | 0.93, 0.93, 1.05 | power for Λ (Mauduit–Rivat 2015) |
| e(0.1·s₂(n)) | 0.530 | 0.98, 1.01, 0.96 | x^{1−σ}, σ ≤ c(0.1) = 0.072 (T1) |
| Sturmian (golden) | 0.514 | 0.66, 0.40, 0.61 | power (rotation bound) |
| random ±1 (not zero entropy) | 0.434 | 0.57, 1.31, 1.02 | x^{1/2+ε} a.s. |

All exponents are 1/2 within the fit noise (±0.1 on 15 windows; the random control gives 0.43). The constant stays at the low end at every K (0.43–0.57), while the unimodular twists range over 0.55–2.1, as random signs do (0.57–3.0). This is a first sign of the frequency-0 rigidity measured in T2. (The Sturmian twist is smaller only because its amplitude is ≤ 0.62.) Every proven rate for a mean-zero twist is a power; for the constant, nothing beats the strip input. The numerics separate no zero-entropy twist from the constant.

## T1 — digit twists through the constant → Mertens (N7, N10) (`t1_analyze.py`, `t1_out.txt`, derivations D1–D2)
**Idea tested.** F(α) = Σ_{n<x} μ(n)e(α s₂(n)) is a polynomial of degree λ = log₂ x in z = e(α), with F(0) = M(x). For fixed α ≠ 0 it has a proven power saving. Interpolate to α = 0 from bounds on the arc ‖α‖ ≥ δ.
- **Sharp interpolation.** The two-constants inequality on the arc. The Green function at 1 is log((1+sin πδ)/cos πδ), and the harmonic measure is explicit (conformal map in D1). The sup form is Bernstein–Walsh.
- **Barrier identity.** For g ≡ 1, F₁ = (1+z)^λ has its only zero on the arc, so the two-constants bound is an equality: checked 1.0000 at every δ from 0.002 to 0.3. The sup form loses exactly a factor x^{log₂(1+sin πδ)}.
- **Closed loop.**
  - Write F_μ = M·((1+z)/2)^λ + R with R(1) = 0.
  - A proven bound |F_μ| ≤ x^{1−σ(α)} gives Θ_out ≤ 1 − ∫(σ − c) dω₁.
  - The frequency-0 part forces σ ≤ 1 − Θ + c wherever it dominates R, hence Θ_out ≥ Θ: the transfer returns the exponent fed in.
  - Type I/II proofs bound the m = 1 Type I term (the g ≡ 1 sum) in absolute value. So σ(α) ≤ c(α), and the bound returned is the trivial one (≥ x).

| δ | g/ln 2 | c(δ) | sup-form slack | (1+z)^λ, two-constants | μ (x = 2^27, Θ_true = 0.310), true F: two-constants Θ | sup-form Θ | random model (Θ_true 0.481), two-constants Θ |
|---|---|---|---|---|---|---|---|
| 0.002 | 0.0091 | 0.0000 | 0.0090 | 1.0000 | 0.333 | 0.503 | 0.489 |
| 0.02 | 0.0907 | 0.0028 | 0.0879 | 1.0000 | 0.494 | 0.585 | 0.562 |
| 0.05 | 0.2276 | 0.0179 | 0.2097 | 1.0000 | 0.675 | 0.722 | 0.685 |
| 0.10 | 0.4609 | 0.0724 | 0.3885 | 1.0000 | 0.920 | 0.955 | 0.898 |

- With the true, square-root-size values of F_μ the interpolation works well (0.33 at δ = 0.002). The obstacle is the input: proven bounds sit at or above x^{1−c(α)}.
- To beat B_2 the sup form needs σ(δ) > 0.172, 0.354, 0.588 at δ = 0.01, 0.05, 0.10. The method cap is c(δ) = 0.0007, 0.018, 0.072.
- **No hidden rigidity in the digit classes.** Exact identity: c_k(λ) = Σ μ(m) over odd m in [2^{λ−1}, 2^λ) with s₂(m) = k, because μ(2m) = −μ(m) and s₂(2m) = s₂(m). With the matching random model the residuals have rms 0.64–1.06 for λ = 15..27. A naive model suggested 0.4–0.7, which was this identity, not rigidity.
- **DH.** The interpolation is pure function theory. The input uses Vaughan's identity for the bounded inverse of 1, which DH lacks. **Power saving.** The family contains 1 at α = 0, and the method stops at the frequency-0 component, which carries M itself.
- **Verdict: dead.** By-product: square-root cancellation on ‖α‖ ≥ δ for every δ > 0 would imply RH, but that hypothesis is stronger than RH and open even under GRH.

## T2 — two-point Chowla (now a theorem) → Mertens (N10, N7) (`t2_chowla_fft.py`, `t2b_scaling.py`, `t2_out.txt`, `t2b_out.txt`, D3)
**Idea tested.** OpenAI's 2026 theorem gives a log-power saving for each fixed shift. Can correlations bound M? Identities: M(X)² = C(0) + 2Σ_{h<X} C_X(h), and ΣM(x)² = D(0) + 2ΣD(h). Here Q(H) = 2Σ_{h≤H}C/C(0), R(H) = (D(0) + 2Σ_{h≤H}D)/D(0), and V(H) = window variance/(6H/π²).

| H = X/2^j | j = 12 | 10 | 8 | 6 | 4 | 3 | 2 | 1 |
|---|---|---|---|---|---|---|---|---|
| V(H), μ, X = 2^18 | 0.96 | 0.91 | 1.04 | 1.04 | 0.85 | 0.44 | 0.14 | 0.05 |
| V(H), μ, X = 2^20 | 0.95 | 0.98 | 1.06 | 1.46 | 0.82 | 0.41 | 0.08 | 0.14 |
| V(H), μ, X = 2^22 | 0.97 | 1.00 | 0.95 | 0.97 | 0.94 | 0.51 | 0.24 | 0.11 |
| V(H), random signs, 32 seeds: mean (sd) | 1.01 (0.03) | 1.01 (0.08) | 1.02 (0.15) | 1.05 (0.32) | 1.13 (0.68) | 1.12 (0.83) | 1.35 (1.60) | 1.55 (2.02) |
| Q(H), μ, X = 2^22 | −0.02 | −0.04 | 0.02 | 0.02 | −0.49 | −1.36 | −1.24 | −0.70 |
| Q(H), random: mean (sd) | 0.00 (0.03) | 0.00 (0.05) | −0.01 (0.11) | 0.02 (0.23) | 0.03 (0.57) | 0.22 (0.81) | 0.32 (1.28) | 0.27 (1.83) |
| R(H), μ, X = 2^22 | 1.00 | 0.97 | 0.97 | 0.76 | 0.37 | −0.40 | −0.10 | 0.17 |
| R(H), random: mean (sd) | 1.01 (0.03) | 1.00 (0.06) | 1.00 (0.13) | 1.03 (0.28) | 0.99 (0.57) | 1.09 (0.81) | 1.24 (1.17) | 1.25 (1.58) |

- **Locally random, globally rigid.** For H ≤ X/64, μ matches random signs, as the Good–Churchhouse prediction 6H/π² says. Below 0.6 it drops only at H ≥ X/8, at all three X. So the location scales with X, not with a fixed H. At H = X − 1, Q = −0.980 and R = 0.040: ΣM² falls to 4% of its diagonal, through shifts h > X/64.
- **Mean square.** The log-mean of M(x)²/x is 0.0294, 0.0301, 0.0295 at X = 2^18, 2^20, 2^22. The sum Σ_ρ 1/|ρζ'(ρ)|² over 300 zero pairs (mpmath) is 0.0287, the prediction under RH and simple zeros. Random signs give 0.56 (sd 0.35), against 6/π² = 0.608. So M's mean square is about 20 times below random, a property of the zeros, not of local statistics. One sequence at large H is noisy (control sd 0.8–2.0). The evidence is the agreement across X and with the zero sum.
- **Triangle route.** √(C(0) + 2Σ|C(h)|) = X^{0.7355} for μ and for random signs alike, while |M(X)| = X^{0.356}. Uniform square-root Chowla over all h ≤ X gives only Θ ≤ 3/4. Uniform |C_X(h)| ≪ X^{1−η} gives Θ ≤ 1 − η/2, so beating B_2 needs η > 0.2532 for every h up to X. That would imply a strip, and no known argument derives it even from GRH.
- **DH.** The identities hold for every sequence; the Chowla inputs use multiplicativity, so they pass. **Power saving.** Fixed-shift results stop at log savings, and M needs signed cancellation over shifts h ≍ X.
- **Verdict: dead (range barrier).** The OpenAI theorem changes nothing for M.

## T3 — rotations and the Sarnak toolbox at the N7 frontier (N7) (`t3_davenport.py`, `t3b_universal.py`, `t3_out.txt`, `t3b_out.txt`, D4)
- **Rotations, a family through the constant.**

| x | 2^12 | 2^14 | 2^16 | 2^18 | 2^20 | 2^22 |
|---|---|---|---|---|---|---|
| sup_α\|S_x(α)\|/√x | 3.13 | 3.06 | 2.74 | 2.85 | 2.80 | 2.96 |
| sup/√(x log x) | 1.09 | 0.98 | 0.82 | 0.81 | 0.75 | 0.76 |
| \|M(x)\|/√x | 0.30 | 0.25 | 0.05 | 0.05 | 0.25 | 0.11 |
| argmax, nearest a/q (q ≤ 64), offset·x | 11/30, 0.1 | 1/30, 0.1 | 1/11, 28 | 17/43, 29 | 1/29, 213 | 8/35, 131 |

  - From x = 2^16 on, the maximum is not at a rational with small denominator (28–213/x away from the nearest a/q, q ≤ 64), as for random polynomials.
  - Proven: classically (log x)^{−A}. With Re s > Θ zero-free for all Dirichlet L, the textbook major/minor split gives x^{(2+Θ)/3+ε}: x^{0.958} at Θ = B_2 (D4).
  - The family contains α = 0, so its exponent is at least Θ. It can be fed by N7 but cannot feed it.
- **Universality cap.** Kátai–BSZ, Daboussi and Montgomery–Vaughan use only "multiplicative, |g| ≤ 1". Take g(p) = e(−pφ) and g(p^k) = 0. Then |Σ_{n≤x} g(n)e(nφ)|/π(x) = 1.0048, 0.9927, 1.0010, 1.0013 at x = 10^4 … 10^7, while μ gives 0.15–0.93√x. Such criteria therefore cap at x/log x even for the golden rotation.
- **At BeyondB2's levers.** These need a power saving Z^{−λ}: λ ≥ 0.000382 for the direct lever, σ ≥ 0.0035 for the dual lever, to reach 0.873. The Sarnak toolbox offers three tools, and none gives a power:
  - Kátai–BSZ: rate at most 1/log, capped as above.
  - Type I/II bounds for μ against a sum of characters Σ_u b_u ψ_u: the bilinear form factors over u and reduces to the family large sieve, the route BeyondB2 measured to lose Z^{0.0333} at the contact (its LS1).
  - OpenAI's divisor-graph trace method: it proves log-power savings only, and I found no step in it that would give a power of Z.
- **DH.** Davenport's major arcs need zero-free regions for L(s, χ), which DH lacks; the minor arcs use Vaughan's identity for μ. **Power saving.** It stops at the major arcs: Σ μχ is the strip itself.
- **Verdict: dead.** Side remark (routine, off-tree, novelty not checked): sup_α|Σ_{n≤x} μ(n)e(nα)| ≪ x^{(2+B_2)/3+ε} unconditionally in this tree's setting.

## Ranking (RH relevance × chance × gain)
| transfer | node | chance | gain if it worked | verdict |
|---|---|---|---|---|
| T1 digit twists → M (two-constants on the arc) | N7, N10 | none: the loop is closed; Type I/II savings ≤ c(α) | a strip below B_2 needs 0.127 of extra saving in harmonic mean | dead (barrier identity) |
| T2 two-point Chowla → M | N10, N7 | none with fixed-shift theorems; needs power-saving Chowla uniform up to h ≍ X (open even under GRH) | uniform η > 0.2532 beats B_2; square-root gives 3/4 | dead (range barrier) |
| T3 rotations, Kátai–BSZ, Type II at BeyondB2's levers | N7 | none: the universality cap is x/log x; Type II is the large sieve | BeyondB2 needs Z^{−λ}, λ ≥ 0.000382 | dead; reverse corollary (2+B_2)/3 only |

**New relative to R32, R58 and R317** (which said "no power saving", "circular" and "controls no prefix sum"):
- the exact two-constants loop and the Bernstein–Walsh slack log₂(1+sin πδ);
- the 2026 two-point Chowla theorem and its range barrier (measured location h ≳ X/16, triangle cap 3/4);
- the g ≡ 1 test;
- the Sarnak toolbox mapped onto BeyondB2's levers.

## Suggested NOGO lines
- explore6 (Orca: Möbius disjointness / Sarnak / Chowla) NO IDEA. Mean-zero disjointness bounds also hold for g ≡ 1. Families through the constant return the exponent fed in: for digit twists the two-constants inequality on the arc is an identity at (1+z)^λ, with sup-form slack log₂(1+sin πδ); for rotations (2+Θ)/3 ≥ Θ. M(X)'s size sits at shifts h ≳ X/16 (measured to 2^22), beyond the 2026 two-point Chowla theorem, and uniform square-root Chowla would give only 3/4. [ExploreMobius](Riemann/Deep/ExploreMobius/REPORT.md)
- Possible new FILTER line next to the DH filter: inputs that g ≡ 1 (or χ, n^{it}) also satisfies cannot bound M(x).

## Files
`Riemann/Deep/ExploreMobius/`:
- `t1_sieve_twists.py` → `t1_sieve_out.json`; `t1_analyze.py` → `t1_out.txt`
- `t2_chowla_fft.py` → `t2_out.txt`; `t2b_scaling.py` → `t2b_out.txt`
- `t3_davenport.py` → `t3_out.txt`; `t3b_universal.py` → `t3b_out.txt`
- `derivations.md` (D1–D4)

Each script ran alone, single-threaded, under `perl -e 'alarm 1800'` with `uv run --no-project --with numpy [--with mpmath]`. Longest run 103 s; peak 692 MB (`t2_chowla_fft.py`). Scratch: `.scratch/explore6/STATE.md`; downloads deleted.

RH relevance: N7 (zero-free strips: no route from disjointness; the strip only feeds the field), N10 (Mertens: the size of M lives at frequency 0 and at shifts h ≍ X); N5 not touched.
