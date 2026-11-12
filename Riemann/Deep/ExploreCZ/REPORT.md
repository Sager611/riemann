NO IDEA: real-variable harmonic analysis (Calderón–Zygmund, maximal, Carleson embedding, sparse, stopping time) sees a Dirichlet polynomial only through positivity and local masses of |D|^{2k}; that class is exactly the mean-value (operator-norm + power) method. Explicit band-limited profiles attain its bound (Theorem P), and its zero-density output is exactly Ingham's 3/(2 − σ). That is already the record on [1/2, 7/10]; every gain beyond 7/10 uses the phase of Σ n^{iu}, whose spikes break the CZ size condition by factors up to 6·10³ at N = 1000 (node N7).

# explore46: Calderón–Zygmund, maximal functions and sparse domination → N7 (Orca exploration, 12 Nov 2026, 19:18–20:20 UTC)

## Transfer
Tool: (a) CZ decomposition of |D|² at height V² and the weak-(1,1) maximal inequality, plus (b)–(d), run on the large-values step of zero density (N7, through the zero detector).
D(t) = Σ_{N<n≤2N} a_n n^{−it}, |a_n| ≤ 1, G_k = coefficient norm of D^k (≤ N^{k+o(1)}), V = N^σ, T = N^τ, W ⊂ [0,T] 1-separated with |D| ≥ V, R = |W|.

## Statements (DERIVATION.md §2–5)
- Lemma L (localization at the Montgomery scale): every CZ stopping interval of |D|² at height V² > G has length < 2c·N·G/(V² − G), c = 3π. This is about 6π N^{2−2σ} when G = N. Proof: Montgomery–Vaughan's local MVT. The weak-(1,1) count of stopping intervals is the MVT. Classical.
- Class 𝒫 (what these tools consume): w ≥ 0, ŵ supported in [−log 2, log 2], and ∫_I w^k ≤ |I|G_k + c2^kN^kG_k for every interval I and every k ≤ K. Every |D|² is in 𝒫. CZ stopping, M f, weak (1,1), good-λ, John–Nirenberg, Carleson embedding of Paley–Wiener measures (the large sieve with measures) and sparse bounds of averaging operators all use only these properties.
- Theorem P: take Fejér F, ρ_k = 1.13∫F^k and B = min_k (T + c2^kN^k)G_k/(ρ_kV^{2k}). Then w = V² Σ_{r≤B/2} F(t − 2rT/B) lies in 𝒫 (under the side conditions 2T/B ≥ 20 and 4ρ_kV^{2k} ≤ c2^kN^kG_k) and has R(w) ≥ ⌊B/2⌋. Conversely every w ∈ 𝒫 has R ≲ B (Akhiezer–Krein factorization plus a reproducing kernel). So the best bound any class argument can prove is the k-th power MVT bound, R_class = N^{LV_c+o(1)} with LV_c(σ,τ) = min_k [max(τ,k) + k(1−2σ)].
- Corollaries:
  - (a) No class argument proves Huxley or Guth–Maynard where they beat LV_c.
  - (b) Through the GM §13 / Ingham zero detector, the class gives A(σ) = 3/(2−σ) exactly.
  - (c) The Halász–Montgomery kernel S(u) = Σ n^{iu} is a modulated Hilbert (CZ) kernel only for |u| < 2πN, where Montgomery–Vaughan is already lossless. Beyond 2πN it is a √N-size almost-periodic kernel, so no CZ or sparse theory applies.

## Tests (scripts and outputs in this directory; all runs under memcap 290 MB, peaks 220/129/<10 MB)
**T1 kernel** (t1_kernel.py):
- On [1, 0.9·2πN], |S(u)| ≤ 12N/u (11.6 at N = 300, 11.95 at N = 1000).
- In the bands u ≈ 2πNm, m ≤ 12: max|S|/√N rises from 1.26 to 2.52, and rms|S|/√N = 0.79 (m = 1), then 0.96–1.08.
- The CZ size ratio max|S|·u/N grows from 374 (m = 1) to 5977 (m = 12) at N = 1000.

**T2 CZ on real polynomials** (t2_cz_real.py), N = 1000, families ζ-block (n/N)^{−1/2}, μ(n), random ±1, Davenport–Heilbronn-type (period 5), χ₄, and the detector Σ_{d|n,d≤31} μ(d).
- Lemma L holds everywhere. For the ζ-block, the only stopping interval for σ ≥ 0.65 is the t ≈ 0 peak, of length 410/205/102 against L* = 2715/1269/616. So the lemma is sharp up to a factor of about 6.
- At τ = 1.6, R(σ = 0.6) for points with t > 2N:

  | family | μ | random | DH | χ₄ | detector | ζ-block |
  |---|---|---|---|---|---|---|
  | R(σ = 0.6), t > 2N | 1554 | 1285 | 974 | 693 | 2161 | 148 |

  At σ = 0.7, R = 0 for every family except the detector (17) and the ζ-block's t ≈ 0 peak (17). The no-constant bounds are MVT 4044, Huxley 15912 and GM 4295.
- Real polynomials therefore sit far below every bound. The extremal configurations are adversarial, not arithmetic.

**T3 class exponents** (t3_class.py; the density reduction is copied from ExploreKakeya/t1_energy_lp.py).
- LV_c against the best of {Huxley, GM}: (σ,τ) = (0.75,1.2): 0.7 vs 0.6; (0.8,1.2): 0.6 vs 0.4; (0.75,1.5): 1.0 vs 0.9; (0.9,1.5): 0.4 vs 0.2. These listed bounds never beat LV_c for σ ≤ 0.7 or τ ≥ 2.
- Class density at σ = 0.55, 0.60, 0.65, 0.70, 0.75, 0.80: 2.0689, 2.1420, 2.2211, 2.3067, 2.4000, 2.4980. Ingham gives 2.0690, 2.1429, 2.2222, 2.3077, 2.4000, 2.5000 (the small gaps are a grid effect).

**T4 smallest cases of Theorem P** (t4_profile.py), with the minimal G_k = N^k, K = 3, σ = 0.8:
- (N,τ) = (10,3): max over 60 window lengths of mass/(|I|G_k + c2^kN^kG_k) = 0.893.
- (30,2): 0.863; (100,1.5): 0.852; (300,1.2): 0.881.
- All four are admissible.

**T5 short range** (t5_short_range.py), N = 1000, W ⊂ [0, 0.9·2πN]:
- λ_max(Gram)/N = 4.93 (random, R = 400), 1.46 (arithmetic progression) and 6.27 (dense, R = 900).
- The large sieve bound (|I| + N)/N is 6.65, and the triangle-inequality (Schur) bound is 8.97–13.75.
- So the CZ range is attained to within 6% already. Sparse domination has nothing to gain there.

## Adversarial (ζ vs Davenport–Heilbronn, GUE/Poisson)
- The class bound depends on (N, T, V, G_k) only. The DH-type, χ₄, μ and random families give the same CZ statistics within factors 1–3 (T2), and non-arithmetic profiles attain the optimum (T4). The transferred statement carries no RH information.
- In a density argument, the ζ-specific input is the detector identity (Möbius coefficients bounded, from the Euler product). DH lacks it: it has ≫ T zeros in a strip σ₁ < σ < σ₂ inside (1/2, 1).
- Zeros (GUE or Poisson) do not enter the large-values step at all.

## Novelty
- Guth's survey 2503.07410 §5 (read in TeX):
  - The operator-norm and power methods are the best known in this range. "For σ < 7/10 and T ≤ N^{6/5} it remains an open problem to improve on this one page argument."
  - The Fu–Guth–Maldague "almost counterexample" (frequencies √(n/N)) makes the orthogonality bound sharp for σ ≤ 3/4. Decoupling methods cannot distinguish it from Dirichlet polynomials.
- Theorem P is a cheap, formal version for the CZ class. It is known in spirit, routine, and not a discovery.
- arXiv API queries found no CZ-, maximal- or sparse-function large-value estimate: "Calderon-Zygmund" with Dirichlet/zeta gave 0 hits, as did "maximal function" with "Dirichlet polynomial". Also checked:
  - Montgomery–Vaughan 1974 (Hilbert inequality = the CZ part);
  - Ingham 1940;
  - ExploreKakeya (GM's chain is pinned at (7/10, 6/5) by energy-free terms, so a CZ treatment of the spikes cannot move 30/13).

## Other suggested bridges (DERIVATION.md §7)
- **N8 (Levinson bad set).** Jensen's step loses only through small values of |f|. A constant |f| = ‖f‖₂ meets every upper-tail (CZ or weak-type) bound and has Jensen gap 0. Upper-tail tools cannot help; the step needs negative-moment input.
- **N10 (ψ in short intervals).** sup over h ≥ h₀ of the error is dominated by the Hardy–Littlewood maximal function: a factor 4 in L² and no log loss. The exponent is unchanged. Routine.

## Side corollary candidate (off-tree, harvest type; DERIVATION.md §8)
- Lacey–Mousavi–Rahimi 2101.10401 prove sharp ℓ^p-improving and sparse bounds for prime averages, with t = 1 under GRH.
- GRH enters through ψ(x;q,a) = x/φ(q) + O(x^{1/2}log²q) (l.185) inside PNTlemma (l.414/452), through Ñ = N^{1/5} (l.484, l:QQ l.497), and through the absence of exceptional characters.
- With the Lean-checked B₂ half-plane (all χ), the major-arc error Q·N^{B₂−1+ε} is still a power saving for q ≤ N^{1/2000}.
- Reviewer (verdict-only): INCORRECT as a verbatim substitution, not refuted after repair. Four repairs are needed:
  - PNTlemma needs Q ≥ N|ξ − a/q|, so the arcs must narrow to |ξ − a/q| ≤ N^{−θ−δ}, with Vinogradov's bound in between;
  - the separation step at l.535 fails for q₀ = 1;
  - the trivial range ⟨f⟩⟨g⟩ < cN^{−1/4} (l.995) must shrink to cN^{−η}, η ≈ 1/2500;
  - LMR's own GRH proof needs the same repairs.
- Status: plausible, not filed. It needs a repaired proof first.

## Ranking (RH relevance × chance × gain)
| transfer | node | chance | gain if it worked | verdict |
|---|---|---|---|---|
| CZ/maximal on \|D\|² (positivity class) | N7 | none: Theorem P | at most Ingham, already the record on [1/2, 0.7] | dead |
| CZ/sparse on the HM kernel | N7 | none: CZ part already lossless (T5); spikes not CZ (T1) | none | dead |
| CZ of log\|VΨ\| (Levinson) | N8 | none: upper tails only | none (record not Levinson) | dead |
| HL maximal for ψ(x+h) − ψ(x) − h | N10 | certain, routine | constant 4, no log | routine |
| B₂ → Lacey–Mousavi–Rahimi t = 1 | off-tree (N11 harvest) | medium: needs repaired arcs and thresholds (reviewer) | GRH removed from a sharp harmonic-analysis bound | candidate, not filed |

## Suggested NOGO line
- explore46 (Orca, exploration: Calderón–Zygmund, maximal functions and sparse domination → N7, N8, N10) NO IDEA: these tools see |D|² only through positivity and local L^{2k} masses. Admissible Fejér-spike profiles attain the k-th power MVT bound (Theorem P), so the class yields at most Ingham's 3/(2−σ) (class density 2.3067 at σ = 0.7), the record already on [1/2, 7/10]. The HM kernel is CZ only for |u| < 2πN, where Montgomery–Vaughan is already lossless (λ_max/N 6.27 against the large sieve's 6.65); beyond, spikes of size 2.5√N break the size condition by factors up to 6·10³. Upper-tail tools cannot repair Levinson's Jensen step. Side: Lacey–Mousavi–Rahimi's GRH-conditional t = 1 plausibly follows from B₂, but only after repairs that their GRH proof needs too (reviewer). [ExploreCZ](Riemann/Deep/ExploreCZ/REPORT.md)

RH relevance: N7 (zero density via large values); side lines N8, N10, N11.

## Files
REPORT.md, DERIVATION.md, t1_kernel.py/t1_out.txt, t2_cz_real.py/t2_out.txt, t3_class.py/t3_out.txt, t4_profile.py/t4_out.txt, t5_short_range.py/t5_out.txt.
