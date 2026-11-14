BARRIER: geometry of numbers sees a rank-2 lattice zeta only below its smoothing height T_s = 2π/λ₁² = 2πy. There RH holds for every lattice, Euler product or not. Off-line zeros of non-Euler lattices are born just above the multiples kT_s, at the Airy extrema of the k-th dual frequency. Every Euler-product lattice has T_s ≤ π√163 ≈ 40.1 (Heegner–Stark), so smoothing and transference arguments say nothing about ζ above height 40 (node N12).

## Transfer
- Tool (b), the smoothing parameter: the Gaussian tail of the dual lattice, as in Banaszczyk's mass lemma and Micciancio–Regev.
- For the unimodular lattice L_z (z = x + iy reduced, λ₁ = y^{−1/2}), Chowla–Selberg is Poisson summation along the shortest vector:
  E*(z,s) = a_y(s) + 4√y Σ_{k≥1} k^{s−1/2} σ_{1−2s}(k) K_{s−1/2}(2πky) cos 2πkx, with constant term a_y(s) = ζ*(2s)y^s + ζ*(2s−1)y^{1−s}.
  The k-th Bessel term is the k-th dual frequency, i.e. the smoothing error of the line sums. At height t it is negligible until t = 2πky, the complex-frequency smoothing parameter.
- Node: N12, which lattice zetas satisfy RH. Existing lines: R046 (along deformations, zeros leave the line only by collision), R073 (Epstein never easier), explore21 (Davenport–Heilbronn zeros for D = −15). None says where the off-line zeros start.

## Statements
1. **Lemma K** (rigorous; contour shift v → v + iα with sin α = t/u in K_ν(u) = ½∫e^{−u cosh v + νv}dv). For 0 ≤ t < u and |μ| ≤ ½:
   e^{πt/2}|K_{μ+it}(u)| ≤ e^{t·arccos(t/u)} K_μ(√(u²−t²)) ≤ √(π/2w)·e^{−u·g(t/u)}, where w = √(u²−t²) and g(a) = √(1−a²) − a·arccos a > 0.
2. **Proposition** (proof sketch in DERIVATION.md §4; routine, possibly folklore). For every δ > 0 there is y₀(δ) such that for y ≥ y₀ every non-real zero of E(z,s) with |Im s| ≤ (1−δ)·2πy lies on Re s = ½ and is simple. The proof combines three things:
   - Lemma K;
   - Rouché against a_y on the boxes [0,1] × [t_j, t_{j+1}], where t_j are the extrema of cos θ_a;
   - RH for a_y when y ≥ 1 (Hejhal 1990, Lagarias–Suzuki 2006, Ki), and the sign changes of F = E*/(2√y|ζ*(1+2it)|) = cos θ_a + β.
3. **Birth law** (numerical, every case verified). The first off-line zero ρ₁ has Im ρ₁ ≈ kT_s + |a′_j|·(kT_s/2)^{1/3}, where a′_j are the zeros of Ai′ (|a′_j| = 1.019, 3.248, 4.820): the extrema of the Airy turning point of K_{it}(2πky). Later births cluster in the same windows. In 20 of 24 verified births, Im ρ lies in [k − 0.06, k + 0.30]·T_s; that window covers 36 % of the line, so the binomial tail is 3·10⁻⁶. In every case Re ρ₁ − ½ ≤ 0.44 for y ≥ 1.9; y·(Re ρ₁ − ½) ranges over 0.21–3.19, median ≈ 0.9.

## Tests
- **Smallest cases, validation** (t1_validate.py). Chowla–Selberg agrees with the ζ·L products to 10⁻¹⁵ for D = −20 and −15. It reproduces Potter–Titchmarsh's zero of x² + 5y², 0.93296970 + 15.66824953i (1.115·T_s), and explore21's zero of x² + xy + 4y², 0.80001099 + 12.03859863i (0.989·T_s).
- **Exact counting** (count.py). Argument principle for s(s−1)E*, against sign changes of F. Integers to 10⁻⁴ at every checkpoint. Pitfall: ζ·L(χ_D) has nearly coincident zeros; missed close pairs look like off-line zeros. Fixed by local refinement and by Newton localisation.
- **Scale** (t4_onset.py, analysis.py, logs/analysis.log): 45 non-Euler lattices, first off-line zero located in 44.

  | family | y range | Im ρ₁ / T_s |
  |---|---|---|
  | 18 principal forms of h = 2 (incl. D = −15, −20) | 1.94–10.33 | 0.982–1.166; D = −403: 2.059 (k = 2) |
  | D = −20 non-principal form 2x² + 2xy + 3y² | 1.118 | 1.086 (zero at 2.079 + 7.629i, Re > 1) |
  | rectangular, z = iy | 2.5–40 | 0.954–1.119; y = 15, 30, 40: 2.229, 2.025, 2.078 (k = 2) |
  | rhombic, z = ½ + iy | 2.5–30 | 0.9999–1.322; y = 40: none up to 2.58·T_s |

  - Airy fit at y ≥ 15, with ξ = (Im ρ₁ − kT_s)(2/kT_s)^{1/3}:
    - ξ = 1.222, 1.016, 0.968, 0.943, 1.119, 0.877, 0.835, against |a′₁| = 1.019;
    - ξ = 3.109 against |a′₂| = 3.248;
    - ξ = 4.734 against |a′₃| = 4.820.
  - At small y, T_s < 40 and ξ spreads over −0.36 to 1.56: the hump is far above threshold, so births come early.
- **All births to 5·T_s** (t6_dist.py, then t7_verify.py: exact count per interval, fine on-line scan, Newton; all located):
  - rect y = 10: Im/T_s = 1.029, 1.066 | 1.971, 2.069, 2.120, 2.290 | 2.616 | 3.065 | 4.106, 4.122, 4.299. Nothing in (1.07, 1.97), (2.62, 3.06) or (3.07, 4.10).
  - D = −427: Im/T_s = 1.034, 1.305, 1.427 | 2.079, 2.177 | 3.048, 3.086, 3.102, 3.121 | 3.804, 3.947 | 4.040, 4.284.
- **Comparison with the node line.** TREE N12 has no Epstein onset. R046's "leave the line only by collision" now has a location: the collisions sit at kT_s + |a′_j|(kT_s/2)^{1/3}.

## Adversarial
- **ζ against Davenport–Heilbronn.** The h = 1 lattices D = −163, −67, −43 are Euler products ζ·L(χ_D) on the same family x = −½ as the h = 2 forms. |β| reaches 2.7, 3.5 and 3.2 above T_s, yet there are no off-line zeros up to 4, 3.39 and 3.49·T_s. The one apparent birth (D = −163 at 3.46·T_s) is 3 zeros, 3 of them on the line.
  - So the geometric criterion (Lemma K, Proposition, |β|) is identical for Euler and Davenport–Heilbronn lattices. It does not distinguish ζ from Davenport–Heilbronn above T_s.
  - The birth law separates them only through whether births exist, which is the arithmetic itself.
  - Off-line zeros also occur where |β| on the line stays below 1 (D = −15: 0.83; rect y = 2.5: 0.73), so the on-line smoothing test is not even necessary.
- **What ζ gains.** By Heegner–Stark, h = 1 forces y ≤ √163/2, so T_s ≤ π√163 = 40.11. The Proposition would give RH for ζ·L(χ_D) only below height ≈ 40 (Gram 1903 level). GUE and Poisson models are not relevant here: the statement is about where the off-line zeros are born, not about spacing statistics.

## Novelty (no reviewer verdict yet)
- arXiv API searches: Epstein + critical line / off-critical / zeros, Chowla–Selberg + zeros, au:Ki_Haseo, Eisenstein constant term.
- TeX read:
  - Bétermin–Šamaj–Travěnec 2110.09368: rectangular lattices, Δ ≤ 1, ρ_y ≤ 21. Their edge points 6a, 8a, 8b, 9a, 14a sit at 0.94–1.01·T_s, but no 2π/Δ law is stated.
  - McPhedran 1601.01724: trajectory of the Potter–Titchmarsh zero.
  - Lagarias–Suzuki math/0412039: a_y only; its summary of Ki (PLMS 90, 2005): truncations have all but finitely many zeros on the line.
  - Bogomolny–Leboeuf chao-dyn/9409004: statistics only.
- Abstracts: Lee 1204.6297, Gonek–Lee, Lamzouri 1907.06387, Rezvyakova 2411.18492 (counts as T → ∞, no onset).
- arXiv Oct–Nov 2026 listings (math.NT/CA/CV/SP/MP) and Zenodo "Epstein zeta" (81 hits, incl. 24032536, 24048053): no overlap.
- Not read (paywalled): Stark 1967 (real zeros), Hejhal 1986 ICM, Arenstorf–Brewer 1993, Ki 2005 in full.
- Verdict: the Proposition is probably known in spirit (Hejhal, Ki). The onset at T_s, the quantisation at kT_s and the Airy-extremum positions were not found.

## Other bridges (measured, rejected)
| bridge | what the tool computes | why no gain |
|---|---|---|
| N5: d_N as a closest-vector problem (LLL/BKZ/Babai) | integer-coefficient approximants | d_N² is a real least-squares problem; QR/Cholesky solves it exactly; integrality only worsens it; certification needs interval Gram entries, not reduction |
| N1: Gaussian test functions, Gaussian mass lemma | 2^{−n} tail bounds in dimension n | Weil windows are 1-D (Chekmachev's periodised lattice minorant is already the 1-D Poisson step); a margin needs the eigenvalue of a continuous form |
| N8/N13: Gaussian heuristic, LWE torus Fourier analysis | lattice-point and successive-minima statistics | no map to zeros of ζ; LWE smoothing is tool (b) again |

## Files
- Code: epstein_cs.py (Chowla–Selberg, arb, radius-checked), count.py (exact N(T)), t1_validate.py, t2_sweep.py (batches), t4_onset.py, t6_dist.py, t7_verify.py, analysis.py, run_all.sh, run2.sh.
- Logs: logs/*.log, logs/*.json.
- Derivation: DERIVATION.md. Peak memory ≤ 18 MB in every run.

## Suggested NOGO line
- explore48 (Orca, exploration: lattice cryptography and algorithmic geometry of numbers — Banaszczyk mass/transference, smoothing parameter, LLL/BKZ, LWE → N12, N5, N1) BARRIER: Chowla–Selberg = Poisson along the shortest vector, and Lemma K bounds the dual-frequency tails, so every rank-2 lattice zeta satisfies RH below T_s = 2π/λ₁² (Euler or not). Off-line zeros are born just above kT_s at the Airy extrema kT_s + |a′_j|(kT_s/2)^{1/3}: 44 lattices, 24 verified births. h = 1 lattices pass the same test with no births and have T_s ≤ π√163 ≈ 40, so geometry says nothing about ζ above 40. CVP/LLL is useless for d_N (a real least-squares problem). [ExploreLattice](Riemann/Deep/ExploreLattice/REPORT.md)

RH relevance: N12 (barrier: lattice-geometric control stops at the smoothing height; the birth law locates where Davenport–Heilbronn-type off-line zeros of lattice zetas come from).
