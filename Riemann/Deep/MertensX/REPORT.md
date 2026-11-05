NO RECORD: Bi–Hou (HAL hal-05838582, 26 Oct 2026) already lowered exp(1.96·10¹⁹) to exp(3.485·10¹⁸), with Gallego (Zenodo 23700212, 27 Oct) reaching the same bound, by the very parameter re-optimisation proposed here (Gaussian width k 3·10⁻⁹ → 1.5·10⁻¹⁰, 2.3·10⁶ zeros). Every candidate below 3.485·10¹⁸ published by KN or Bi–Hou falls short: the best is |K2| = 0.990 at y = 2.79·10¹⁸. Gallego's best point is the same as Bi–Hou's. A new point needs a lattice search stronger than theirs, plus more than 10⁶ certified zeros with ζ′(ρ) (about 13 core-hours in their setup). That is far beyond 2 h on one core.

Node: N10 (explicit Mertens bounds). Agent: mertensx (Orca), 2026-11-05, 21:27–21:50 UTC. No Lean, no commits.

## 1. How exp(1.96·10¹⁹) arises (Kim–Nguyen, arXiv 2502.21021)
- **Criterion** (Pintz 1987, refined by Saouter–te Riele 2014): let h(y) = 2 Σ_{0<γ<T0} α e^{−kγ²} cos(γy − ψ), where α = 1/|ρζ′(ρ)| and ψ = arg ρζ′(ρ). If |h(y)| > 1 + E, then the first counterexample satisfies x < exp(y + 2√(ky)); KN quote the weaker exp(y + √y).
  - h(y) is the Gaussian average (variance 2k in log x) of M(x)/√x, written as a sum over zeros.
  - The error E is dominated by the contour above T0, which needs a bound on 1/ζ just right of the critical line (factor u·t¹¹·e^M against e^{−kT0²}).
  - Zeros off the line are handled by Pintz's Theorem A (Bi–Hou Thm 4.8 gives the full statement).
- **KN parameters**: k = 3·10⁻⁹, threshold 1 + 6·10⁻⁸. Their code sums 98,625 zeros (γ < 74,000). Bi–Hou's README notes that Saouter–te Riele's T0 is 3·10⁵. Bi–Hou's float estimate puts the zeros KN left out near 10⁻⁹. That is not certified (it needs lower bounds on |ζ′(ρ)| up to 3·10⁵), but it is far below KN's margin of 0.0075.
- **Lattice data**: an (N+1)-dimensional lattice built from the N = 140 zeros of largest weight α*, with ν = 130, ν_y = 100, ν_t = 30.
  - It was BKZ-88 reduced, then pruned enumeration with radius factor 1.19 around the target ψ + π gave 17,406 points.
  - The best point is y = 19571878850562201959.215107, with h = −1.00745 → x < exp(1.95719·10¹⁹).
  - All data is public on Zenodo 10775723: zeros, α and ψ to 300 digits, the basis, the enumeration output, and 11,754 rescored points.
- **Free parameters**: k, T0 (E must stay below the margin), the lattice (N, ν, ν_y, ν_t, radius, target), and local polishing of y. KN §4.1 suggests re-optimising these.

## 2. Why it cannot be beaten here
1. **The target is already met by others.**
   - Bi–Hou (code: github.com/bzyzz/mertens_bound) use k = 1.5·10⁻¹⁰ and T0 = 1.3·10⁶, with 2,325,572 zeros, each with a certified Arb enclosure.
   - At v = 3484971332733132396.5383341489440461 they get K2 = 1.0003052504 ± 2.3·10⁻¹¹ with E ≤ 2.4·10⁻¹¹. Hence x < exp(3 484 971 332 733 178 124).
   - Gallego reaches the same bound, x < 10^{1.51351·10¹⁸} = exp(3.485·10¹⁸), with 16,622,627 zeros up to 7,999,035.
   - A record now needs X < 3.485·10¹⁸.
2. **No published candidate below 3.485·10¹⁸ is close enough.** Our recomputation is check (4) of `check.log`, with k = 3·10⁻⁹ after local maximisation in y.
   - KN's 9 clusters below 3.485·10¹⁸ reach at most |h| = 0.9663.
   - Bi–Hou's other 4 sub-record points reach at most 0.98632 (y = 2.789·10¹⁸). At k = 1.5·10⁻¹⁰ (their 600k-zero float values) the best is 0.990143, a shortfall of 0.0099.
3. **Lowering k helps too little.**
   - Across Bi–Hou's 25 sub-record candidates, going from k = 3·10⁻⁹ to 1.5·10⁻¹⁰ changes |K2| by a median of +0.0034 and at most +0.0064 (one outlier is −0.199).
   - Heuristically (Gonek–Hejhal, Σ_{γ>T} α² ≈ 0.19/T), zeros above 8·10⁴ add noise with standard deviation ≈ 0.002, so a +0.01 jump is about a 5σ event. Even smaller k (Bi–Hou §9, T0 ≈ 10⁷) does not close a 0.01 gap.
4. **A new point is out of budget.**
   - KN used C++/NTL enumeration with BKZ-88. Bi–Hou ran GPU stream enumeration over many basis ensembles and kept 1,919 clusters with |h| ≥ 0.96; at small k, their best point below 1.957·10¹⁹ is the record itself.
   - Scoring candidates at small k needs about 6·10⁵–2·10⁶ zeros with ζ′(ρ). Bi–Hou's zero cache takes ~50 min on 16 workers.
   - Here python-flint takes 0.09 s per zero (`acb.zeta_zero`, height 7.4·10⁴) plus 0.025 s per ζ′(ρ). That is about 3 h for the first 98,625 zeros alone.

## 3. Own verification (`check.py` → `check.log`, 125 s, 322 MB, one core)
- N(74000) = 98625 (Arb `zeta_nzeros`). KN's data matches Arb at 40 indices (1–20 plus 20 spread to 98,625) to 4·10⁻⁵⁶ in γ and 4·10⁻⁵⁵ in α and ψ.
- h_StR(19571878850562201959.215107) = −1.00745157222 ± 2·10⁻¹², in Arb ball arithmetic on KN's data. This matches KN (−1.007) and Bi–Hou (−1.0074515722).
- Our polished values reproduce Bi–Hou's k = 3·10⁻⁹ column to 6 digits: 0.995703, 0.986321, 0.982871, 0.981097, 0.980674.
- Not certified by us: Bi–Hou's k = 1.5·10⁻¹⁰ values, which are taken from their CSV.

## 4. Inputs
- KN data: Zenodo 10775723 (`zeroes.txt`, `alpha.txt`, `beta.txt`: first 100,000 zeros to 300 digits; `140-119.OUT`).
- Bi–Hou: `search/results/subrecord_peaks_summary.csv` and `certification/k2_certificate.json` (github.com/bzyzz/mertens_bound, HEAD on 2026-11-05).
- Tools: python-flint 0.9.0 (FLINT/Arb) and numpy via uv. RH height used: none beyond 74,000 (Arb-certified zero indices).

## 5. Rerun (≈2.5 min including the 91 MB download)
`OMP_NUM_THREADS=1 uv run --no-project --with python-flint==0.9.0 --with numpy python3 Riemann/Deep/MertensX/check.py /tmp/mertensx-data`, then delete `/tmp/mertensx-data`.

## 6. Novelty pass (2026-11-05)
- `lst.py 'Mertens' math.NT --zenodo 'Mertens conjecture'`: 0 arXiv title hits for 2026-10/11. Its Zenodo branch listed unrelated recent records, so the Zenodo API was queried directly.
- arXiv API (2025–26: `all:"Mertens conjecture"`, `abs:Mertens AND abs:counterexample`, `ti:Mertens`, `abs:"Mertens function"`): the only counterexample bound is 2502.21021 (KN). Also 2607.07566 (M(10²⁵)), already in TREE.
- HAL hal-05838582 (Bi–Hou, 26 Oct 2026): x < exp(3.485·10¹⁸). **TREE N10 is stale here.**
- Zenodo 23700212 (Gallego, 27 Oct 2026): the same bound, independently deposited, with 833 MB of data.
- Zenodo 23651459 (Weisz, 25 Oct 2026; unrefereed, AI-produced, two verifiers): lim sup M(x)/√x > 2.06175 and lim inf < −2.06471. This beats Hurst's 1.826054 / −1.837625, which TREE N10 still lists. **Stale too.**
- Zenodo 23405928 (Gaur, 15 Oct 2026): the naive Mertens analogue fails for Q(√−3) and Q(√5). This is a side result for L-function analogues.
- HAL hal-05675375 (Yarman et al., Jun 2026) says the 1985 disproof does not replicate (their LLL search reaches only 0.53). Check (3) above and Bi–Hou's certified K2 contradict that claim. It is a candidate for REFUTATIONS.md.

RH relevance: N10 (size of M(x)/√x; the Mertens conjecture is a stronger-than-RH statement).
