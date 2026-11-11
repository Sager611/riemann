RESULT: CCM's first missing step (QW_λ has a simple, even ground state) certified independently on the continuum λ² ∈ [20, 30] (log half-width 1.498–1.701; with WeilSimple λ² ∈ [2, 30]) and at the windows λ² = 35, 40 (1.778, 1.844), with explicit truncation thresholds at 25 windows, every process ≤ 180 MB; nothing beyond Botha's L ≤ 2 (λ² ≤ 54.6): an independent audit of the first stretch of his imported-trust Piece IV, no new window; λ² = 47, 54 not reached (time, not memory).

# CCMFMM — CCM's simple even ground state beyond λ² = 20, in the particle form (ccmfmm, 11 Nov 2026)

RH relevance: N13 (and N1).

## 1. Statement
Notation: c = λ², half-width a = log λ = (ln c)/2 (Botha's and TREE's "L"), Q = QW_λ = the Weil form of ../WeilZeta on L²(−a, a)
(../WeilSimple/DERIVATION.md sec. 1); λ_k^even, λ_k^odd = min-max values of Q on even, odd functions.

**Theorem (computer-assisted: python-flint 0.9.0 Arb, FLINT 3 arb in C, the proved node scan of ../WeilZeta/RIGOR.md).**
(i) For every c in [20, 30] the operator A_λ of QW_λ (CCM sec. 3.3) has a simple smallest eigenvalue μ_λ, with an even eigenvector:
spec(A_λ) meets (−∞, L(c)) only in μ_λ ≤ U(c) (tables; between chain points U of the left point and L of the right point). With
../WeilSimple's chain: every c in [2, 30], i.e. λ in [√2, √30] = [1.414, 5.477]. The same holds at the single windows c = 35 and c = 40 (a = 1.778, 1.844).
(ii) At every window of the tables (integer c = 21, …, 30, 35, 40 and the 13 further chain points) QW_λ^N is even-simple (CCM Def. 5.3) for
every N ≥ NG (200 for c ≤ 25, 240 for 25 < c ≤ 30, 320 at c = 35, 400 at c = 40), so CCM's Theorem 5.10 holds there for all those N (explicit N₀, ../WeilSimple sec. 3).
(iii) Not covered: c in (30, 54.6] other than 35, 40 (Botha's claim, unaudited by us), c > 54.6 (beyond every claim), the second missing step (λ → ∞).

## 2. Method (parents: [WeilSimple](../WeilSimple/REPORT.md), [ExploreFMM](../ExploreFMM/REPORT.md); proofs: DERIVATION.md)
Three certified inequalities per window, as in WeilSimple: (a) U = Q(g)/‖g‖² ≥ λ₁^even (gstate.py = cert_simple.py's ground_vec + rayleigh
on the closed-form even block of weil_matrix.py, copied unchanged); (b) Q ≥ L_odd on odd f; (c) Q − L_2 + ⟨·, g⟩² ≥ 0 on even f, so λ₂^even ≥ L_2.
U < min(L_odd, L_2) gives (i) at c. (b) and (c) are decided in explore28's particle form instead of WeilSimple's dense Woodbury–Cholesky:
- (b): ExploreFMM's odd certificate unchanged (prep_ccm.py = prep_slt.py; slt.c point mode: inertia of J″(ν) = G + diag(δ), G the
  1/(y − y′) interaction of particles at y = k², 4j², y_v with two charges, plus the backward-error bound of Lemma B).
- (c), new: the deflation column g = Σ_{l ≤ NG} g_l φ_l is one dense U index (a "border") of J″, with r_b = ⟨g, P0⁻¹c_b⟩ = β_b Φ(y_b),
  Φ(y) = Σ_l g_l A_l/(p_l(y − 4l²)) (Lemma P′). It is eliminated first; its Schur complement is Cauchy-like of displacement rank 4
  (second generator pair u₂ = r, w₂ = y r/J″_gg, Lemma R), and the GKO recursion, pivoting, balancing and Lemma B run on both pairs
  (Lemma B′: border input radii, σ includes δ_g = 1/M, ‖Y_s‖ includes the g row). Memory 7n numbers.
- A level is certified when some ν > κη gives inertia (n_U + [1], n_B) (then λ_min(Id − W) ≥ ν − κη > 0); each run walks a ladder ν = ε·10^k.
- Levels: L_2 = (e₁e₂)^{1/2}, L_odd = (e₁o₁)^{1/2} from approximate eigenvalues of the truncated blocks (WeilSimple's rule; e₁ = U);
  chain: L(c_{i+1}) = 2.00001 U(c_i) for both sectors; c = 35, 40: L_2 = L_odd = (U o₁)^{1/2}, o₁ by inverse iteration on the odd block at
  N = NG = 320, 400 (o1fix.py; no full eigendecomposition).
- Scan premise at every window: rscan.c (proved rounding bound, RIGOR.md; flagged nodes re-checked in Arb), K2 = max(401, ⌈5c²/3⌉) as in
  WeilSimple; c = 35: ../WeilZeta/scan_zeta.py with the early-exit rscan2.c (RIGOR.md, last section): 977,745,180 nodes, 26 s, 65 MB,
  E = 4.4e-9 < margin/10, selected 664 odd and 1531 even nodes beyond K2 = 2043 (logs/scan35.out); c = 40: 3,325,526,502 nodes, 52 s, 72 MB, E = 1.6e-8,
  selected 1393 odd and 2254 even beyond K2 = 2667 (logs/scan40.out). The scan caches are in logs/scan/.

## 3. Validation: WeilSimple's L_2 and L_odd at c = 10, 16, 20 (logs/val/, val_ccm.sh)
Same NG, bits, K2, ε and the same g (same seed: U equals WeilSimple's to all 12 printed digits, 1.94449796975e-43, 3.74299185237e-75,
1.65743315856e-96). m = WeilSimple's margin (its λ_min(Id − W_up)); point mode at ν = m(1 − 10⁻⁹) must pass and at m(1 + 10⁻⁹) must fail:
| c | sector | ε (WeilSimple) | WeilSimple margin m | ν = m(1 − 1e-9) | ν = m(1 + 1e-9) | κη | n | RSS MB (prep, LDL) | s (prep, LDL) |
|---|---|---|---|---|---|---|---|---|---|
| 10 | defl | 7.348097e-40 | 9.55135345279e-37 | pass | fails | 2.1e-79 | 166 | 43, 3 | 0.8, 0.3 |
| 10 | odd | 1.349542e-41 | 7.19952806309e-40 | pass | fails | 1.2e-78 | 164 | 43, 3 | 0.9, 0.1 |
| 16 | defl | 4.085821e-71 | 1.12391899385e-67 | pass | fails | 8.5e-137 | 414 | 62, 4 | 2.3, 0.8 |
| 16 | odd | 4.580865e-73 | 2.93457950161e-71 | pass | fails | 4.4e-130 | 412 | 44, 3 | 1.7, 0.7 |
| 20 | defl | 2.841006e-92 | 1.62050936675e-88 | pass | fails | 5.2e-153 | 1101 | 63, 5 | 6.2, 7.6 |
| 20 | odd | 2.497138e-94 | 6.60901632398e-92 | pass | fails | 9.9e-145 | 1099 | 47, 5 | 3.2, 5.1 |
All six: λ_min(Id − W) = m(1 ± 10⁻⁹) ± κη, i.e. WeilSimple's margins to 9 digits, from independent linear algebra (dense Woodbury +
Cholesky there, structured inertia here). Regression: without a border slt.c reproduces explore28's c = 10 even bracket
[5.026355097, 5.026355107]e-45 (WeilZeta: 5.02635510181e-45). Cost at c = 20: WeilSimple 130 s + 116 s, 676 MB; here 14 s + 8 s, 64 MB.

## 4. Certificates at the integer windows (logs/int.jsonl, run_ccm.py)
| c = λ² | λ | half-width log λ | U ≥ λ₁ᵉᵛᵉⁿ | L_odd | L_2 | U < min | NG | bits | n (defl / odd) | certified ν (defl, odd) | max κη | peak RSS MB (g, prep, LDL) | s (g + certificates) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 21 | 4.5826 | 1.5223 | 7.024e-102 | 1.619e-99 | 1.195e-97 | yes | 200 | 640 | 1184 / 1183 | 1.20e-94, 1.62e-97 | 2.6e-162 | 63, 64, 6 | 96 |
| 22 | 4.6904 | 1.5455 | 2.824e-107 | 1.297e-104 | 5.826e-103 | yes | 200 | 640 | 1275 / 1274 | 5.83e-100, 1.30e-103 | 1.7e-165 | 63, 64, 6 | 94 |
| 23 | 4.7958 | 1.5677 | 1.426e-112 | 1.564e-109 | 3.038e-108 | yes | 200 | 704 | 1081 / 1080 | 3.04e-105, 1.56e-108 | 2.5e-182 | 69, 63, 8 | 119 |
| 24 | 4.8990 | 1.5890 | 5.579e-118 | 2.096e-114 | 1.383e-113 | yes | 200 | 704 | 1697 / 1696 | 1.38e-110, 2.10e-114 | 1.9e-177 | 54, 64, 8 | 226 |
| 25 | 5.0000 | 1.6094 | 2.984e-123 | 5.471e-121 | 7.277e-119 | yes | 200 | 768 | 1616 / 1615 | 7.28e-116, 5.47e-119 | 1.4e-199 | 70, 64, 8 | 168 |
| 26 | 5.0990 | 1.6290 | 9.599e-129 | 1.966e-126 | 2.861e-124 | yes | 240 | 768 | 2237 / 2236 | 2.86e-122, 1.97e-124 | 7.0e-200 | 87, 65, 9 | 241 |
| 27 | 5.1962 | 1.6479 | 4.544e-134 | 9.024e-132 | 1.316e-129 | yes | 240 | 832 | 1975 / 1974 | 1.32e-127, 9.02e-130 | 1.6e-217 | 86, 64, 10 | 429 |
| 28 | 5.2915 | 1.6661 | 1.727e-139 | 3.767e-137 | 5.967e-135 | yes | 240 | 832 | 2467 / 2466 | 5.97e-133, 3.77e-135 | 6.4e-217 | 75, 57, 11 | 385 |
| 29 | 5.3852 | 1.6836 | 8.691e-145 | 1.908e-142 | 3.067e-140 | yes | 240 | 832 | 2289 / 2288 | 3.07e-138, 1.91e-140 | 9.8e-218 | 91, 72, 16 | 423 |
| 30 | 5.4772 | 1.7006 | 3.274e-150 | 7.493e-148 | 1.244e-145 | yes | 240 | 896 | 4004 / 4003 | 1.24e-141, 7.49e-146 | 1.3e-232 | 97, 81, 16 | 792 |
| 35 | 5.9161 | 1.7777 | 3.355e-177 | 8.973e-175 | 8.973e-175 | yes | 320 | 1024 | 5500 / 5499 | 8.97e-174, 8.97e-174 | 3.6e-269 | 120, 84, 21 | 646 |
| 40 | 6.3246 | 1.8444 | 2.898e-204 | 8.944e-202 | 8.944e-202 | yes | 400 | 1152 | 8276 / 8275 | 8.94e-201, 8.94e-201 | 2.0e-269 | 174, 93, 37 | 1783 |
- "certified ν": the largest ladder value ν = ε·10^k with inertia (n_U [+1], n_B) and ν > κη, so λ_min(Id − W) ≥ ν − κη > 0; in rows 21–29
  the next ladder value fails in both sectors (λ_min bracketed within 10× for c ≤ 25, 100× after); at c = 30 (defl) and 35 the top passes. Peak RSS: gstate (Python), prep (Python), LDL (C), /usr/bin/time -l.
- L_odd at c = 21–24 came from an unconverged odd block (N = 120: o₁ 2× to 400× too high); those certificates passed anyway (a
  certificate at a higher level is stronger). At c = 25 the same slip put L_odd above λ₁^odd and the odd run failed (row kept in
  logs/int_failed.jsonl); with o₁ at N = 200 (o1fix.py, logs/o1.jsonl) it passes (fix25.py). From c = 26 on o₁ uses N = NG.
- Converged ratios (c = 25–30, approximate spectra): o₁/e₁ = 3.4–5.2·10⁴, e₂/e₁ = 6·10⁸–1.5·10⁹; certified L_odd/U ≈ 180–230,
  L_2/U ≈ 2.4–3.8·10⁴. Node sets are Chekmachev's (c = 29: 2289 = 1205 + 1084 particles, as in ExploreFMM, WeilZeta and his manifest).
- c = 35: o₁ = 2.40e-172 (o₁/U = 7.2·10⁴); c = 40: o₁ = 2.76e-199 (o₁/U = 9.5·10⁴); both certificates pass at the top of the
  ladder ν = 10ε at both windows (margin > 10ε).
- Consistency with N1: U(29) = 8.69e-145 ≥ Chekmachev's certified even lower bound 3.389e-145 (Zenodo 23848597); U(25) = 2.98e-123 ≥ 2^-416.

## 5. Chain (logs/chain.jsonl, chain_ccm.py)
At c_{i+1} both certificates at ε = L(c_{i+1}) = 2.00001 U(c_i); U(c_{i+1}) from gstate.py at c_{i+1}; steps from the interpolated
approximate odd eigenvalues (first step attempt 62/3 → 257/12 failed: its heuristic still used the unconverged N = 120 o₁; halved).
| c_i → c_{i+1} | U(c_i) | L(c_{i+1}) = ε | U(c_{i+1}) | certified ν (defl, odd) | n (defl / odd) | NG | bits | ok | s |
|---|---|---|---|---|---|---|---|---|---|
| 20 → 62/3 | 1.66e-96 | 3.31e-96 | 3.87e-100 | 3.3e-94, 3.3e-96 | 1018 / 1017 | 200 | 640 | ok | 47 |
| 62/3 → 257/12 | 3.87e-100 | 7.74e-100 | 4.01e-104 | 7.7e-98, - | 1093 / 1092 | 200 | 640 | FAIL | 44 |
| 62/3 → 511/24 | 3.87e-100 | 7.74e-100 | 1.85e-103 | 7.7e-98, 7.7e-100 | 1105 / 1104 | 200 | 640 | ok | 56 |
| 511/24 → 263/12 | 1.85e-103 | 3.71e-103 | 7.87e-107 | 3.7e-101, 3.7e-103 | 1121 / 1120 | 200 | 640 | ok | 59 |
| 263/12 → 271/12 | 7.87e-107 | 1.57e-106 | 2.22e-110 | 1.6e-104, 1.6e-106 | 1161 / 1160 | 200 | 704 | ok | 56 |
| 271/12 → 93/4 | 2.22e-110 | 4.45e-110 | 5.77e-114 | 4.4e-108, 4.4e-110 | 1505 / 1504 | 200 | 704 | ok | 82 |
| 93/4 → 287/12 | 5.77e-114 | 1.15e-113 | 1.55e-117 | 1.2e-111, 1.2e-113 | 1705 / 1704 | 200 | 704 | ok | 106 |
| 287/12 → 295/12 | 1.55e-117 | 3.11e-117 | 4.34e-121 | 3.1e-115, 3.1e-117 | 1584 / 1583 | 200 | 768 | ok | 98 |
| 295/12 → 25 | 4.34e-121 | 8.67e-121 | 2.98e-123 | 8.7e-119, 8.7e-121 | 1616 / 1615 | 200 | 768 | ok | 58 |
| 25 → 77/3 | 2.98e-123 | 5.97e-123 | 6.27e-127 | 6.0e-122, 6.0e-124 | 1894 / 1893 | 240 | 768 | ok | 103 |
| 77/3 → 79/3 | 6.27e-127 | 1.25e-126 | 1.70e-130 | 1.3e-125, 1.3e-127 | 1902 / 1901 | 240 | 832 | ok | 182 |
| 79/3 → 27 | 1.70e-130 | 3.40e-130 | 4.54e-134 | 3.4e-129, 3.4e-131 | 1975 / 1974 | 240 | 832 | ok | 96 |
| 27 → 83/3 | 4.54e-134 | 9.09e-134 | 1.18e-137 | 9.1e-133, 9.1e-135 | 2253 / 2252 | 240 | 832 | ok | 329 |
| 83/3 → 227/8 | 1.18e-137 | 2.36e-137 | 1.82e-141 | 2.4e-136, 2.4e-138 | 2327 / 2326 | 240 | 832 | ok | 256 |
| 227/8 → 349/12 | 1.82e-141 | 3.64e-141 | 2.80e-145 | 3.6e-140, 3.6e-142 | 3171 / 3170 | 240 | 896 | ok | 405 |
| 349/12 → 715/24 | 2.80e-145 | 5.59e-145 | 4.51e-149 | 5.6e-144, 5.6e-146 | 3241 / 3240 | 240 | 896 | ok | 316 |
| 715/24 → 30 | 4.51e-149 | 9.03e-149 | 3.27e-150 | 9.0e-148, 9.0e-148 | 4004 / 4003 | 240 | 896 | ok | 306 |

chain check: 16 ok records, 1 failed attempts; covered: [20, 30]; max RSS 88 MB

## 6. Adversarial test (logs/adv/; explore27's moved pair, as in ExploreFMM sec. 7)
γ₁, γ₂ → (γ₁ + γ₂)/2 ± 0.3i (zero count kept; prep_ccm.py adv = 1 adds explore27's exact Δμ(r) to μ in Arb), at c = 10 (WeilSimple's
levels) and c = 21 (ours), both sectors; the defl run keeps the true-Q vector g. Cells: inertia at ν (number of positive pivots).
| c | sector | ε | ν = −1 | ν = −0.1 | ν = 0 | ν = ε | n_U (expected #pos) | verdict |
|---|---|---|---|---|---|---|---|---|
| 10 | defl | 7.348097e-40 | ok (132 pos) | wrong (133 pos) | wrong (133 pos) | wrong (133 pos) | 132 | not certified (correct) |
| 10 | odd | 1.349542e-41 | ok (130 pos) | wrong (131 pos) | wrong (131 pos) | wrong (131 pos) | 130 | not certified (correct) |
| 21 | defl | 1.195093e-97 | wrong (881 pos) | wrong (881 pos) | wrong (881 pos) | wrong (881 pos) | 880 | not certified (correct) |
| 21 | odd | 1.618853e-99 | wrong (879 pos) | wrong (879 pos) | wrong (879 pos) | wrong (879 pos) | 878 | not certified (correct) |
Nothing is certified: at ν = 0 and ν = ε every run has one positive pivot too many, so λ_min(Id − W_moved) < κη; at c = 10 the inertia is
right at ν = −1 (λ_min ∈ (−1, −0.1)), at c = 21 not even there. The rank-one deflation by Q's ground vector cannot remove the moved
pair's negative direction, and the linear algebra adds no false positive.

## 7. Cost (one core, background priority on a shared machine: wall time 2–7× CPU time in the busiest runs)
| c | WeilSimple (dense Woodbury–Cholesky) | here (particle form + border): peak RSS (gstate / prep / LDL), wall time per window |
|---|---|---|
| 20 | 315 s (certificates 246 s), 676 MB | 64 / 63 / 6 MB; 6 s gstate + 14 s defl + 8 s odd (logs/val) |
| 25 | its limit, ≈ 1.5 GB (WeilZeta's dense c = 25: 1.08 GB per sector) | 70 / 64 / 8 MB; 168 s |
| 29 | does not fit 1.5 GB (WeilZeta's dense c = 29 stopped at 1.54 GB) | 91 / 72 / 16 MB; 423 s |
| 30 | does not fit | 97 / 81 / 16 MB; 792 s (194 s gstate incl. two eigendecompositions) |
| 35 | does not fit | 120 / 84 / 21 MB; 646 s (scan 26 s, gstate 70 s, o1fix 54 s, defl 306 s, odd 215 s) |
| 40 | does not fit | 174 / 93 / 37 MB; 1783 s (scan 52 s, gstate 111 s, o1fix 104 s, defl 862 s, odd 706 s); n = 8276, 87,891 sources |
Beyond: c = 47 per sector ≈ 40–60 min (ExploreFMM's measured c = 47 runs: 371,000 sources × 13,234 particles, O(n²) LDL) plus gstate at
NG ≈ 400 and the border (≈ 1.5× LDL); c > 54.6 (new) needs ExploreFMM's c = 53 estimate (≈ 10¹⁰ Cauchy terms, n ≈ 2·10⁴, a 10¹¹-node
scan with rscan2): hours per sector, all under 300 MB. The wall is time, not memory.

## 8. Novelty (2026-11-11)
- python3 .scratch/lst.py 'prolate|ground state|Weil form|Connes' --zenodo '"ground state" Weil': 5 arXiv title hits (2026-10, 2026-11):
  ACLG 2611.13129 (numerics), 2610.32082 (multiband prolate matrices), three off-topic; Zenodo: 24076153 (positivity through (log 11)/2,
  not simplicity), 24086102 (Thorin computations), rest off-topic. A scout's Zenodo/arXiv API queries ("Weil quadratic form",
  "Connes Consani Moscovici", "zeta spectral triples", "even ground state", "ground state" "Weil form"): nothing else certifies simplicity.
- Botha, Zenodo 23852617 v2.0.0 (1 Nov 2026; still the latest version): claims every a = log λ in (0, 2], i.e. c ≤ e⁴ = 54.6. All our windows
  lie inside: NOT new mathematics. By his Limitations (R415) Piece III (a in (1.19, 1.4035]) is unaudited and Piece IV (a in (1.4035, 2])
  is "conditional on the imported trust" of rStation AI-agent preprints. Ours (with WeilSimple) independently certifies a in [0.347, 1.701] and a = 1.778, 1.844,
  i.e. all of his Piece III, the first stretch (1.4035, 1.701] of Piece IV and two windows inside it, with independent code (closed-form Galerkin U,
  Chekmachev-type minorant, structured inertia).
- Kramarenko-Byrd 23007100 (a ≤ 4/5), Zhu 2608.24827 (one window, a = 0.8), Siche (April 2026, all windows, unverified): below or not
  checkable. ACLG 2611.13129: simple even ground state for 5 ≤ λ² ≤ 300 from truncated Arb matrices, truncation error not controlled
  (numerics only); our certificates confirm it rigorously for 5 ≤ λ² ≤ 30 and λ² = 35, 40.
- Our files: WeilSimple ([2, 20], the method), ExploreFMM (particle form; c = 47 for N1), NOVELTY of WeilSimple/ExploreFMM.
- Verdict: an audit result (independent confirmation of a published, partly unaudited claim on a new range) plus engineering (the
  deflation as a displacement-rank-4 border); no discovery. The explicit truncation thresholds N₀ = NG at the new windows extend
  WeilSimple's list (RECORDS-type).

## 9. Limits and trust
- Trust: Arb (python-flint 0.9.0; FLINT 3 via Homebrew in slt.c); IEEE binary64 +, −, × in rscan.c / rscan2.c (RIGOR.md); our Python
  and C, unchecked by a proof assistant. g, the pivot order, shears and the ν ladder are floating choices that enter only as exact inputs.
- weil_matrix.py's closed forms are WeilSimple's (checked there by quadrature to 3.5e-49); the minorant is Chekmachev's (audited SOUND).
- Point-mode η has a floor from double underflow (ExploreFMM): a valid overestimate. table_ccm.py re-checks all 206 point-mode LDL lines of
  the 61 certificate runs (logs/runs): σ > 2η (so S_s − η > 0, as Lemma B needs), κ finite; max η/σ = 5.7e-175.
- Reviewer (verdict on DERIVATION.md + the ldl() border code): Lemma R's sign and the two-pair recursion exact, B′'s border bounds, σ and
  ‖Y_s‖ correct; GAP: slt.c's double-precision lower bound of |y_a − y_b| for a non-integer point can fail for two close huge
  non-integers. Closed by the point domain (DERIVATION.md, after Lemma B′): the only non-integer point is the pole y_v ∈ (−0.35, −0.22) (c ∈ [20, 40]),
  all others are integers ≥ 0, so every such distance is ≥ 0.22 and its binary64 evaluation errs by ≤ 3·2⁻⁵³ relative, inside the 10⁻¹²
  margin. The two-pair updates (slt.c lines 168, 172) and κ (unchanged from ExploreFMM) were not shown to the reviewer; checked by us.
- Between chain points the truncation statement holds for N ≥ some non-explicit N₀(c) only (as in WeilSimple).
- Nothing here bears on λ → ∞ (the second missing step); μ_λ ≥ 0 for all λ is RH.

## 10. Files and replay
New (ccmfmm): gstate.py (g, U, approximate eigenvalues); prep_ccm.py (= ../ExploreFMM/prep_slt.py + the defl border, lines marked CCM);
slt.c (= ../ExploreFMM/slt.c md5 63d228c2f525a58a4592b936cb60fb36 + border and second generator pair, lines marked CCM); run_ccm.py
(integer windows), chain_ccm.py (chain), o1fix.py, fix25.py; one.sh, val_ccm.sh, q3.sh, q4.sh, q5.sh, q6.sh (the exact run sequences); table_ccm.py
(tables and the re-checks of every chain inequality and every LDL line); DERIVATION.md. Copied unchanged: weil_matrix.py (../WeilSimple),
rscan.c (../WeilZeta, md5 bd77d55f0a0bff7e7955e6c734aab9a9), prep_slt.py (../ExploreFMM, for the diff). logs/: int.jsonl, int_failed.jsonl, chain.jsonl, o1.jsonl, val/, adv/,
runs/ (every prep and LDL JSON line and /usr/bin/time record), scan/ (rigorous scan caches: selected node sets), scan35.out/.time, scan40.out/.time.
Replay (from this folder; Homebrew flint; Python only via uv; one process at a time):
    cc -O2 -o ../../../.scratch/ccmfmm/slt slt.c -lflint -lgmp -lmpfr -lm -I/opt/homebrew/include -L/opt/homebrew/lib
    cc -O2 -ffp-contract=off -o ../../../.scratch/ccmfmm/rscan rscan.c
    ./val_ccm.sh                       # c = 10, 16, 20: six certificates, about 1 min, < 70 MB
    uv run --no-project --with mpmath --with numpy --with python-flint python3 run_ccm.py /tmp/int.jsonl 21     # about 2 min, < 70 MB
    CCM_U0=1.65743315856e-96 CCM_O0=3.76226271626722e-92 uv run --no-project --with mpmath --with numpy --with python-flint python3 chain_ccm.py logs/int.jsonl /tmp/ch.jsonl 20 21
    python3 table_ccm.py               # tables and checks from logs/

## Suggested TREE line (N13)
- … independently confirmed for λ ∈ [√2, √30] (weilsimple [2, 20]; ccmfmm [20, 30]: weilsimple's deflation in explore28's particle form as a
  displacement-rank-4 border, every process ≤ 180 MB, explicit N₀ at 25 more windows) and at λ² = 35, 40: all of Botha's Piece III and
  (1.4035, 1.701] of his imported-trust Piece IV [CCMFMM](Riemann/Deep/CCMFMM/REPORT.md).
