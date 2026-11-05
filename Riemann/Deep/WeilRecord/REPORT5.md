CERTIFIED: half-width 1.6 (independent certificate; not a record: Chekmachev, Zenodo 23848597, reaches 1.925, audited sound in ChekAudit)

# REPORT5 -- Weil positivity on [-1.6, 1.6] with a memory-lean pipeline (weilooc, 2026-11-05)

Node N1 of the RH approach tree (TREE.md): positivity of Weil's quadratic form Q on a finite window.

## Result
For every complex C^2 function f supported in [-1.6, 1.6] (window length 3.2, x = e^3.2 = 24.5, prime powers 2..23):
  Q(f) >= 2^-430 ||f||^2  (2^-430 = 3.6e-130),
both parities, computer-assisted (exact integer arithmetic + Arb balls), proof chain of PROOF.md sections 1-5 unchanged;
the new pipeline's steps are proved in ooc/DERIVATION.md.  Certified margins after every quadrature, truncation and
rounding error: lambda_min(W^-1 + Z) >= 3.481e-124 (even), >= 1.148e-119 (odd).
A function supported in [-L', L'] with L' < 1.6 is supported in [-1.6, 1.6] and Q is the same functional, so this also
covers 1.56 and 1.58; they were not run separately.

| L | primes | norm of C_p <= | gamma | Omega | u | N (N/c) | F (bits) | tau | margin even | margin odd |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.55 (re-run) | 2..19 | 4.591281212 | 1395/256 | 1640 | (913, 303) | 3400 (1.338) | 448 | 2^-372 = 1.04e-112 | 1.648e-111 | 2.717e-107 |
| 1.6 | 2..23 | 4.903480029 | 1475/256 | 2242 | (1300, 432) | 4700 (1.310) | 480 | 2^-430 = 3.6e-130 | 3.481e-124 | 1.148e-119 |

1.6 (even / odd): smallest eigenvalue of the computed H' = m^2 (W + W Z W) ~ 8.38e-123 / 2.97e-118, sigma' = 2^-407 / 2^-392,
||H' - sigma' I - L L^T||_F <= 2.1e-142 (exact residual), ||H'_true - H'~||_2 <= 2.8e-132 / 2.9e-132,
epsJ (band, poles, U) 7.4e-136, ||W'||_2 <= 2.948 / 2.939.  H1: weighted Schur test on all 32000 cells (weight in
ooc/runs/weights).  Full numbers: ooc/runs/r160b/state.json and log.txt.

## Consistency checks
- 1.55 re-certified (parameters of REPORT3, with gamma, s, tau rounded to dyadics): margins 1.648e-111 / 2.717e-107
  (cert4: 1.215e-111 / 4.004e-107); the matrix H is the same: lambda_min(H~')/m^2 = 4.60e-110 / 1.16e-105
  against cert4's lambda_min(H_mid) = 4.598e-110 / 1.162e-105.
- mu = 5: the same H as cert4 (lambda_min 6.705e-20 / 4.825e-16 in both codes); mu = 16 (log 4): certified 2.541e-78 / 4.311e-74.
- Gram matrices (gvals vs grams_c, mu = 5 N = 140 and L = 1.55 N = 800) and band matrices (bvals vs cert4.band_lean,
  mu = 5) agree within the old balls' radii; fixed-point products checked bit-exact against Python integers.

## The memory-lean pipeline (ooc/, README.md)
- Every matrix is a fixed-point disk file (F fractional bits); products are exact (45-48 primes < 2^21, float64 BLAS,
  explicit CRT), one rounding per stored entry, rigorous 2-norm error bounds carried in Arb.  Blocks of 480 x 480,
  primes in 3 groups, files read by seek/read (no memmap).
- Legendre values by a norm-one rotation recurrence (D_j = P_j - e^{i theta} P_{j-1}): F + 64 working bits instead of
  the ball recurrence's F + 1.28 N + 64 bits.
- Division-free W' = m W(X), H' = m^2 H with dyadic m; X by iterative refinement (5-7 exact products); block Cholesky
  with an exactly computed residual.
- One process at a time, single-threaded, fresh process per stage, time budget 2100 s inside perl alarm 2400, restart
  when RSS > 600 MB; every stage resumes at block level.
- Resources (shared machine, load 6-27): 1.55: ~48 min wall, about 40 min CPU (old: ~2.8 h CPU, 4-5 GB), max RSS 563 MB
  after the RSS guard (one early process reached 1.1 GB before it existed; killed, resumed).  1.6: 3.5 h wall on the
  shared machine (load 6-45), 1.8 h user CPU + 0.2 h system in 24 processes, max RSS 566 MB, max footprint 632 MB.
  Scratch disk peak ~10 GB (value files), deleted.

## What failed (and why)
- 1.6 at F = 448, tau = 2^-400 (ooc/runs/r160_F448): even lambda_min(H~') ~ 5.5e-124 but the error bound e_H' was
  1.2e-122, and tau = 2^-400 had absorbed most of the margin.  Precision-limited, not a positivity failure: F = 480,
  tau = 2^-430 certified.
- The margin falls fast: lambda_min(H') ~ 3.4e-110 at 1.55 and ~8.4e-123 at 1.6 (even), about 12 decades per 0.05
  in L; the precision needed grows by about 40 bits per 0.05.

## Cost estimate for the next steps (this pipeline, one core)
- log 5 = 1.6094 (same prime set; Schur 4.960357039, gamma 1490/256, Omega 2377, u (1385, 461), N 5012): ~1.2 x the
  1.6 run (~2.2 h CPU, 3-4 h wall here), F = 480, extrapolated even margin ~1e-125 (errors ~3e-132).
  Ready: make_config.py 25 w_25_32000.txt <dir> 1.31 440; chain.sh <dir>.
- 1.7 (25, 27, 29 enter; Schur 5.564634257, Omega 4355, N 9698, d 4849): ~9 x the 1.6 cost (about 16-20 h CPU), margin
  ~1e-147 if the decay rate holds, so F ~ 560; ~18 GB of scratch value files if the parities' values are generated
  one after the other.

## Files
- ooc/: fxmm.py, oocla.py, oocert.py, gvals.c, bvals.c, bandplan.py, fxkern.c/.py, wgen.py, make_config.py,
  chain.sh, README.md, DERIVATION.md; runs/ (config, state, log of mu = 5, log 4, 1.55, 1.6 at F = 448 and 480);
  runs/weights (Schur weights for 1.55 and 1.6, 32000 cells).  Existing certificate files untouched.

## Status
External computer-assisted certificate, not refereed, not in Lean, not progress toward RH itself.  Trust: FLINT/Arb and
GMP, Accelerate dgemm on integers < 2^51 (exact in IEEE double), this code, and the lemmas of PROOF.md and DERIVATION.md
(new: rotation-recurrence error bound, fixed-point Gram error budget, mirrored-symmetric error accounting).
Novelty (light check): arXiv math.NT and math.CA listings for 2026-10 and 2026-11 show no certified window beyond 1.55
(2610.32704 assumes RH; 2611.02944 is about Polya frequency of Xi); R284 found no newer certified window.
