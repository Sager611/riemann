NO RECORD: a deflated certificate (high precision only on the k ~ c near-null directions, all O(n^3) work in binary64) re-proves Q(f) >= eps ||f||^2 for all f at c = 16, 20, 25, 29 and 31 with the same margins as our full-precision pipeline and Chekmachev's; at equal precision it takes the same time (c = 16: 7-8 s at 320 bits for both; Chekmachev's run: order 133 at 512 bits) but 2-10x less memory, which brings c = 29 and 31 (new for our code) under 0.6 GB; the tail-integral ladder of explore17 cannot carry a full-space proof; the gain is a constant factor (c = 47 estimated at 10-12 h and ~1.3 GB, against about a day and 35 GB), and the node scan (2.5e13 nodes at c = 71, 1.2e16 at c = 100) blocks the record range.

# WeilDefl: deflating the near-null ladder in Weil window certificates (weildefl, 2026-11-06)

RH relevance: N1 (Weil positivity on windows; cost of window certificates).

Window and form exactly as TREE.md N1 (Chekmachev's convention, also ../WeilZeta): W = [-(ln c)/2, (ln c)/2],
Q(f) = 2<f,cosh(x/2)>^2 - 2<f,sinh(x/2)>^2 + (1/2 pi) int (Re psi(1/4 + it/2) - log pi)|F(t)|^2 dt - 2 sum_{n<c} Lambda(n) n^{-1/2} Re (f * f~)(log n), F(t) = int f e^{itx}.

## Claim
1. Theorem (computer-assisted: Arb balls, the proved node scan of ../WeilZeta/RIGOR.md, binary64 with a proved rounding bound). For every f in
   L^2(W) of the stated parity, Q(f) >= eps ||f||^2, with (c, parity, eps) as in the results table; for each c both parities hold, so the statement
   holds for all f with the smaller eps. c = 16: eps = 1.512800554e-75 for all f.
   Not new mathematically: Chekmachev (Zenodo 23848597) certifies every c <= 47. New for our own code: c = 29 and 31 (../WeilZeta stopped at 25,
   its c = 29 run passed 1.5 GB); our c = 29 margins agree with Chekmachev's published margin_lower values to his 7 digits, and at c = 31 with
   0.1 Delta from his eps_choice.csv (his rule eps = 0.9 Delta/slope leaves the margin 0.1 Delta; logs/chek_manifests.log). c = 31 is the end of the
   range his own Arb checker covers in full (Zenodo 23986313); our two c = 31 runs re-certify it independently (both parities, so all f with
   eps = 5.7945727624820630486e-156), 26 min and 0.58 GB each in lean mode.
2. Method: the k ~ c tiny eigenvalues of the certificate matrix are split off by a Schur/Cauchy-Schwarz lemma, so that only k columns need
   high precision; margins equal ../WeilZeta's to all 12 printed digits (c = 5, 10, 16, 20, 25). At equal precision it is as fast as the
   full-precision pipeline and needs 2-10x less memory (lean mode: no high-precision matrix is ever stored).
3. Negative: the brief's version (V = span{Phi^(j)|_W}, B and Q_V from the tail integrals, Q_{V-perp} by Galerkin plus a tail bound) cannot give a
   full-space certificate; it would need frequency cut-offs of about 1e38 at c = 16 (table below).

## Method (DERIVATION.md: every step and constant)
- Full-space reduction, unchanged from ../WeilZeta (cert_zeta.py imported as is): Chekmachev's lattice minorant with the rigorous node scan gives
  Q - eps ||.||^2 >= P - N on all of L^2(W), and P - N >= 0 <=> T = I - X + Y^T S^{-1} Y >= 0 (T of order m; S of order nU). The scan bound beyond
  kscan (h increasing, tail decreasing) is re-derived in DERIVATION.md (A1'); it plays the role of the brief's "high-frequency tail bound".
- Deflation lemma: for any m x k matrix Yt and delta > 0, T >= 0 if (i) T >= delta on Yt-perp and (ii) the k x k matrix
  Yt^T T Yt - (P_perp T Yt)^T (P_perp T Yt)/delta >= 0. The loss is quadratic in the residual P_perp T Yt.
- Yt: the invariant subspace of the k eigenvalues of T below 1e-2 from a binary64 eigensolve (the cluster below 1e-15 is noise there, its subspace is
  not: gap delta0 = 0.06-0.25), refined at p/2 + 64 bits by a correction equation solved in binary64 (12-15 digits per step).
- (i): T >= T_lo (Chekmachev Lemma 5.3), R_aug = T_lo + Yt Yt^T - delta I >= 0 by a shifted Cholesky in binary64 with a rigorous rounding-error bound
  (or, WD_COMPL=arb, in 96-bit Arb). (ii): in Arb at p bits on k columns: a lower bound of Yt^T T Yt that needs no S^{-1}, and an explicit bound
  for the S-solve residual (S >= D^{-1}). Lean mode (WD_LEAN=1): X, Y, S are never stored; every product re-evaluates the closed forms in row blocks.
- Bits: p = log2(1/margin) + 64 is still needed (c = 16: 320, the margin is 2^-253), but only for O(n^2 k) work; O(n^3) runs in 53 bits.

## Certificate (one command per c and parity)
  cd Riemann/Deep/WeilDefl && cc -O2 -ffp-contract=off -o /tmp/rscan ../WeilZeta/rscan.c && \
  WZ_RSCAN=/tmp/rscan WZ_SCAN=/tmp uv run --no-project --with python-flint --with mpmath --with numpy python3 cert.py 16 even 1.512800554e-75 320
  (odd: cert.py 16 odd 2.105612047e-71 320; c = 29: cert.py 29 even 3.3887641037680338350e-145 560 and cert.py 29 odd 1.5329466444386175069e-140 560;
  add WD_LEAN=1 for the lean mode, and WD_CKPT=<file> to keep the Cauchy sums across a restart). One JSON line; "rigorous": true iff (i) and
  (ii) pass. Logs: logs/final_validation.log, final_lean2.log, rerun_lean.log, c29_*.log, c31_*.log.

## Results and cost (one core, single-threaded, shared machine; wall clock; RSS = peak resident memory)
| c | parity | eps | margin (ours) | Chekmachev margin_lower | m / nU | k | delta | bits p / p_mid | ours: time, RSS | WeilZeta: bits, time, RSS |
|---|---|---|---|---|---|---|---|---|---|---|
| 16 | even | 1.512800554e-75 | 7.2936230421e-77 | 7.293616e-77 | 133 / 280 | 15 | 0.117 | 320 / 224 | 7.4 s, 69 MB (lean 11 s, 57 MB) | 512: 22 s, 117 MB; 320: 7.8 s, 120 MB |
| 16 | odd | 2.105612047e-71 | 2.99318593786e-72 | 2.993183e-72 | 133 / 279 | 14 | 0.039 | 320 / 224 | 10.5 s incl. a fresh 1.7 s scan, 84 MB (lean 14 s, 59 MB) | 512: 17 s, 144 MB; 320: 7.4 s, 120 MB |
| 20 | even | 2^-327 | 3.1057122801e-97 | (other eps) | 311 / 789 | 19 | 0.127 | 384 / 256 | lean 77 s, 96 MB | 448, 82 s, 0.70 GB |
| 20 | odd | 2^-327 | 6.70422678283e-92 | (other eps) | 311 / 788 | 18 | 0.039 | 384 / 256 | lean 88 s, 96 MB | 448, 81 s, 0.68 GB |
| 25 | even | 2^-416 | 5.21043141861e-124 | (other eps) | 558 / 1058 | 24 | 0.108 | 480 / 304 | lean 311 s, 0.15 GB (stored: 322 s, 0.70 GB) | 576, 295 s, 1.08 GB |
| 25 | odd | 4.0546561810437569560e-119 | 8.07878305192e-120 | 8.078775e-120 | 558 / 1057 | 23 | 0.030 | 480 / 304 | lean 284 s, 0.17 GB (stored: 314 s, 0.69 GB) | 576, 414 s, 1.08 GB |
| 29 | even | 3.3887641037680338350e-145 | 1.42611291013e-146 | 1.426111e-146 | 1084 / 1205 | 28 | 0.092 | 560 / 344 | lean 850 s, 0.27 GB (stored: 765 s, 0.87 GB) | stopped at 1.5 GB |
| 29 | odd | 1.5329466444386175069e-140 | 2.25792756276e-141 | 2.257925e-141 | 1084 / 1204 | 27 | 0.033 | 560 / 344 | lean 843 s, 0.26 GB (stored: 806 s, 0.87 GB) | not run |
| 31 | even | 5.7945727624820630486e-156 | 2.28614928017e-157 | 0.1 Delta = 2.286149e-157 (his eps_choice) | 1663 / 1628 | 30 | 0.115 | 592 / 360 | lean 1550 s, 0.58 GB | does not fit 1.5 GB |
| 31 | odd | 3.1251202725938626635e-151 | 6.38218067472e-152 | 0.1 Delta = 6.382181e-152 (his eps_choice) | 1663 / 1627 | 29 | 0.034 | 592 / 360 | lean 1532 s, 0.58 GB | does not fit 1.5 GB |
Margin = smallest eigenvalue of the k x k block Sk (= smallest eigenvalue of T to all digits shown); the certified lower bound is half of it (nu2).
Chekmachev's margin_lower (his manifests, Zenodo 23847003; a rigorous bound about 1e-6 below his computed margin) agrees wherever the eps agree.
All "lean" entries and the c = 16 stored entries are the final cert.py (logs/final_lean2.log, rerun_lean.log, c31_*.log, final_validation.log);
the "stored" entries at c = 25, 29 came from the version before the last speed patch (refinement tolerance and symmetric lean products only;
no bound changed), with identical margins. WeilZeta at 320 bits (logs/weilzeta_c16_320bits.log): same margins, 7.8 s / 7.4 s, 120 MB, so at
equal precision the time is the same and the gain is memory. Time split at c = 29 (lean, even): Cauchy sums over |E| = 48,670 (shared with
WeilZeta) 261 s, refinement 479 s (8 steps), k x k block 77 s, binary64 copies 12 s, complement 1 s.

## Why the explore17 tail-integral ladder fails in full space (DERIVATION.md sec. 7; obstruct.py; estimates, not theorems)
| c | rung edge value Phi^(L/2) | coupling L^2 mass above Omega = 1000 | its Schur share | lambda_min | Omega with share < lambda | minorant drop R(v,v) | part beyond kscan |
|---|---|---|---|---|---|---|---|
| 5 | 3.6e-4 | 6.2e-11 | 7.7e-12 | 1.2e-17 | 2.5e8 | 6.4e-10 | 6.4e-10 |
| 10 | 2.7e-10 | 9.3e-23 | 7.6e-24 | 1.9e-43 | 8.5e21 | 4.6e-22 | 1.6e-23 |
| 16 | 5.1e-18 | 5.1e-38 | 3.4e-39 | 3.7e-75 | 1.4e38 | 1.2e-37 | 4.7e-40 |
| 20 | 3.0e-23 | 2.3e-48 | 1.3e-49 | 3.7e-99 | 4.7e51 | 3.0e-48 | 9.5e-52 |
| 25 | 7.4e-30 | 1.6e-61 | 8.1e-63 | 5.9e-126 | 1.6e65 | 1.4e-61 | 1.8e-65 |
- The coupling beta = Q(Phi^(j)|_W, .) jumps by w_n Phi^(j)(L/2) at +-(L/2 - log n) and is log-singular at +-L/2, so its high-frequency part is
  ~ e^{-2 pi c}, far above lambda_min ~ e^{-10c}: it must be computed to relative accuracy ~ e^{-3.7c}, which no Galerkin cut-off plus bound does
  (Omega column). In the lattice minorant (the only full-space reduction at hand) the raw rung loses R(v, v) >> lambda_min, and the part beyond kscan
  is known only to O(1). Chekmachev's log-Laplacian paper (Zenodo 23978194) shows the true minimiser behaves like beta (ln 1/d)^{-1/2} at the edge with
  beta^2 = -(1/2) d lambda/da ~ e^{-10c}: its edge coefficient is ~ e^{-5c}, the rung's ~ e^{-pi c}; closing that gap is the high-precision step.
- Boundary-adapted rungs P(D)Phi|_W vanishing to order p at the edge (e1_adapted.py, logs/e1_adapted.log) become nearly dependent (Gram condition 6e17 at c = 10, p = 4) and
  still lose 5e-13 (p = 4) or 7e-23 (p = 8) against lambda_min = 1.9e-43 in the worst direction (a direction-by-direction analysis was not done).
- What survives of the lead: the deflated subspace is the ladder's. k = 4, 9, 15, 19, 24, 28 at c = 5, 10, 16, 20, 25, 29 (even; tau = 1e-2), close to
  explore17's k ~ c; at c = 16 the principal-angle cosines between span{B^T Phi^(2i)|_W, i < 19} and the deflated subspace of T are all >= 0.84 (the
  first ones 0.9997, 0.9990, 0.9983; logs/ladder.log). The certificate reads the subspace off T itself.

## Scaling: c = 31, 47, beyond (measured scaling plus estimates)
- Node scan: kscan = 2 L e^{2 sigma(c) + 1} exactly matches the code (sigma = sum_{n<c} Lambda(n)/sqrt n): 3.1e5 (c = 16), 8.8e7 (29), 5.1e10 (47:
  658 s with ../WeilZeta's early-exit scan), 2.8e11 (50), 2.5e12 (60), 2.5e13 (71: ~4 core-days), 1.2e16 (100: years). The lattice minorant (Chekmachev's
  method and ours) cannot reach c = 100.06 whatever the linear algebra.
- c = 47 for all f (an independent record-level check; Chekmachev's sizes m = 7170, nU = 6063, |E| = 371,000, margin 6.5e-244, p ~ 870),
  extrapolated from the measured lean runs (c = 29: 12.9 us per closed-form entry per product including the block product; 60 s per refinement step):
  Cauchy sums 4.9e9 terms (~4 h), binary64 copies (~20 min), ~10 refinement steps of ~1.5e8 entries each (~35 min each, of which ~20 min are the
  Arb n x n by n x k products, 1.1e10 multiply-adds at ~500 bits), final pass (~1 h): about 10-12 h on one core, to be split into checkpointed
  pieces under 30 min (one refinement step already exceeds that). With the closed forms and the Cauchy sums in C: ~4 h, bounded by the Arb
  products. Memory: O(n) high precision; the binary64 complement (orders 7170, 6063) needs ~1.3 GB with in-place kernels (the present compl_double
  holds ~15 matrices of order m: ~6 GB at c = 47), about half in binary32 with the same kind of rounding bound (gamma_n ~ 4e-4 << delta).
  Against WeilZeta's estimate (35 GB, about a day): ~2x in time, ~25x in memory. Not done here.
- Beyond 47 the bad-node count grows (|E| = 371,000 at 47), n and the O(n^2) binary64 complement grow with it; c ~ 50-55 is the limit of this route in
  1 GB, c = 71 needs days of scanning.
- Cost per step, measured: Arb multiplies n x n by n x k at ~100 ns per multiply-add (256 bits), closed forms ~5 us each; WeilZeta's dense
  O(n^3) stage at full precision is the part removed. High-precision flop count: ~n^2 k p/50 against ~n^3, a factor ~3-5 at c = 47.

## Checks
- Validation: margins equal ../WeilZeta/logs/val.log, z20.log, z25.log to all 12 printed digits at c = 5 (both), 10 (even), 16, 20, 25 (both)
  (logs/final_validation.log, final_lean2.log, rerun_lean.log); the Schur loss 2|P_perp T Yt|^2/delta is below 1e-34 of the margin in every run.
- Adversarial 1, one prime weight (logs/run_c16pert_c20.log): c = 16 even, w_2 -> (1 + 1e-70) w_2: margin -2.7486e-72, FAILS as it must;
  (1 - 1e-70): +2.7459e-72, passes.
- Adversarial 2, Davenport-Heilbronn (WD_FORM=dh; conductor 5, odd gamma factor, signed -f'/f coefficients, no pole; zeros off the line):
  c = 16, eps = 0: certified >= 0 with margins 2.55e-14 (even), 5.39e-11 (odd), below independent float64 Galerkin values 7.2e-14, 9.2e-11 (upper
  bounds), as they must be; with eps = 8e-14 > 7.2e-14: margin -3.96e-14, FAILS (logs/dh_c16.log, sharpness.log, dhscan.log). DH's window form stays
  positive (float64, to its noise floor) up to c = 30: its off-line zero 0.8085 + 85.70i is not seen by windows this small, so no DH window that must
  fail is within reach.
- Sharpness: c = 16 even with eps = 3.75e-75 > U(16) = 3.7430e-75 (an explicit Rayleigh quotient, ../WeilSimple): margin -8.98e-76, FAILS.
- Both implementations of step (i) pass (binary64 + bound: Etot = 4e-11 at c = 16, 1.5e-8 at c = 29, against s0 - nu >= 0.008); lean and stored
  modes give identical margins (c = 16, 20, 25, 29).
- Reviewer (DERIVATION.md sec. 1-5 and cert.py): SOUND under the stated arithmetic assumptions; one note (missing |kappa| in the general-kappa
  interface) fixed (logs/review1.md). Second reviewer (results table against the logs): all entries match except the c = 16 odd lean cost, which
  came from an earlier log; rerun with the final code and corrected (logs/review2.md, final_lean2.log).
- No float enters the proof without a bound: Yt, W, Ws, the Cholesky factors and s0 are arbitrary exact inputs of proved inequalities; the binary64
  complement carries the rounding bound of DERIVATION.md (i'); the scan carries RIGOR.md's bound.
- Constants of DERIVATION.md sec. 6 rg-checked in cert.py: logs/rg_constants.log.

## Caveat
- Rests on ../WeilZeta (minorant, closed forms, scan; audited in ../ChekAudit) and python-flint 0.9 Arb; step (i) by default also on IEEE binary64
  and a classical dgemm (the WD_COMPL=arb variant avoids that and passed at c = 5-25).
- Time equals WeilZeta's at equal precision (c = 16: 7.4/8.2 s against 7.8/7.4 s at 320 bits; its 22 s came from its default 512 bits); the
  gain is memory (lean mode) and the removal of O(n^3) work at full precision, which only starts to pay at c >~ 30; a constant factor, not a change
  of scaling.
- The obstruction for the tail-integral ladder is quantitative estimation, not a theorem.

## Novelty
arXiv math.NT, math.SP listings 2026-10, 2026-11 (show=2000; titles Weil|positivity|explicit formula|Connes|deflat|mixed precision|Schur): nothing on
deflated or mixed-precision Weil window certificates (nearest 2610.32704, 2611.02944, known). Zenodo ("truncated Weil form", "Weil positivity",
Chekmachev, Wiktorowicz): Chekmachev 23848597, 23986313 (verifier), 23847003 (data), 23978194 (log-Laplacian boundary balance, explains sec. 7; no
certificate method); Wiktorowicz 23869032 (LDL at 7300 bits). Connes-Consani-Moscovici 2511.22755, Connes 2602.04022: no deflation (explore17).
The lemma is standard linear algebra (Schur complement with Cauchy-Schwarz; invariant-subspace refinement): engineering, no DISCOVERIES entry.

## Suggested NOGO line (for the orchestrator; this agent edits only WeilDefl/)
- weildefl (Orca, N1 certificate cost; explore17's ladder lead) NO RECORD: deflating the k ~ c near-null directions of the lattice-minorant matrix
  (Schur/Cauchy-Schwarz lemma; O(n^3) in binary64, high precision on k columns) re-proves c = 16, 20, 25, 29, 31 for all f with WeilZeta's margins
  (c = 29, 31 for the first time with our code, under 0.6 GB; Chekmachev's margins to 7 digits); same time as full precision at equal bits,
  2-10x less memory; constant-factor gain only. The
  tail-integral Phi-ladder (B, Q_V from Phi^(j) 1_{|u|>L/2}) cannot give a full-space proof: the couplings jump by w_n Phi(L/2) at L/2 - log n, their
  high-frequency Schur share e^{-2 pi c} must be resolved to e^{-10c} (cut-off ~1e38 at c = 16); the minimiser's edge coefficient is e^{-5c}
  (Chekmachev 23978194) against the rung's e^{-pi c}. Scan kscan = 2L e^{2 sigma + 1}: 2.5e13 nodes at c = 71, 1.2e16 at c = 100: the record range is
  out of reach of the lattice minorant [WeilDefl](Riemann/Deep/WeilDefl/REPORT.md)

## Files (Riemann/Deep/WeilDefl/)
REPORT.md; DERIVATION.md; cert.py (the certificate; imports ../WeilZeta/cert_zeta.py unchanged); obstruct.py (sec. 7 estimates); ladder.py
(principal angles of the Phi-ladder; not part of the proof); e1_adapted.py (boundary-adapted ladders, float64 + Arb jets); dhscan.py (DH float64
Galerkin, explore17's gw.py); logs/: one JSON line per run (final_validation, final_lean*, rerun_lean, c25_*, c29_*, c31_*, run_c16pert_c20, dh_c16,
sharpness, weilzeta_c16_320bits; the c = 20 even line of run_c16pert_c20.log is the earlier plain ball-Cholesky version, "rigorous": false,
see DERIVATION.md sec. 4 (ii)), /usr/bin/time outputs (*.time; final.time also holds three runs refused by the checkpoint guard after a leaked
WD_CKPT variable, rerun in final2.time), obstruct.log, ladder.log, e1_adapted.log, dhscan.log, chek_manifests.log, rg_constants.log, review1.md
(derivation and code: SOUND), review2.md (table against logs: one unsupported entry, c = 16 odd lean, since rerun and corrected).
Scratch: .scratch/weildefl/STATE.md only.
