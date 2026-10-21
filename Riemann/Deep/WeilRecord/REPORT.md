NO DISCOVERY: Q(f) >= 1e-76 ||f||^2 on [-log 4, log 4] is now fully checked (all five review gaps closed, end-to-end re-run green: PROOF.md, run2/), but it is not a record. Elston (Zenodo 23174403, 2026-10-06; github zacharyelston/weil-windows) already claims a certificate on the larger half-width 1.495 (Q >= 3.5e-96 ||f||^2). Open gap: none in our proof; Elston's certificate itself is unchecked by us.

## Conventions (the brief's "Chuk" paper is by Xuefeng Zhu)
- Zhu, arXiv 2608.24827v2: supp f in [-L, L] (L = half-width), Q(f) = 2F(i/2)^2 + (1/2pi) int |F|^2 Psi_L, n < e^{2L}. Certified L = 0.8 ("support 1.6"), even f and (Sect. 6) all complex f; L = 1.19 only exploratory (earlier claim retracted).
- Our old files (Riemann/Deep/Weil, Yoshida) use the FULL width W: f in [-W/2, W/2]. Old cert 1.75 = Zhu L 0.875. Weil form identical in all three papers (pole 2Re(A conj B) = 2|int f cosh|^2 - 2|int f sinh|^2).
- Published rigorous claims, half-width L: Zhu 0.8; Mori Lean 4 ½log5 = 0.805 (formal); Liu 1 and 17/16 (MCOM submission, github luciferyu666); Mori VI ½log10, VII ½log12, VIII ½log15 = 1.3540 (Arb certificates, Zenodo 23238480/23265728/23281236, unrefereed). Not comparable: Mittermeier (Suzuki screw function, one function per t), Desogus (log 7, gaps). Scout + Zenodo/arXiv/GitHub API checks (2026-10-20): nothing beyond ½log15. CORRECTION (re-check, 2026-10-20): wrong. Elston (Zenodo 23174403) has half-width 1.495 > log 4; see PROOF.md section 8.
- Correction for DISCOVERIES.md: "Chuk" is Zhu; Mori's log 15 is an interval certificate, not heuristic; Mori's Lean log 5 proof (github moriryota/weil-positivity-log5-lean) is far beyond our Lean window 3/4 — re-judge weilPos_le_75.

## New result (bc.py, Arb/python-flint 0.9, ~25 min CPU)
- mu = e^{2L} = 16, prime powers 2..13. Proof chain (bc.py docstring), every lemma re-derived: Q >= <(M+V)f,f> with M = gamma - C_p, V = poles + band(a - gamma, |t| <= Omega) + rank-2 tail U; block elimination G^{-1} + J >= alpha + beta D; W(X) >= G; congruence H = W + WZW >= 0; all on M - tau.
- Parameters: ||C_p|| <= 3.642073270 (step-weight Schur test, 32000 cells), gamma = 9/2, Omega = 640 (C = 0.1236), N = 1400 Legendre modes, u = (343,114), 384 bits, tau = 1e-76.
- Margins: lam_min(T) ~ 1.9e-75 (even), 2.8e-71 (odd); quadrature 3e-99, projection alpha ~ 1e-150; Cholesky residual + radii <= 8e-81 < 2^-256. Logs: cert_mu16_tau1e-76.log, cert_mu16_tau1e-80.log.
- Validation: with Mori VI's exact parameters (mu=10, N=640) we get finite margin 8.45e-44 vs his "> 2^-144 = 4.5e-44"; mu = 3, 5: Q - (M+K) >= 0 against the independent x-space matrix (weil.py); mu = 5 certified (>= 1e-20).

## What is new in the method
1. Sharp prime-operator norm: positive step weight on a fine grid (Collatz-Wielandt fixed point) gives ||C_p|| to 5 digits (mu=15: 3.487 vs Mori's 4.125; mu=16: 3.642). Exact rational shifts (n = 2,4,8 when mu = 16) are needed, else boundary cells are miscounted.
2. No prolate defect: shifting M by tau gives an explicit constant 1e-76 (Mori: 2^-2904 type).
3. Gap rule (new NOGO): gamma - ||C_p|| must be about 0.85 at N/c ~ 1.6 (0.51 fails, slack 3e-60 > margin 2e-75, even with a two-level W at N2 = 1400); a small gap needs larger N/c (mu=10: gap 0.58 passes at N/c 3.5).

## Old pointwise certificate (brief item 2)
- Bug at window 1.875 was not assembly: arb's J_nu(x) start values (nu ~ 1230, x ~ 600-840) at fixed precision have radius up to 1e+97; the recurrence then overflows double range ("inf"). Fix (cert2.py): raise precision until the start value is accurate; also big Bernstein panels (2372 nodes, was 64800). Certified Q >= 9.99e-26 on W = 1.875 (Zhu L 0.9375). This method cannot pass Zhu's barrier R ~ 2pi e^{A_L} (W = 2 needs R ~ 2600); superseded by bc.py.

## Status and next
- External interval certificate, not Lean; unrefereed; trust = Arb + this code + derivations; reviewer (logic of steps 1-4, 2026-10-20): SOUND. Reviewer's verification gaps (Schur on all cells, e/n/h, quadrature propagation, Cholesky norm bound, end-to-end re-run) closed in PROOF.md and run2/ (2026-10-20). Window 3 (half-width 1.5, mu = e^3, gamma = 5.15, Omega = 1216, N = 3000): Gram matrices running.
- Reproduce (copy this folder's .py files next to Riemann/Deep/Weil/cert; venv with python-flint 0.9; numpy for schur_w.py): python3 schur_w.py 16 32000; bc.py schur 16 w_16_32000.txt; bc.py grams 16 1400 384; bc.py band 16 1400 640 9/2 384; bc.py cert 16 grams_16_1400_384.bin band_16_1400_640_9o2_384.bin 36421/10000 1/2 1e-76 343 114 384 - 1400.
