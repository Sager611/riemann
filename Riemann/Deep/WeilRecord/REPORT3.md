DISCOVERY: certified half-width 1.55 (window length 3.1, x = e^3.1 = 22.2): Q(f) >= 1e-112 ||f||^2 for every f supported in [-1.55, 1.55], both parities (Arb ball arithmetic, PROOF.md chain). Beyond our 3/2 (REPORT2) and Elston's 1.495; no other posted certificate above 1.495 found (2026-10-21).

## Result (2026-10-21, .scratch/weilrec4, logs in run4/)
- L = 31/20, prime powers 2..19. H1 ||C_p|| <= 4.591281212 <= s = 4.5913 (weighted Schur, all 32000 cells). gamma = 109/20 (gap 0.8587), Omega = 1640 (C = 0.11457), u = (913, 303), N = 3400 Legendre modes (N/c = 1.34, was 1.65), 448 bits (Cholesky 576), tau = 1e-112.
- Truncation (rho = 2.2): e = 1.1e-196, n = 6.9e-397. Quadrature: epsJ <= 8.2e-120 (even), 2.9e-120 (odd). ||rad H||_F <= 1.2e-127.
- Even: H_mid - 2^-365 I = LL^T + Res, ||Res||_F + ||rad H||_F <= 1.4e-122, so lambda_min(W^-1 + Z) >= 1.2e-111. Odd: >= 4.0e-107.
- Slack of the chain (floating diagnostics, not proof): lambda_min(T) ~ 1.52e-108 vs Galerkin lambda_min(A + J) ~ 1.67e-108 (even, ratio 0.91); 3.95e-104 vs 4.39e-104 (odd).
- Cost: grams 2 x ~2 CPU h (grams_c, bound mode), certificate 55 min CPU / 4 GB peak; machine heavily shared (wall 2.6 h).

## Why the first 1.55 run failed, and the fix (cert4.py)
- cert3 flushed Bessel values below 2^-400 to zero and charged 1.9e-112 per band entry; with the U-term weight u0 = 913 that gave epsJ = 7.2e-106, above the even margin (~1e-108): the even FAIL was bookkeeping, not positivity.
- cert4.py: flush threshold 2^-(prec+160); column error sqrt(d) errB (was d errB); beta D replaced by beta (b - tau)^2 I (D never formed); W(X) = I/m + X^T Q X + SX + (SX)^T (3 products); band keeps only K and Bb's first column; one parity in memory; grams parts summed per parity with a check that the pieces cover 0..npieces-1 exactly once; part files integrity-checked (a full disk once truncated a part silently; grams_c.c now checks ferror and removes bad output).
- Checked: cert4 reproduces the mu = 5 smoke margins of bc.py exactly (2.777e-21, 2.412e-17); split and unsplit grams give identical results.

## Negative result
- Gap 0.65 (gamma = 131/25, Omega = 1340, N = 2825): even certified (slack ratio 0.52) but odd H_mid not PD. Gap ~0.86 stays necessary (NOGO).

## log 5 (x = 25, half-width 1.6094)
- STOPPED (orchestrator, machine out of RAM): planned gamma = 291/50, Omega = 2380, u = (1405, 468), N = 5210 (N/c 1.36), 480 bits, Cholesky 640, tau = 1e-127; predicted even margin ~1e-122. Estimated certificate peak ~10 GB (N = 3400 peaked at 4 GB, memory ~ d^2): needs an out-of-core/blockwise cert4 under 2 GB first.

## Novelty checks (2026-10-21)
- Zenodo API ("Weil positivity", "Weil quadratic form"): no record after 2026-10-19. GitHub (weil-windows last commit 2026-10-05; moriryota max log 15; new repos: jose-paul-b <= 31/32, davidangularme log 2, darreal44 mu 11, aidilamrym-ops finite matrix only). arXiv 2609.04908 (Shi): finite Prime-Weil matrices, no window certificate.

## Status
- External interval certificate, unrefereed; trust = Arb + this code + PROOF.md lemmas. Not a Lean proof.
