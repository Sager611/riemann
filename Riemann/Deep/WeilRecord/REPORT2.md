DISCOVERY: certified half-width 3/2 (window length 3, x = e^3): Q(f) >= 1e-100 ||f||^2 for every f supported in [-3/2, 3/2], both parities (Arb ball arithmetic, chain of PROOF.md). Strictly beyond Elston's 1.495 (Zenodo 23174403), the largest posted certificate found; the gain is small (+0.005 in half-width).

## Result (2026-10-20, .scratch/weilrec3, logs copied to run3/)
- L = 3/2, mu = e^3, prime powers 2..19. H1: ||C_p|| <= 4.295809746 (weighted Schur test, all 32000 cells) <= s = 4.2959. gamma = 103/20, Omega = 1216 (C = a(Omega) - gamma = 0.11544), u = (660, 219), N = 3000 Legendre modes (N/c = 1.645), 384 bits, tau = 1e-100.
- Truncation (rho = 3): e = 2.6e-371, n = 5e-746. Quadrature: epsJ <= 4.1e-100 (even), 1.4e-100 (odd).
- Even: H_mid - 2^-326 I = LL^T + Res, ||Res||_F + ||rad H||_F <= 1.7e-108 < 7.3e-99, so lambda_min(W^-1 + Z) >= 5.95e-100. Odd: lambda_min >= 9.97e-96.
- Cost: Gram matrices 3 parallel bc.py processes (~7 CPU h, exact Gauss), band 25 min, certificate 21 min / 2.3 GB peak (cert3.py).
- Elston's convention checked: f in [-a, a], window length 2a = log x, n < x, same functional; a = 1.495 certified >= 3.50e-96 (even), 8.26e-92 (odd), Rayleigh-Ritz upper bound 6.27e-96; 28 workers x ~6 h per sector; 2.99 run "regression-tested", not audited.

## New tools (Riemann/Deep/WeilRecord)
- cert3.py: same proof chain as bc.py, but the floating Cholesky runs at prec + 128 bits (bc.py's leaves ||Res||_F ~ 1e-80 at 384 bits, above the ~1e-99 margin here: the old e3 chain would have failed), sigma from 1/||L^-1||_F^2, one parity in memory, H1 mandatory, raw grams_c input, band computed in-process (no 0.5 GB file).
- grams_c.c: C/FLINT port of the Gram step; bit-identical to bc.py at mu = 5, 16, e^3 (small N). "bound" mode: per-piece n-point Gauss with a rigorous Bernstein-Walsh bound (|P_j(z)| <= rho(z)^j, convex a(z) maxed on a circumscribed polygon, ATAP Thm 19.3) added to every entry: 4.7x fewer nodes at e^3.
- bc.py Setup: new tags L<p>o<q> (rational half-width), backward compatible.

## Novelty checks (2026-10-20)
- Zenodo API (phrase searches "Weil positivity", "Weil quadratic form", window/certificate): newest relevant records (23515504 crosses log 5/2; 23447303 is an upper bound) are all below 1.495. GitHub API (weil-windows last commit 2026-10-06; moriryota log 15; luciferyu666; teal-sea/zeta-lab) and arXiv API: nothing beyond 1.495.

## Status and caveats
- External interval certificate, unrefereed; trust = Arb + this code + PROOF.md lemmas (unchanged). Not a Lean proof.
- Running (next2.sh): cert3 replay of the log 4 certificate; bound-mode grams at e^3 checked against the exact ones; then half-width 1.55 (||C_p|| <= 4.591281212, gamma = 109/20, Omega = 1640, N = 4200, 448 bits, tau = 1e-112). Next after that: log 5 (||C_p|| <= 4.960357039; N ~ 6300).
