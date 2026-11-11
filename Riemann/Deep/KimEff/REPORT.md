RECORD: J^{d,n} is hyperbolic for every n >= 998,616,377,771 and every d <= 6.45 n^2/log^2 n, an explicit Kim–Lee wedge n >~ 0.394 sqrt(d) log n (for example c = 0.55, C = 0.8716: every d and every n >= max(n_*, C d^c)); it contains the two-thirds wedge n^2 log(n+2) >= 24 d^3 for every d > 1.0471e8 and uses no RH verification.

# Explicit Kim–Lee-shape hyperbolicity region for Riemann's Jensen polynomials (node N6)

## Claim

J^{d,n}(X) = sum_j binom(d,j) gamma(n+j) X^j, 8 xi(1/2+z) = sum gamma(m) z^{2m}/m! (the N6 family; Theta(w) = 8 xi(1/2 + sqrt w),
Kim–Lee's Xi_0(-w)). Let u(n) solve u(pi e^{2u} - 9/4) = n, n_* = ceil(Z(12)) = 998,616,377,771, and
R_n = (9/10)(n - 1/2)^2/(u(n) + 1)^2. With a = 0.37338196661320752 (the Lean-checked half-plane |Re rho - 1/2| <= a):

1. For every n >= n_*, Theta^{(n)} has no zero w in the parabola region (Im w)^2 <= 4a^2(a^2 - Re w), Re w >= -R_n.
2. Hence J^{d,n} is hyperbolic for every d <= D_n = R_n^2/(4a^2(R_n + a^2)), in particular for every d <= 6.45 n^2/log^2 n.
3. Kim–Lee shape, explicit: J^{d,n} is hyperbolic for every d >= 1 and every n >= max(n_*, C d^c) with
   (c, C) = (0.51, 7.6890), (0.55, 0.87160), (0.6, 0.069829), (1, 1.1853e-10) (C rounded up).
4. With only the classical strip (a = 1/2): d <= 3.59 n^2/log^2 n, and (c, C) = (0.55, 1.2031), (0.6, 0.099247).

Comparison (combined with the all-n range D_all = 1.6141e25): D_n > D_all for all n >= 4.6974173e13; at d = D_all + 1
the two-thirds wedge needs n >= 3.41725e37, ours n >= 4.69742e13; for every d > D_all our region contains the two-thirds region. Without the all-n range it contains the two-thirds
region at every d > 1.0471e8.
At n = 1e30: d <= 1.46e57 (two-thirds wedge: 1.42e20; all-n range: 1.61e25).

## Before

- Explicit: J^{d,n} hyperbolic whenever n^2 log(n+2) >= 24 d^3 (ours, savvak, two-thirds wedge, exponent 3/2), plus the all-n
  range d <= 1.6141e25; self-contained 2.76e25 d^3.
- Best shape: Kim–Lee (arXiv 2105.05386, Thm 1): N(Xi_0; d) = O(d^c) for every c > 1/2, ineffective through Kim's 1996
  derivative-zero theorems (PAMS 124, 819-830, Thm 2; JKMS 33, 657-667, Thm 1). No effective version was known (R413).

## Method

Kim–Lee's chain is: zero location of Xi_0^{(n)} (Kim 1996, ineffective) -> Obreschkoff's sector theorem (their Thm 3,
explicit) -> J^{d,n} hyperbolic. We keep the second half (it is exactly the Chasse/Obreschkoff step of the N6 all-n
record) and replace Kim's theorem, for Xi only, by an explicit zero-free region:

1. Gauss–Lucas (genus 0): all zeros of Theta^{(n)} lie in the convex parabola P = {(Im w)^2 <= 4a^2(a^2 - Re w)}, which
   holds the zeros of Theta. So if Theta^{(n)} has no zero in P with Re w >= -R, every zero has angle
   |Im w|/|w| <= 2a sqrt(R + a^2)/R, and Obreschkoff gives hyperbolicity for d <= R^2/(4a^2(R + a^2)).
2. Mixture identity (from Riemann's integral xi(1/2+z) = 4 int phi(t) cosh(zt) dt):
   Theta^{(n)}(w) = gamma(n) E_n[ 0F1(; n + 1/2; t^2 w/4) ], E_n over the probability density p_n ∝ t^{2n} phi(t).
   0F1(; n+1/2; x) is a Bessel function whose zeros are -j_{n-1/2,k}^2/4 <= -(n-1/2)^2/4, and p_n is concentrated at
   t = u(n) (relative width about n^{-1/2}). So Theta^{(n)} is an average of dilates of one Laguerre–Pólya function,
   each with a zero-free disc of radius about n^2/u^2, over an extremely narrow window of dilations.
3. Explicit Laplace bounds on the tilted density e^{psi}, psi = log(t^{2n} phi_1) + log 0F1(..; -t^2 s/4), which is concave:
   peak bracket [u - 17/25, u + 1/n] (uniform in w), a window of half-width 2 eta (eta <= 6.86e-6) around the peak, a
   phase bound (the argument of 0F1 turns by <= 1.26e-4 rad on the window, from Im(w G'/G) = v sum c_k/|x + c_k|^2),
   tangent-line bounds for the side tails, and a crude bound e^{|x|/b} beyond t = U = 1.0435(u+1) where p_n is doubly
   exponentially small. Bessel input: product formula, Rayleigh sums sum 1/c_k = 1/b, sum 1/c_k^2 = 1/(b^2(b+1)), and
   j_{nu,1} > nu (DLMF 10.21.15, 10.16.9, 10.21.3). theta = 9/10 is the fraction of the first Bessel zero used.
4. Every inequality reduced to monotone majorants in u = u(n) (Sec. 9 of DERIVATION.md) and certified at u = 12 in Arb.

The method does not make Kim's general theorems effective; it uses the special structure of Xi (Pólya's positive kernel).
Kim's papers were not needed (AMS served a bot check to curl; noted in NOVELTY.md).

## Certificate

    cd Riemann/Deep/KimEff && uv run --no-project --with python-flint python3 certify.py          # a = 0.37338196661320752
    cd Riemann/Deep/KimEff && uv run --no-project --with python-flint python3 certify.py 0.5      # a = 1/2

Under 2 s each; prints every constant with directed decimal rounding (lower bounds down, upper bounds and sufficient
thresholds up, each verified against the Arb ball) and "ALL CHECKS PASS" (logs/certify.log, logs/certify_half.log).
Key values at u = 12: eta <= 6.8561470e-6, phase <= 1.2507254e-4, E+ <= 7.8678191e-9, c6 >= 0.78582827,
crude-tail log ratio <= -7.84740e11, FINAL >= 1.8747279 (> 0 needed).

## Checks (sanity, not part of the proof)

- norm.py (logs/norm.log): gamma(m) = 32 m! Mt(m)/(2m)! for m = 0..6 to 15 digits (fixes the normalization of the kernel).
- sanity1.py (logs/sanity1.log): mixture identity at n = 3, w = -50 + 10i: Taylor sum and mixture agree to 12 digits.
- sanity2.py (logs/sanity2.log): Theta^{(n)}(-s)/gamma(n) is positive and decreasing on a 60-point grid of (0, R_n] for
  n = 50 and 200, and its first sign change on a 3% grid lies well beyond R_n (values in the log; the proof needs u >= 12,
  these small n only test the scale). n = 2000 (logs/sanity2_big.log): positive and decreasing on the grid,
  value 2.9e-110 at R_n; an n = 10000 run exceeded its 600 s limit and was dropped.
- sanity3.py (logs/sanity3.log): j_{nu,1} > nu, the Rayleigh sums, and the Bessel bounds (5.4)-(5.8) on 120 random points
  (b = 5.5, 20.5, 60.5): 0 failures.
- The certificate's grid u = 12..300 shows every majorant moving in the monotone direction claimed in DERIVATION Sec. 9.
- rg of every DERIVATION constant in certify.py: Sec. 11 table lists the variable for each.

## Novelty

arXiv listings 2026-10 and 2026-11 (math.NT, math.CV, math.CA, lst.py 'Jensen|hyperbolic|derivatives|Laguerre|Polya'):
nothing on explicit large-n Jensen regions. arXiv API (nonreal zeros + derivatives; Jensen + effective/explicit; Jensen +
Bessel; repeated differentiation of xi): no effective Kim–Lee. OpenAlex: Kim–Lee 2105.05386 has 0 citers; Kim PAMS 1996
has 13 citers, none effective for Xi_0 (checked titles); Kim JKMS 1996 0 citers. Zenodo "Jensen polynomials" (60 records,
newest 5 Nov 2026) and HAL: nothing of this shape (Planat's hal-05675741 fits a constant numerically). GXL evidence API:
quota exhausted (429). Our files: SavvaK (two-thirds wedge), HollandK, Jensen9 - weaker. Rechecked 11 Nov 2026 (listings through 2611.13932,
Zenodo through 10 Nov): nothing new. Details: NOVELTY.md.

## Caveat

- Cited, not re-proved: Obreschkoff's theorem (as Kim–Lee Thm 3), Riemann's integral for xi, the Bessel product formula,
  0F1–Bessel relation and j_{nu,1} > nu (DLMF), Hadamard, Hurwitz; the half-plane a = 0.37338... is a repo input
  (a = 1/2 version is classical-only).
- Proof is paper + Arb, not Lean. Review so far: four verdict-only reviewer subagents on parts A-D of DERIVATION.md
  (GAP/ERROR on presentation: an "Equivalently" that should be "Consequently", a rounded-down printed constant
  (C for c = 1, M1 = 4.0000000), missing hypotheses nu >= 0 and c > 1/2, the X < 0 case of two series bounds, the Hc
  concavity, a containment monotonicity claim that was false for n < 8); all repaired (directed rounding in the script,
  added derivations); no error in the mathematics of the chain was reported. Second round, two verdict-only reviewers
  on the full revised DERIVATION.md (Secs. 0-10 pasted inline, Arb enclosures assumed): CORRECT, CORRECT.
  The orchestrator's own reviewer and a rerun of the certificate are still to come.
- Only hyperbolicity is claimed (not distinct zeros). theta < 1 is a limit of the method (the frozen Bessel zero); the
  true first zero of Theta^{(n)} is further out, so the constant 6.45 is not sharp, but the exponent 1/2 (up to log) is
  Kim–Lee's.

RH relevance: N6 (RH iff every J^{d,n} is hyperbolic; the explicit hyperbolic region grows from n >~ d^{3/2} to
n >~ 0.394 sqrt(d) log n, Kim–Lee's shape made effective; uniformity in d, which is RH itself, stays out of reach:
the method needs n >~ sqrt(d) log d)

## Files

certify.py (certificate), logs/certify.log, logs/certify_half.log, DERIVATION.md (every step and constant), NOVELTY.md,
norm.py, sanity1.py, sanity2.py, sanity3.py and their logs.
