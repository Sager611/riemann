# Explicit constant in Holland's hyperbolicity wedge — derivation

Source: G. Holland, arXiv 2608.08682 (TeX `jensen.tex`, line numbers below). Certificate: `certify.py` (Arb ball
arithmetic, python-flint). Every numbered inequality below is evaluated there; the script aborts on any failure.

## 0. Statement and notation

`J^{d,n}(X) = sum_{j<=d} binom(d,j) gamma(n+j) X^j`, `xi(1/2+z) = sum gamma(n) z^{2n}/n!` (l. 57-77).
`d0 := 16,140,857,488,599,731,839,205,676` (our all-n degree record, COROLLARIES/TREE N6).

**Theorem.** Put `L* = 90`, `N* := Z(L*/2) = 45 (pi e^{90} - 9/4) ≈ 1.725e41`.
(i) For every `D_all >= 9e24` up to which J^{d,n} is hyperbolic for all n (e.g. `d0`, or T^2 with T = 3.000175e12), all integers `d > D_all` and `n >= 0` with `n^3 log^2(n+2) >= 7000 d^5`, `J^{d,n}` has `d` distinct negative zeros.
(ii) Hence `J^{d,n}` is hyperbolic for every `(d,n)` with `d <= D_all` or `n^3 log^2(n+2) >= 7000 d^5`.
(iii) Self-contained form of Holland's Theorem 1.1 (l. 81-89), not using any all-n range: with L* = 80,
`K = N*^3 log^2(N*+2) <= 2.430e114` works for all `d >= 1` (for d <= 4, P_F = p_F, l. 1825-1828; `certify_self.log`).

Notation. `Z(u) := u (pi e^{2u} - 9/4)` (increasing on `u > 0`), `u_0(x) := Z^{-1}(x)`, `L := 2 u_0(n)`,
`eps := 1/L`, `nu := 1/n`. Then `n >= N*  <=>  L >= L*  <=>  eps <= 1/90`. Also `L <= log n` and
`log(n+2) <= L + log(pi L/2) + 2/n`.

Kernel normalization. With `phi(u) = sum_k (2 pi^2 k^4 e^{9u/2} - 3 pi k^2 e^{5u/2}) e^{-pi k^2 e^{2u}}` put
`Mt(z) = int_0^oo phi(u) u^{2z} du`. Then `gamma(z) = C * 4^{c z} * Gamma(z+1) Mt(z)/Gamma(2z+1)` for constants
(checked numerically, `sanity.py`: the moment ratios are exactly 4, 16, 64 for n = 0,1,2). Only
`h^{(k)}`, `k >= 2`, of `h := log Mt` enter the proof, so the constants never matter (l. 503-516 make the same remark).

## 1. Architecture (Holland, Prop. 2.2, l. 228-315) and the one change we make

Prop. 2.2 with Cauchy radius `lambda r` instead of `2r` (identical proof, l. 280-300): if `c` is holomorphic near
`Omega := {dist(z,[0,d]) <= lambda r}`, `c(0..4) = 1`, `c(d) > 0`, `sup_Omega |c-1| <= e < lambda^4 (lambda-1)`, and
`|y^k p^{(k)}(y)/p(y)| <= r^k` at every critical point of `p`, then `|P/p - 1| <= e sum_{k>=5} lambda^{-k}
= e / (lambda^4 (lambda-1)) < 1` there, and `P` has `d` simple positive zeros. We use `lambda = 1.4`
(`lambda^4(lambda-1) = 1.53664` exactly); Holland's `lambda = 2`, `e < 16` also passes (K about 4x larger).

Everything else follows Holland's chain: Laguerre defect (l. 622) -> Jacobi model (l. 705-963) -> four-coefficient
finite-free model (l. 965-1278) -> root localization (l. 1371) -> fifth-order residual (l. 1523) -> derivative ratios
(l. 1636) -> assembly (l. 1812). Each "sufficiently large n" becomes `n >= N*`; each `O(.)` becomes an Arb enclosure.

## 2. Part A — explicit complex saddle (replaces Prop. 4.1, l. 471-553, i.e. the complex GORTTW (3.2))

Finding: Holland's analytic input (l. 508-520) is "the complex form of [GORTTW, (3.2)], stated in the sentence following
that equation". In GORTTW (arXiv 1910.01227, TeX l. 514-518) that is one sentence ("the analytic continuation of L_M and
Stirling's formula imply that even when M is complex, we may keep the existing error term") with no proof and no
constant. We prove an explicit version.

Fix `x0 >= N*`, `a0 := u_0(x0) >= 45`, `D(x0) := {|z - x0| <= x0/10}`.

A1 (saddle). For `z in D(x0)`, `T_z(u) := (1/2) log((z/u + 9/4)/pi)` maps `{|u-a0| <= 0.06}` into itself with
`|T_z'| <= 1/(2(a0-0.06)(1-w1)) < 0.012`, `w1 := 9(1+0.06/a0)/(3.6(pi e^{2a0}-9/4))`. Its fixed point `u_0(z)` solves
`Z(u) = z`, is holomorphic, real for real z, and `|u_0(z) - a0| <= beta := (1/2)(-log 0.9 - log(1-0.06/a0) + 4 w1)
<= 0.0533474` (a0 >= 45, L* = 90; `<= 0.0534309` for a0 >= 40, L* = 80; outward-rounded, `certify.log`). (|log(z/x0)| <= -log 0.9; |log(1+w)| <= 2|w|.)

A2 (integral). `f(u) := 2z log u + 9u/2 - pi e^{2u}`, `Psi(z) := f(u_0(z))`, `a(u) := 1 - (3/2pi) e^{-2u}`.
`Mt(z) = 2pi^2 int_0^oo e^{f} a du + sum_{k>=2} int phi_k u^{2z} du`. Shift `[0,oo)` to
`[0, al/2] + [al/2, al/2 + i beta'] + {u_0 + s: s >= -al/2}` (`u_0 = al + i beta'`; the integrand is analytic in
`Re u > 0` and decays doubly exponentially in `|Im u| < pi/4`). On the ray `F(s) := Re(f(u_0+s) - f(u_0))` is
concave: `Re f'' <= -4 pi cos(2 beta) e^{2(al+s)}` because `Re(z/u^2) >= 0` (`asin 0.1 + 2 atan(2beta/al) < pi/2`).
With `s1 = 0.3`, `C3 := sup|f'''|/6 <= c3 e^{2al}`, `c3 := [4(al+beta)(pi + (9/4)e^{-2al})/(al-s1)^3 + 8 pi e^{2 s1}]/6`
at `al = a* - beta` (decreasing in al): `c3 <= 7.63354` (L* = 90), `<= 7.63382` (L* = 80); note `c3 > (4pi/3) e^{0.6} = 7.63247..`,
the real-axis limit of `|f'''(a+0.3)|/(6e^{2a})`. And `C3 s1 <= kappa0/4`
(`kappa0 := Re(-f''(u_0)) >= 4 pi cos(2beta) e^{2al}`):
  core error `<= 16 C3/kappa0^2`, Gaussian tails `<= 2e^{-kappa0 s1^2/2}/(kappa0 s1)`, ray tails (tangent line of the
  concave F) `<= 2 e^{-kappa_min s1^2/2}/(kappa_min s1)`, amplitude `|a-1| <= (3/2pi) e^{-al}`;
  the pieces near 0 are `<= e^{5.2} x0 e^{b x0}` relative, `b = 0.9 log(1/4 + beta^2/al^2) + 0.6 beta/al + 1.1/al <= -1.2224` (L* = 90; `-1.2193` for L* = 80);
  the theta terms k >= 2 are `<= C_theta 2^{-4.5}/(2pi^2) * exp(-g (pi e^{2a0} - 9/4)) * (e^{2beta+1} sqrt(phi2) + 1)`
  relative, `C_theta = sum_{k>=2} (2pi^2 k^4 + 3pi k^2) e^{-pi(k^2-4)}`, `g = 0.9 log2/(1+beta/a0) - q2max/(360(1-beta/a0)) >= 0.61746` (`0.61736` for L* = 80; `q2max <= 2.0225`)
  (deficit `Psi(x) - Re Psi(z) <= (y^2/2) sup|Psi''|` versus the gain `x log2/u_0(x)` from `(w - log 2)^{2x} <= w^{2x} e^{-2x log2/w}`).
Result: `Mt(z) = 2 pi^2 e^{Psi(z)} sqrt(2pi/(-f''(u_0(z)))) (1 + delta(z))`, `|delta| <= deltaA`, where
`deltaA = rel_core + rel_gt + rel_rt + rel_am + rel_pc + rel_k` with `rel_core = c3 sqrt(2 phi2) e^{-al}/(pi^2 cos^2(2beta))`,
`rel_am = (3/2pi) e^{-al} (sqrt(phi2 e^{2 s1}/cos(2beta)) + rel_rt)` and four super-exponentially small terms, all at
`al = a* - beta`: `deltaA <= 5.33e-20` (L* = 90: rel_core <= 3.36e-20, rel_am <= 1.97e-20, the rest < 1e-300) and
`deltaA <= 7.91e-18` (L* = 80). Inputs: `phi2 <= 1.01113`, `cos(2beta) >= 0.994313` (L* = 90).
Checked against quadrature (`sanity.py`): relative error 2e-3 at x0 = 2000, also for complex z (consistent with O(e^{-2a})).

A3 (logarithm). `Mt != 0` on `D(x0)`; `h := log Mt = const + Psi(z) + G(z)` with
`G := -(1/2) log(-f''(u_0(z))) + log(1+delta)`, `-f''(u_0) = (4z/u_0)(1 + 1/(2u_0) + 9u_0/(4z))`, and
`sup_{D(x0)} |G - G(x0)| <= etaG := (1/2)(-log0.9 - log(1-beta/a0) + 2(beta/(2 al a0) + 2 w1)) + 4 deltaA <= 0.053287`
(L* = 90; `<= 0.053366` for L* = 80).
Cauchy: `|G^{(k)}(w)| <= k! etaG / (x0/10 - |w-x0|)^k`.

## 3. Part B — the main term

`Psi'(z) = 2 log u_0(z)` (envelope), `Psi'' = 2/(u Z'(u))`, `d/dz = (1/Z'(u)) d/du`. With `t = 1/u_0`,
`s = 9/(4 pi e^{2u_0})`: `Psi^{(k)}(z) = Q_k(t,s) / (z^{k-1} 2u_0(z))`, `Q_k` rational (sympy, exact integer
coefficients): `Q_2 = 4(1-s)/(2+t(1-s)) = 2 - t + ..`, `Q_3 = -2 + t^2 + ..`, `Q_4 = 4 + t + ..`, `Q_5 = -12 - 5t + ..`
(Holland's eq. hsecond..hfifth, l. 563-586, with explicit errors). `|Q_5| <= 14.968` (L* = 90; `15.526` for L* = 80) on the complex polydisc
`|t| <= 1/(a* - 0.07)`, `|s| <= tiny`. So `h^{(k)} = Psi^{(k)} + G^{(k)}`, `k >= 2`, with explicit enclosures.

## 4. Defects and Jacobi parameters (l. 622-963), normalized and enclosed in Arb per eps-cell

`X0 := nL delta0`, `delta_j := Delta^2 h(n+j)` (averages of `h''` on `[n+j, n+j+2]`), `Y := n^2 L Delta^3 h(n)`,
`W := n^3 L h''''(n+z)` (z in [0,4]), each = `Q_k * [(n/t')^{k-1} u_0(n)/u_0(t')] +- g_k`, `g_k = k! etaG (10/(1-50nu))^k L/n`.
`lam0 := n log q0 = n ell(b) - eps X0 > 0` (q0 > 1, l. 624), `n ell(b) in [1-nu, 1]`.
`n^2 L D_n = Y + 2 n/(b+2) X0 + (1/2)F''(xi) eps X0^2`, `(1/2)F'' in [-1, -1+3nu]` (F = log(2-e^{-x}); exact
Taylor form of l. 689-701) `> 0` (l. 653).
`nL(T-q0) = n(q0-1) n^2L(q1^L - q1) / (n(q1-1) n(q1^L-1))` (identity l. 763), `aJ := A_J/(nL) = 2q0/(nL(T-q0))`,
`YB := nL(1/b - 1/B_J) = (1+1/b) X0 (1-e^{-delta0})/delta0 - q0/aJ > 0`, `betaJ := B_J/n = 1/(n/b - eps YB)`,
`bJ := (B_J-b)L/n = (b/n) betaJ YB`.
`n^3 L E_J''''(z) = W - M [1, 1+4nu] - [0, 2 eps^2/aJ^3 (1+3nu)]`, `M = 2 bJ p q (p^2+pq+q^2)`, `p = n/(b+z)`,
`q = 1/(1/p + bJ eps)` (exact integral of `-6/x^4`; `psi'''' in -(6/x^4)[1, 1+4/x]`). `tQ := n^3 L Q_n = -avg >= 1.9116`
(L* = 90; `>= 1.9008` for L* = 80) (l. 921-926); `n^2 L D_n >= 1.9129` (L* = 90; `1.9023` for L* = 80).
Cross-check at n ~ 1e45 from the main term at 600 digits (`sanity2.py b`): `X0 = 1.9802, n^2LD_n = 1.9216, aJ = 1.0202,
tQ = 1.9214` at L = 100, inside the certified enclosures.

## 5. Four-coefficient matching (Lemma 6.1, l. 1012-1220): Brouwer instead of the 1/30 contraction

Unknowns `A = a A_J`, `B = B_J + b n/L`, `C = D + c n/L`, `D = B_J/2` (l. 999-1003, 1052-1056). With
`g(x) = log(1+1/x)`, `g' = -1/(t(t+1))`, `(Delta g)' = 2/(t(t+1)(t+2))`, `(Delta^2 g)' = -6/(t(t+1)(t+2)(t+3))`
(exact), the rescaled system `Phi_n = Phi_inf + Delta` with `Phi_inf = (1 - 1/a - b + 4c, 2b - 16c, -6b + 96c - 2)`
(l. 1110-1117) and
  Delta_0 = (1/aJ - 1)(1-1/a) + (1/aJ)((lam-1) - (lam-1)/a) + b(1 - avg 1/P1) + c(avg' 1/P1 - 4),
  Delta_1 = +-2|a-1|eps/(min(a,1)^3 aJ^2) + b(avg 2/P2 - 2) + c(16 - avg' 2/P2),
  Delta_2 = +-6|a-1|eps^2/(min(a,1)^4 aJ^3) + b(6 - avg 6/P3) + c(avg' 6/P3 - 96) + (2 - tQ),
`Pm(tau) = tau(tau+nu)...(tau+m nu)`, avg over `tau in [betaJ, betaJ + b eps]`, avg' over `[betaJ/2, betaJ/2 + c eps]`.
`Phi_n(x) = 0` iff `x = S(x)`, `S = (1/(5/6 + 3/4 D1 + D2/12 + D0), 1/3 - D1 - D2/6, 1/24 - D1/16 - D2/48)`.
The certificate checks `S(box) ⊂ int(box)`, box `[1.10,1.30]x[0.25,0.42]x[0.015,0.07]` (Brouwer => solution), then
shrinks the enclosure by `encl := S(encl) ∩ encl` (6 times). Ranges (l. 1035-1039) checked: `A >= 8B`, `B, D >= 32d`,
`C - D >= 4d` (as `xc > 4 log(N*) d0/N*`), `C - D <= D/4`.

## 6. Root localization (Lemma 7.3, l. 1371-1513), sharpened Gershgorin

Jacobi-matrix entries (l. 1403-1410) and the diagonal-difference identity (l. 1421) verified numerically
(`check_jacobi.py`). Bounds (k <= d-1, `H = U-d-1`):
p1 = q_{A,B}: `|diag - B| <= (2(1+2d/H)^2 + 2B/A) d`, off-diagonal `<= (1+2d/H) sqrt(1+d/B) sqrt(Bd)` (using
`(k+beta)(k+H) <= (2k+H)^2-1`), so roots in `B +- CJ1 sqrt(Bd)`, `CJ1 = 2(1+2d/H)sqrt(1+d/B) + cdiag1 sqrt(d/B) <= 2.00000002` (L* = 90)
(Holland: `C_J = 8`, l. 1455).
p2 = q_{C,D}(D.): `|diag - D| <= (3 + 2d/D) d`, off-diagonal `<= (1+2d/D)(1+d/D) sqrt(Gd)`, `G = C-D`, so
`eta = (3+2d/D) d/D + 2(1+2d/D)(1+d/D) sqrt(Gd)/D < 1/2`.
Lemma 7.1 (MSS Thm 1.13, l. 1302-1328): every root and critical point y of p_F has `|y - B| <= Cloc sqrt(Bd)`,
`Cloc = CJ1 (1+eta) + eta sqrt(B/d) <= 2.0874` (L* = 90; `<= 2.0930` for L* = 80; `CJ1 <= 2.00000002`, resp. `2.00000306`)
(Holland: unnamed absolute constant, l. 1379).

## 7. Derivative ratios (Lemma 9.1, l. 1636-1810)

Identities (l. 1676-1680 factorization, l. 1717-1722 V_m expansion, l. 1745-1749 Jacobi recurrence) verified symbolically
(`check_identities.py`). With `r = K6 sqrt(Bd)`, `K6 = 2.6` (Holland: "sufficiently large", l. 1535):
`a_m <= (Cloc + sqrt(d/B)(1+y/A)) / (K6 (1-y/A))`, `b_m <= (1+d/A)(1 + Cloc sqrt(d/B)) / ((1-y/A) K6^2)`,
`theta <= (B/D)(1+Cloc sqrt(d/B)) eps_p C_V / ((1-y/A)(1-r/D))`, `eps_p = (C-D)/C`,
`C_V = 1 + r/A + 3d/A + 2 sqrt(d/B)/K6 + 10 (d/A) sqrt(d/B)/K6 + d/(4K6^2 B) + ..` (explicit form of C_*, l. 1723-1726).
The contradiction step (l. 1797-1807) needs only `a + b + theta < 1` (Holland asks each <= 1/4): certified `<= 0.97084`
(L* = 90: a <= 0.81815, b <= 0.15075, theta <= 0.0019431; `<= 0.97568` for L* = 80); `C_V <= 1.00000001` (L* = 90).

## 8. Fifth-order residual (Lemma 8.1, l. 1523-1627), signed

`n^4 L E_F^{(5)}(z) = Q5(t) (n/(n+z))^4 (u_0(n)/u_0(n+z)) + 6(bJ+xb) P1 P2 (P1^3+P1^2P2+P1P2^2+P2^3)
 - 6 xc P3 P4 (P3^3+P3^2P4+P3P4^2+P4^3) + err`, `P1 = n/(b+z)`, `P2 = 1/(1/P1 + eps(bJ+xb))`, `P3 = n/(D+z)`,
`P4 = 1/(1/P3 + xc eps)` (exact integrals of `24 w^{-5}`; `|psi^{(5)}(w) - 24 w^{-5}| <= 240 (Re w)^{-6}` by
Euler-Maclaurin), err = G-term + psi-remainders + `psi''''(A+z)`. Evaluated in complex balls over
`|z| <= R0 = d + lambda r`: `|n^4 L E_F^{(5)}| <= ce5 <= 12.769` (L* = 90; `12.981` for L* = 80) (the three main terms nearly cancel: -12 + 32 - 32;
the triangle inequality, as in l. 1586-1603, gives 84). Hermite-Genocchi (l. 1607-1614):
`sup_Omega |E_F| <= (ce5/120)(d + lambda r)^5/(n^4 L) <= Cfin d^{5/2}/(n^{3/2} L)`,
`Cfin = (ce5/120)(lambda K6)^5 (B/n)^{5/2} (1 + sqrt(d/B)/(lambda K6))^5 <= 70.600` (L* = 90; `72.118` for L* = 80).

## 8a. Domains, branches, uniformity (answers to the first review)

* `E_F(z) = g(n+z) - g(n) - z log S + z log A + z log C - z log D + [log Gamma(B+z) - log Gamma(B)] - [log Gamma(A+z) - ..]
  + [log Gamma(D+z) - ..] - [log Gamma(C+z) - ..]` (l. 1554-1561), `g(w) = const + c w + h(w) - log Gamma(w + 1/2)`,
  `h = log Mt` the branch on `D(n)` real on the real axis (exists: `Mt != 0` on the convex disk `D(n)`, A3), `log Gamma`
  principal on `Re w > 0`. For `|z| <= n/20` (certified `R0 = d + lambda r <= n/20`): `n+z in D(n)`, and
  `Re(A+z), Re(B+z), Re(C+z), Re(D+z), Re(b+z) >= D - R0 > 0`. So `E_F` is holomorphic on the disk `|z| < n/20`,
  which contains the convex set `Omega` and every `conv{0,..,4,z}`, `z in Omega`.
* Precisely: for `|z| < n/20`, `Re(w + z) >= min(A,B,C,D,b) - n/20 >= n(betaJ/2 - 1/20) > 0.4 n` for each shifted
  argument w in {A, B, C, D, b}; and `|n + z - n| < n/20 < n/10`.
* `s := 9/(4 pi e^{2 u_0(n+z)})` (the second variable of `Q_k`): `|s| = (9/4pi) e^{-2 Re u_0(n+z)} <=
  (9/4pi) e^{-2u_0(n) + 2 rho1} = (9/4pi) e^{-1/eps + 2 rho1}`; the certificate uses the complex ball of that radius
  (`sc`), and `t = 1/u_0(n+z)` in the complex ball of radius `2 rho1 (2/L*)^2` around the real cell interval (`tc`).
* Cauchy for `G^{(5)}` is applied on the circle of radius n/20 centred at each point `n + z`, `|z| <= R0 <= n/20`;
  that circle lies in `D(n)`, where `|G - G(n)| <= etaG` (`<= 0.053287`), so `|G^{(5)}(n+z)| <= 120 etaG (20/n)^5` for every such z.
* Scaled polygamma remainders: `n^4 L |B-b| * 240 (n - R0)^{-6} = (bJ+xb) 240/(n(1-R0/n)^6)` and
  `n^4 L (C-D) * 240 (D - R0)^{-6} = xc 240/(n (D/n - R0/n)^6)`; both enter `ce5` (variable `errs`).
* At integers `0 <= j <= d`, `E_F(j)` is real and `exp E_F(j) = R_j/R_j^F > 0`, so `E_F(j) = log(R_j/R_j^F)`
  (real log); matching (Lemma 6.2) gives `E_F(0..4) = 0`. `c_F := exp E_F`, `c_F(d) > 0`.
* Hermite-Genocchi: `E_F(z) = prod_{j<=4}(z-j) * int_{simplex} E_F^{(5)}`, so `|E_F(z)| <= |prod| sup_{conv} |E_F^{(5)}|/5!`.
* Uniformity of the interval evaluation: in each eps-cell the certificate evaluates `n^4 L E_F^{(5)}(z)` with `z/n`
  ranging over the complex ball of radius `R0/n + 1/n` (all of `Omega`), `t = 1/u_0(n+z)` over the complex ball of radius
  `2 rho1 (2/L*)^2` around the real cell interval, and all parameters (eps, bJ, betaJ, xa, xb, xc) over their certified
  enclosures. Ball arithmetic encloses the range over the whole product set, so `ce5` bounds `sup_{Omega}` for every n
  in the cell; the cells cover `0 < eps <= 1/L*`, i.e. all `n >= N*`. A denominator that could vanish would make Arb
  return a non-finite ball and the check would fail; none did.
* The theta-tail term "+1" in rel_k (A2) bounds `[2e^{-kappa1/2}/kappa1 + e^{-Psi(x)} int_{2 u_0(x)}^oo e^{g_1}] *
  e^{x log2/u_0(x)} / sqrt(2pi/|f''|)`, which is `<= e^{-x}` since `g_1(2u_0(x)) - Psi(x) <= 2x log2 + 5u_0 + x/u_0 + 3 -
  x^2/(pi u_0^2)`.

## 9. Assembly (l. 1812-1853) and K

Need `sup|E_F| <= thetaE := log(1 + lambda^4(lambda-1)) - 0.001` (`>= 0.929840`): holds if `n^3 L^2 >= (Cfin/thetaE)^2 d^5`;
since `log(n+2) <= (1+lamL) L`, `lamL = (log(pi L*/2) + 2/N*)/L*`, it suffices that
`n^3 log^2(n+2) >= Kres d^5`, `Kres = (Cfin/thetaE)^2 (1+lamL)^2`. For `d > d0` the wedge with `K >= KN :=
N*^3 log^2(N*+2)/d0^5` forces `n > N*` and `d/n <= d0/N*` (used as `omega^2` in all "tiny" terms).
`K := max(Kres, KN(9e24))`. Certificate (L* = 90, K6 = 2.6, lambda = 1.4, 1000 cells, `certify.log`, all values rounded
outward): `1 + lamL <= 1.05502`, `Kres <= 6416.59`, `KN(d0) <= 42.27`, `KN(9e24) <= 784.13`, so `K <= 6416.59 <= 7000`.
With L* = 80 (`certify_self.log`): `1 + lamL <= 1.06043`, `Kres <= 6764.23 <= 7000`, `KSELF <= 2.430e114`.
Independent check (`sanity3.py`): at n = 1e43, d = 6.75e25 on the wedge boundary, the true max |E_F| on the boundary of
Omega is 0.784, our certified bound there is <= 0.8887 (Cfin <= 70.600), the requirement is <= 0.92984.

## 10. What is not re-proved here (cited, as in Holland)

MSS Thm 1.13 (largest root of multiplicative convolution), MMP Props 2.7, 2.11, 2.17 (real-rootedness, interlacing,
log-mesh), Szego's Jacobi matrix, Polya's integral; our d <= d0 record (COROLLARIES) for part (ii).
