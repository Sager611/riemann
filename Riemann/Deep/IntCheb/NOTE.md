# Integer Chebyshev constant of [0,1]: rigorous re-check and a new lower bound

## Object
t_Z([0,1]) = lim_n (min over nonzero integer polynomials P of degree <= n of sup_[0,1] |P|)^(1/n).
arXiv 2610.00386 (D. Niedbala Giraudin, Oct 2026, v3) proves 0.4222286000 <= t_Z([0,1]) <= 0.4226846976.
Their upper bound comes from an explicit integer polynomial (92 factors in w = x(1-x), degree 2*10^10).
Their lower bound comes from Theorem "cert": for a probability measure mu on [0,1], pairwise distinct
irreducible primitive Q_j with leading coefficient a_j >= 1, c_j >= 0, B = sum c_j d_j <= 1,
t_Z >= exp(min{inf_[0,1] A, -sum c_j log a_j, min_i G_i}), where A = U_mu - sum c_j log|Q_j| and
G_i = L_i(1-B+c_i d_i) + int U_mu dnu_i - sum_{j!=i} c_j int log|Q_j| dnu_i.

## Proof check of Theorem cert (by hand): correct
- Res(Q, Q_j) is a nonzero integer for coprime factors, so d_j L_Q + int log|Q_j| dnu_Q >= 0.
- A is superharmonic off [0,1]; if B < 1 then A -> +inf at infinity, so inf over C = inf over [0,1].
- Useful identity: G_i = (1/d_i)[int log|Q_i| dmu - sum_{j!=i} c_j log|Res(Q_i,Q_j)|] (uses mu(total)=1).
- If the total mass M is not exactly 1 (floating point), the same proof gives
  t_Z >= exp(min(inf A, min G_i)/M) with the G_i formula above, provided B < M. We use this exact form.

## Rigorous verification method (ours)
- lowcert.c: bisection for inf A. On I = [t-H, t+H]: U_mu = sum_k beta_k F(x - e_k), F(u) = u log|u| - u,
  beta_k = density jumps. Edges and real roots within 2H of t: exact monotone bounds (F decreases on [-1,1];
  log|x - r| <= log(|t - r| + H)). All other terms: second-order Taylor bound
  f(t) - |f'(t)| H + min(0, inf_I f'') H^2/2, with inf_I f'' by ball evaluation on I. All arithmetic in arb balls;
  roots from arb_fmpz_poly_complex_roots (rigorous enclosures).
- gcert.c: G_i in arb with exact integer resultants (fmpz_poly_resultant) and U_mu at complex roots through
  Re(u log u - u) with atan2 for the argument.
- upcert.c: the same bisection for the upper certificate's f(w) = sum (e_j/N) log|R_j(w)| on [0,1/4].
- verify.py: exact Fraction checks (edges, masses, multipliers, M, B, exact x -> 1-x symmetry), irreducibility
  (python-flint or sympy), then the C programs, then exp by mpmath interval arithmetic.

## Results
| item | published | ours (rigorous) |
|---|---|---|
| lower, paper's data | 0.422228600019 | inf A >= V with margin 1.2e-11 (27 522 intervals vs their 12.3 M); min G = -0.86220828217; t_Z >= 0.42222860001915 |
| upper, paper's data | 0.422684697505 | max f <= -1.72225754796595 (margin 2.6e-10, 5 752 intervals); t_Z <= 0.422684697505 |
| new lower | - | t_Z >= 0.42223549904149 |

## How the new certificate was found
- Symmetric LP (mu and c invariant under x -> 1-x; 198 library orbits), variables = cell masses, multipliers, V.
  Rows: A(x) >= V at sample points, G_i >= V, B <= 1 - 1e-9, total mass 1.
- HiGHS (scipy) needed 340-480 MB on this dense LP, so we wrote a dense primal-dual interior-point method
  (Mehrotra; memory = matrix + n^2). Cutting planes: after each solve, A is scanned on 200 000 points plus
  golden-section refinement of local minima; violated minima are added, slack points pruned.
- LP optimum on the paper's 947 cells: V = -0.8622082833 (paper: -0.8622084054).
- Cells re-graded: the paper's 23 support intervals on [0,1/2], each with twice as many Chebyshev-graded cells,
  gaps Chebyshev-graded with width <= 0.002, edges rounded to multiples of 2^-50 so that 1 - e is exact.
  K = 1832: V_lp = -0.86219116, fine-grid min of A = -0.862192064, certificate V = -0.862192066.
- Multipliers below 1e-10 set to 0 (17 remain); this only lowers B and raises every G_i, and the final A check
  is rigorous anyway.
- Binding divisibility terms: the polynomials x and x-1 and 2x-1 (min G_i = -0.86219116325).

## Assessment
A strict, certified improvement of the lower bound by 6.9e-6 (1.6e-5 in log). It uses the paper's method and
library; the gain comes from a better-optimized and finer measure. The paper says more cells and more polynomials
should help, so this is incremental. Convergence: 947 -> 1832 cells gained 1.7e-5 in V; the limit for this library is
probably near V = -0.86217 to -0.86219 (t ~ 0.42224).

## Not done
- Upper bound: no new factors searched.
- Library growth: candidates via dual reduced costs RC(Q) = d_Q (lambda + sum tau_i L_i) + sum_p w_p log|Q(p)|
  (points = dual-weighted sample points and roots of tau-weighted library polynomials), searched with LLL.
- Spectral basis: on each support interval, density p(t)/sqrt((t-a)(b-t)) with Chebyshev p; potentials are explicit
  (-T_n/n inside, -zeta^(-n)/n outside), so the LP gets much smaller.
