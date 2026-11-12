# explore45 derivations (appends of <= 60 lines)

## D1. Setup and the profile LP (16:55)
f = V*Psi analytic in a strip, u = log|f| subharmonic. Rescale s = 1/2 + x/L + i t, L = log T.
m(x) = lim (1/T) int_T^{2T} u dt (mean log), C(x) = mean |f|^2, h(x) = (1/2) log C(x).
Littlewood: m(x) - m(x') = 2 pi int_x^{x'} n(v) dv (x < x'), n(v) = rescaled density of zeros with Re > v.
Hence m is convex and decreasing; n(0) = right-zero density (Levinson counts zeros of V with beta >= 1/2).
Normalisation: leading Dirichlet coefficient 1, so m(+inf) = 0 and m >= 0.
Jensen on each line: m(x) <= h(x).
LP: maximise s := -m'(0+)... precisely n(0) subject to {m convex, decreasing, m(+inf)=0, m <= h}.
Claim D1: LP value = inf_{R>0} h(-R)/R (times 1/(2 pi)) = Levinson's bound kappa >= 1 - (1/R) log c(R).
 Proof (upper): m(-R) >= m(0) + R*2pi n(0) >= R*2pi n(0) (convexity/monotone n, m(0) >= 0), and m(-R) <= h(-R).
 Proof (attained): m*(x) = max(0, -s* x), s* = inf_R h(-R)/R; m* convex, decreasing, m*(+inf)=0,
 m* <= h: for x<0 since -s*x = s*|x| <= h(x) by definition of s*; for x >= 0 since h >= 0 (C >= 1).
Key point: u*(x,t) := m*(x) (independent of t) is itself SUBHARMONIC (convex in x), with
 mean m*, E e^{2u*} = e^{2 m*} <= C on every line, Riesz measure (1/2pi) Delta u* = (s*/2pi) dt on x = 0
 (smeared "zeros" on the line, density = Levinson's bound). So every inequality valid for all subharmonic
 functions (star function, Steiner/circular symmetrisation, Polya-Szego, harmonic-measure comparison,
 Carleman, Hadamard-Doetsch three lines) is satisfied by u*, which meets Levinson's bound with equality.
 => no such inequality + linewise second moments can beat Levinson. Escape: only the quantisation of
 the Riesz measure (atoms of mass 2 pi at zeros, finite total density) separates log|f| from u*.

## D2. Finite-T measurement (17:50; .scratch/explore45/levin.py, ana.py)
f = V Psi, V = F + (b/L)F', Psi = sum_{n<=y} F^{-1}(n) n^{a-1/2} (log(y/n)/log y) n^{-s}, theta=1/2, b=1.03217, R=1.24618
(asymptotic optimum for P(x)=x, linear Q: kappa >= 0.347439 = Levinson), T in [1e4, 1.4e4], dt = 0.01, L = log(T_mid/2pi) = 7.555.
Euler-Maclaurin (K=12) zeta checked vs Arb to 1e-10; DH zero 0.808517+85.699348i reproduced (|f| = 6.5e-7).
Counts by winding on [1/2,3] x [T0,T1]: zeta: V 0 right zeros, Psi 0 (under RH, b < 2 forces V's zeros left of 1/2).
Finite-T Levinson bound nu_V(0) <= h_f(-R)/(R kn) = 0.2491 (kappa-equiv 0.502); truth nu_V(0) = 0.
 decomposition: Jensen/concavity 0.1437 + Littlewood weights of V's zeros in (a,1/2) 0.1056 (mollifier zeros 0, m_f(0) ~ 0).
 LP over all 23 lines (only h_f used): 0.2426 at x = -1 (= best single line, as D1 predicts).
 Jensen gap J_f(x) = h_f - m_f: 0.179 at x=-R, 0.164 at -1, 0.096 at 0, -> 0 to the right; J/var = 0.77 at -R (not Gaussian).
 Without mollifier J_V(-R) = 0.734: the mollifier is a Jensen-gap reducer.
DH (same code, a_n = 1, k, -k, -1, 0, k = 0.28408; mollifier = truncated Dirichlet inverse): V has 275 right zeros (nu 0.057),
 bound 0.534, concavity 0.197, weights 0.261, m_f(0) term 0.020. (DH normalised by zeta's count; its own L would be log(5T/2pi).)

## D3. Quantised inverse problem (17:46; .scratch/explore45/invlev.py)
Question: D1's extremal u* is a smeared (non-atomic) Riesz measure. Do genuine log|f| (atoms 2pi, bounded density) with the
same line second moments also reach the LP? Competitor: f(z) = p(w), w = e^{-z/Pz} (period 2pi Pz in y = L(t - t0), Pz = zeta
zeros per period), p = prod_{j<=K} (1 - w/w_j), c_0 = 1 (leading coefficient 1, f -> 1 as x -> +inf), K = 1.61 Pz (= total zero
density of V Psi: 1 + log y / L). Line x <-> circle |w| = e^{-x/Pz}; mean|f|^2 on the line = sum |c_k|^2 r^{2k}, exact with
2K+2 nodes. Constraint: mean|f|^2 <= C(x) = measured mean|V Psi|^2 on every measured line. Maximise M = #{|w_j| <= 1}
(zeros with Re >= 1/2). L-BFGS on a soft-max of log S - log C over roots (u_j = -+a_j^2, v_j), several starts.
Profile 23 lines x in [-4, 18.9] (T in [1e4, 1.4e4]): Pz=30: M=6 feasible (max violation -4e-4), M=7 infeasible (+0.024);
 Pz=60: M=12 feasible, M=14 infeasible (+0.024). LP bound 0.2424 Pz.
Dense profile 42 lines x in [-12, 18.9] (T in [1e4, 1.2e4], LP 0.2407): Pz=60: M=8,10,12 feasible, M=13 +0.0049 (near),
 M=14 +0.028. Feasible competitors put all M right zeros ON the line (x = 0) and a screening layer at x ~ -2.7.
=> an analytic f with the same 42 line second moments has right-zero fraction 0.200 (true V Psi: 0); smeared LP 0.2407;
 single-line Levinson 0.2478. Quantisation can buy at most 0.2407 - 0.200 = 0.041 (17% of the bound) at this T.

## D4. Larger period and independent check (19:08)
Pz = 100 (K = 161, 42 lines): M = 21 +0.0011 (unconverged at the 1500 s cap), 22 +0.0092, 23 +0.023 (LP 24.1): nu_max ~ 0.20-0.21.
verify.py rebuilds the Pz = 60, M = 12 competitor with np.poly: np.roots gives 12 roots in |w| <= 1, max(h_comp - h) = -1.5e-4
on all 42 lines; m_comp = 0.200|x| on [-1.4, 0]; J_comp = 0.030 at x = -1 (true V Psi 0.164), 0.087 at x = 0.
Rigorous: no argument from these line second moments + analyticity + c_0 = 1 + bandwidth 1.61 proves nu < 0.200.
Numerical only: quantised optimum in [0.200, ~0.21] < smeared 0.2407.
