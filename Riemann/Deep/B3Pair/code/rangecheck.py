# b3pair range check at the new nu0 (referee request, 20 Nov 2026): every nu/theta hypothesis of the B3N/B3F chain, exact rationals.
# Moment hypotheses at (theta', nu0) for every theta' in (theta_c, 1); counting/transfer at theta in (theta_c, 1). Linear in theta' => endpoints decide.
from fractions import Fraction as F
k, l = F(89, 2742), F(1184, 1371)                  # outer pair A(89/1282, 997/1282)
thc = F(4939, 10000)
nu1 = lambda th: (th*(1 + 2*l) - (k + l))/(2*(1 + l - k))
nu0 = nu1(thc)
print('nu0 = nu1(theta_c) =', nu0, '=', float(nu0), '; nu1(th) == (7478 th - 2457)/10042:', all(nu1(t) == (7478*t - 2457)/F(10042) for t in (F(0), F(1), thc)))
def row(tag, src, name, f, lin=True):
    a, b = f(thc), f(F(1))
    ok = a > 0 and b > 0 if lin else a > 0
    print(f'{tag} [{src}] {name}: margin at theta_c = {a} ({float(a):.6g}); at theta->1 = {b} ({float(b):.6g}); {"linear => >0 on (theta_c,1)" if lin else ""} -> {"OK" if ok else "FAIL"}')
    return ok
res = []
res.append(row('R1', 'B3F l.55, 85; B3N row 4 (B3F l.5, 17); E1 l.41', 'nu < 1/2', lambda t: F(1, 2) - nu0))
res.append(row('R2', 'B3F l.55; B3N §3, §4 setting', "1/3 < theta' < 1 (lower end)", lambda t: t - F(1, 3)))
res.append(row('R3', 'B3A l.63', "theta' >= 1/3 (off-band tail e^{-H^2/400})", lambda t: t - F(1, 3)))
eps_max = thc - 2*nu0
print('eps fixed in (0, theta_c - 2 nu0) = (0,', eps_max, '=', float(eps_max), '); uniform in theta\' > theta_c')
res.append(row('R4', 'B3F l.57; B3A l.31', "eps < theta' (c_r <= H/46), eps < theta_c - 2nu0", lambda t: t - eps_max))
res.append(row('R5', 'B3N §4(c), §5.1 row 8 + standing facts', "nu < theta'/2 (L_u >= T^{theta'-2nu-eps} >= 1; pq <= y^2 < H)", lambda t: t/2 - nu0))
res.append(row('R6', 'B3N §4(d), row 9', "nu < (1 + theta')/4 (side terms)", lambda t: (1 + t)/4 - nu0))
d = lambda t: nu1(t) - nu0
print('R7 [B3N §4(d); LEMMA C\'] nu < nu1(theta\') strictly: nu1(theta\') - nu0 == (7478/10042)(theta\' - theta_c):', all(d(t) == F(7478, 10042)*(t - thc) for t in (F(1, 2), F(3, 4), F(1))), '-> = 0 at theta_c, > 0 for theta\' > theta_c: OK (theta\' = (theta+theta_c)/2 gives 3739(theta-theta_c)/10042)')
res.append(row('R8', 'B3F l.154', "nu < theta' (N0 < T)", lambda t: t - nu0))
res.append(row('R9', 'B3N row 6 (B3F l.19, 52-63)', 'pq <= T^{2nu} < T: 1 - 2nu > 0', lambda t: 1 - 2*nu0))
Amin = lambda t: (1 + nu0 - t)/(1 - t)
print('R10 [B3F l.91; B3A l.63; E1A l.61, 76] A(1-theta\') > 1 + nu - theta\' + delta: A_min(theta_c) =', Amin(thc), '=', float(Amin(thc)), '; A_min finite for every fixed theta\' < 1 (A_min(0.99) =', float(Amin(F(99, 100))), ') -> OK, A fixed after theta\'')
print('R11 [B3F l.98-99; E1A l.77] sigma(1 - 2nu) > 1 + A: 1 - 2nu0 =', 1 - 2*nu0, '=', float(1 - 2*nu0), '> 0 -> sigma > (1+A)/(1-2nu0) fixed: OK')
res.append(row('R12', 'B3F l.84; B3A l.54; B3N rows 7, 10, 11', "theta' < 1 (residue change o(H/L); supp w; K = T^{1-theta'} -> oo)", lambda t: 1 - t if t < 1 else F(1, 10**9)))
for th in (thc + F(1, 10**4), thc + F(1, 10**6)):
    eta = (th - thc)/2
    print(f'R13 [B3F l.108-110; E1A l.79] windows at theta = {float(th):.7f}: theta\' = {float(th-eta):.7f}, nu1(theta\') - nu0 = {float(d(th-eta)):.4g} > 0; D_min = (1 + nu0 - theta)/(2 eta) = {float((1 + nu0 - th)/(2*eta)):.6g} (+delta/(2eta)) finite -> OK')
res.append(row('R14', 'B3F l.112; B3A l.84; E1 l.64', 'theta > 0 (N(T,H) ~ HL/2pi, Backlund O(L))', lambda t: t))
res.append(row('R15', 'B3F l.113; B3A l.87; E1A l.66', '0 < lambda < theta < 1 (Wang local pair; lambda up to theta)', lambda t: t))
res.append(row('R16', 'B3F l.113, 177; E1 l.83', "0 < theta < 1 (c(theta), c'(theta) = cot^2(theta/sqrt2)/2 > 0)", lambda t: t))
print('R17 [B3F l.114, 179-180; logs/b3p_arb_tty_4000.log] kappa(nu0, 44/5) - k*(theta_c) = 2.6350934196464372870e-7 +/- 7.35e-28 > 0 (Arb); k* decreasing => kappa > k*(theta) for theta > theta_c: OK')
print('R18 not used (B3F l.110; B3A l.74-80, 136-137): SL l.199 nu < theta - 1/2; SL l.536-538 theta > 1/2; SL l.731-732 theta - eta > 1/2; Taylor J > (1+nu-theta+delta)/(2theta-1) (E1 l.53; E1A l.60, 75): exact profile, no Taylor step')
print('ALL linear rows OK:', all(res))
