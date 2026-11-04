# check_residuals.py: interpolation residuals at the 51 finite nodes and origin values H1(0), H2(0)
import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
import core
core.setup(40)
b, h, B = core.consts()
PL = core.laurent_P()
tab = core.read_table()
(L1, C1), (L2, C2) = core.lists_from_table(tab)
t0 = time.time()
M = {}
for g in (48, 96):
    M[(1,g)] = core.Measure(L1, C1, PL, gl_n=g)
    M[(2,g)] = core.Measure(L2, C2, PL, gl_n=g)
print('build', round(time.time()-t0,1), 's')
for s in [0, 1, 7.5, 40, 100, 160]:
    s = mp.mpf(s)
    for i in (1,2):
        a = core.K_derivs(M[(i,48)], s, 1); c = core.K_derivs(M[(i,96)], s, 1)
        print(' K%d(%s)=%s  GL48-GL96: %s %s' % (i, s, mp.nstr(c[0], 15), mp.nstr(a[0]-c[0], 3), mp.nstr(a[1]-c[1], 3)))
Mi = {1: M[(1,96)], 2: M[(2,96)]}
# residuals at nodes in I_f
worst = [0, 0]
rows = []
for m in core.nodes(100):
    Q, D = core.node_QD(m, PL)
    for i in (1, 2):
        K0, K1 = core.K_derivs(Mi[3-i], mp.mpf(m), 1)
        L = L1 if i == 1 else L2
        c, d = L[m]
        rc = c + K0 / Q
        rd = d + (1 if (i == 1 and m == 1) else 0) + (K1 - D * K0) / Q
        worst[0] = max(worst[0], abs(rc)); worst[1] = max(worst[1], abs(rd))
        rows.append((m, i, rc, rd))
print('max |c residual| =', mp.nstr(worst[0], 5), ' max |d residual| =', mp.nstr(worst[1], 5))
for m, i, rc, rd in rows:
    if m <= 4 or abs(rc) > 2e-10 or abs(rd) > 2e-10:
        print('  m=%d i=%d rc=%s rd=%s' % (m, i, mp.nstr(rc, 4), mp.nstr(rd, 4)))
# origin values
p1_0 = core.p_derivs(Mi[1], mp.mpf(0), 0)[0]; p2_0 = core.p_derivs(Mi[2], mp.mpf(0), 0)[0]
K1_0 = core.K_derivs(Mi[1], mp.mpf(0), 0)[0]; K2_0 = core.K_derivs(Mi[2], mp.mpf(0), 0)[0]
H1_0 = p1_0 + K2_0; H2_0 = p2_0 + K1_0
Q1 = core.node_QD(1, PL)[0]
target = 6 * Q1 * mp.exp(-mp.pi * h)
print('p1(0)=%s (pi^2/3=%s)  p2(0)=%s' % (mp.nstr(p1_0, 15), mp.nstr(mp.pi**2/3, 15), mp.nstr(p2_0, 6)))
print('H1(0)=%s  H2(0)=%s  6 Q1 e^{-pi h}=%s' % (mp.nstr(H1_0, 15), mp.nstr(H2_0, 15), mp.nstr(target, 15)))
print('H1(0)-H2(0)=%s  H2(0)-target=%s' % (mp.nstr(H1_0-H2_0, 4), mp.nstr(H2_0-target, 4)))
pickle.dump({'note': 'measures not pickled'}, open('/dev/null', 'wb'))

