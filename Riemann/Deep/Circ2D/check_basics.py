# check_basics.py: P_l, Q_n, D_n, spectral measure == direct P*R, node identities, table norms
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
import core
core.setup(40)
b, h, B = core.consts()
PL = core.laurent_P()
paper = [5, -1+2j*b, 1+2j*b, -2, -mp.mpf(1)/2+1j*b, -1-2j*b, 1]
print('P_l check (max |mine-paper|):', max(abs(PL[l]-paper[l]) for l in range(7)), ' conj sym:', max(abs(PL[-l]-mp.conj(PL[l])) for l in range(7)))
# P identity vs product
for s in [mp.mpf('0.3'), mp.mpf('5.7'), mp.mpf('17.25')]:
    print(' P prod vs Laurent', s, abs(core.Pval(s) - core.P_deriv(s,0,PL).real))
sq3 = mp.sqrt(3); pi = mp.pi
QDpaper = {0:(pi**2/3, -7*pi*sq3/18), 1:(2*pi**2/3*(2-sq3), pi*(3+sq3)/18), 3:(2*pi**2/3*(2-sq3), -pi*(3+sq3)/18),
           4:(pi**2/3, 7*pi*sq3/18), 7:(2*pi**2/3*(2+sq3), -pi*(3-sq3)/18), 9:(2*pi**2/3*(2+sq3), pi*(3-sq3)/18)}
for a in A if False else core.A:
    for n in [a, a+12, a+24]:
        Q, D = core.node_QD(n, PL)
        print(' node', n, 'Q', mp.nstr(Q, 12), 'dQ', mp.nstr(Q-QDpaper[a][0], 3), 'D', mp.nstr(D, 12), 'dD', mp.nstr(D-QDpaper[a][1], 3))
tab = core.read_table()
(L1, C1), (L2, C2) = core.lists_from_table(tab)
print('table nodes', len(tab), min(tab), max(tab))
for L, nm in [(L1,'1'),(L2,'2')]:
    tot = sum(abs(c)+abs(d) for n,(c,d) in L.items())
    big = sum(abs(c)+abs(d) for n,(c,d) in L.items() if n > 40)
    print(' list', nm, 'norm', mp.nstr(tot, 12), 'norm n>40', mp.nstr(big, 8))
t0 = time.time()
M1 = core.Measure(L1, C1, PL, gl_n=48)
print('measure build s', round(time.time()-t0,1))
for s in [mp.mpf('0.37'), mp.mpf('2.5'), mp.mpf('10.1'), mp.mpf('40.7')]:
    pd = core.p_direct(L1, C1, s)
    pm = core.p_derivs(M1, s, 0)[0]
    print(' p direct vs measure', s, mp.nstr(pd, 15), mp.nstr(pd-pm, 3))
# node identities p(n)=Q c, p'(n)=Q(d+Dc)
for n in [0, 1, 3, 13, 40]:
    Q, D = core.node_QD(n, PL)
    c, d = L1[n]
    v = core.p_derivs(M1, mp.mpf(n), 1)
    print(' node', n, 'p(n)-Qc', mp.nstr(v[0]-Q*c, 3), "p'(n)-Q(d+Dc)", mp.nstr(v[1]-Q*(d+D*c), 3))

