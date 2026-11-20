# z6_mesh.py — explore66, P4 smallest case: forward gap rates gdot_k = 4/g_k - 2 g_k S_k at t = 0 (zero ODE), first 100 gaps of
# zeta (z = 2 gamma; all zeros +-x_j up to n = 400, tail by density) and of Davenport–Heilbronn (real zeros + its off-line pairs, gamma < 140).
# P4: for the globally minimal gap m, S <= 2/m^2 (lattice equality) => gdot >= 0. Test: sign of gdot_k and the lattice ratio S_k g_k^2 / 2.
import mpmath as mp, json, time
mp.mp.dps = 20
T0 = time.time()
LOG = '/Users/adriansager/playground/riemann/Riemann/Deep/ExploreConley/logs/'
xz = [2*mp.zetazero(n).imag for n in range(1, 401)]
print('zeta zeros', round(time.time()-T0, 1), 's', flush=True)
def rates(xs, others, nmax, tailrho=None):
    allz = [x for x in xs] + [-x for x in xs] + others
    out = []
    for k in range(nmax):
        a, b = xs[k], xs[k+1]; g = b - a
        S = sum(1/((b - z)*(a - z)) for z in allz if z != a and z != b)
        if tailrho is not None:   # zeros beyond the list, density rho(x) = log(x/4pi)/(4pi), both signs: int_X^inf 2 rho(x)/x^2 dx approx
            X = xs[-1]; S += 2*tailrho(X)/X
        S = mp.re(S)
        out.append((k+1, g, S, 4/g - 2*g*S, S*g**2/2))
    return out
rz = rates(xz, [], 100, lambda X: mp.log(X/(4*mp.pi))/(4*mp.pi))
kmin = min(rz, key=lambda r: r[1])
print('zeta, first 100 gaps: min gap (z units) at k=%d: g=%s, S=%s, gdot=%s, S g^2/2=%s (lattice = 1)' % (kmin[0], mp.nstr(kmin[1], 8), mp.nstr(kmin[2], 8), mp.nstr(kmin[3], 8), mp.nstr(kmin[4], 6)))
print('zeta: #gdot<0 among 100:', sum(1 for r in rz if r[3] < 0), '; max S g^2/2 =', mp.nstr(max(r[4] for r in rz), 6), '; min S g^2/2 =', mp.nstr(min(r[4] for r in rz), 6))
# Davenport–Heilbronn real zeros: sign changes of the real function Lambda_f(1/2 + i gamma) on a grid
kap = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1); c1 = (1 - 1j*kap)/2; c2 = mp.conj(c1)
chi = [0, 1, 1j, -1j, -1]; chib = [0, 1, -1j, 1j, -1]
def Lam(g):
    s = mp.mpf(1)/2 + 1j*g
    return mp.re((5/mp.pi)**((s+1)/2)*mp.gamma((s+1)/2)*(c1*mp.dirichlet(s, chi) + c2*mp.dirichlet(s, chib)))
grid = [mp.mpf(k)/50 for k in range(1, 50*140)]
vals = [Lam(g) for g in grid]
gz = []
for i in range(len(grid)-1):
    if vals[i]*vals[i+1] < 0:
        gz.append(mp.findroot(Lam, (grid[i], grid[i+1]), solver='anderson'))
print('DH real zeros with gamma < 140:', len(gz), ' time', round(time.time()-T0, 1), flush=True)
offl = [(mp.mpf('0.808517182457'), mp.mpf('85.699348485378')), (mp.mpf('0.650830080610'), mp.mpf('114.163342730757'))]
oth = []
for be, ga in offl:
    for bb in (be, 1 - be):
        z = 2*ga - 2j*(bb - mp.mpf(1)/2); oth += [z, mp.conj(z), -z, -mp.conj(z)]
xd = [2*g for g in gz]
nd = min(100, len(xd) - 2)
rd = rates(xd, oth, nd, lambda X: mp.log(5*X/(4*mp.pi))/(4*mp.pi))
kd = min(rd, key=lambda r: r[1])
print('DH, first %d gaps: min gap at k=%d: g=%s, S=%s, gdot=%s, S g^2/2=%s' % (nd, kd[0], mp.nstr(kd[1], 8), mp.nstr(kd[2], 8), mp.nstr(kd[3], 8), mp.nstr(kd[4], 6)))
print('DH: #gdot<0:', sum(1 for r in rd if r[3] < 0), '; max S g^2/2 =', mp.nstr(max(r[4] for r in rd), 6), '; min =', mp.nstr(min(r[4] for r in rd), 6))
json.dump({'zeta': [[r[0], str(r[1]), str(r[2]), str(r[3])] for r in rz], 'dh': [[r[0], str(r[1]), str(r[2]), str(r[3])] for r in rd], 'dh_real_gamma': [str(g) for g in gz]}, open(LOG + 'z6_mesh.json', 'w'))
print('time', round(time.time()-T0, 1))
