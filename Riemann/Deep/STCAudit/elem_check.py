
# Sampled suprema of F_u, F_l and F_u^S, F_l^S (Prop. 4.9 of Casana-Giner) on [1, 1e4] (geometric grid, ratio 1.001)
# and spot values beyond; compare with Table 2; recompute A3 = ceil_3(Cz + sup).
import json, math, sys
import mpmath as mp
mp.mp.dps = 20
K = json.load(open('Riemann/Deep/STCAudit/kernels.json'))
A1s = {'L5': 0.20176, 'L20': 0.10823, 'L60': 0.09473, 'L150': 0.09349, 'C400': 0.08399, 'D6000': 0.08382}
T2 = {'L5': (0.70360002, 0.66536440, 0.13949090, 1.369, 0.844), 'L20': (1.46807534, 0.70207071, -0.18278649, 2.171, 1.286),
      'L60': (2.10426949, 0.72894156, -0.15847575, 2.834, 1.946), 'L150': (2.25578346, 0.73031685, -0.15757634, 2.987, 2.099),
      'C400': (4.46374163, 0.73756551, -0.14582673, 5.202, 4.318), 'D6000': (4.58951035, 0.73573753, -0.14641036, 5.326, 4.444)}
def M(T): return T/(2*mp.pi)*mp.log(T/(2*mp.pi*mp.e))
def A(s): return mp.atan(s.imag/s.real) + mp.atan(s.imag/(s.real-1)) - s.imag/2*mp.log(mp.pi) + mp.im(mp.loggamma(s/2))
def Fs(k, A1, T):
    T = mp.mpf(T); d = k['d']; s1 = mp.mpc(0.5+d, T); s2 = mp.mpc(0.5+2*d, T)
    E = (2*A(s1) - A(s2))/mp.pi - M(T)
    ku = 0; kl = 0
    for (X, j), c in zip(k['B'], k['c']):
        Dj = (1/s1 + 1/(s1-1) - mp.log(mp.pi)/2 + mp.digamma(s1/2)/2) if j == 1 else (s1**(-j) + (s1-1)**(-j) + (-1)**(j-1)/mp.factorial(j-1)/2**j*mp.psi(j-1, s1/2))
        v = mp.re(Dj) if X == 'R' else -mp.im(Dj)
        ku += c*v; kl += (1 if X == 'R' else -1)*c*v
    Eu = E + ku/mp.pi; El = E - kl/mp.pi
    vp = mp.siegeltheta(T)/mp.pi + 1 - M(T) - mp.mpf(7)/8     # = varpi/pi
    L = mp.log(T)
    return (float(Eu - A1*L), float(-El - A1*L), float(Eu - mp.mpf(7)/8 - vp - A1*L), float(-(El - mp.mpf(7)/8 - vp) - A1*L))
out = {}
Ts = []
T = 1.0
while T <= 1e4:
    Ts.append(T); T *= 1.001
Ts += [2e4, 5e4, 1e5, 1e6, 1e8, 1e12]
for name, k in K.items():
    A1 = A1s[name]
    best = [(-9, 0)]*4
    for T in Ts:
        v = Fs(k, A1, T)
        best = [max(best[i], (v[i], T)) for i in range(4)]
    supN = max(best[0][0], best[1][0]); supS = max(best[2][0], best[3][0])
    cz, tN, tS, a3N, a3S = T2[name]
    rN = math.ceil((cz + tN)*1000 - 1e-12)/1000; rS = math.ceil((cz + tS)*1000 - 1e-12)/1000
    out[name] = dict(sup_Fu=best[0], sup_Fl=best[1], sup_FuS=best[2], sup_FlS=best[3])
    print('%-5s sampled sup max(Fu,Fl)=%.8f (Fu at T=%.4g, Fl at T=%.4g)  Table2 <= %.8f | sampled sup S-mode=%.8f (at T=%.4g/%.4g) Table2 <= %.8f | A3N=ceil(Cz+T2)=%.3f (paper %.3f)  A3S=%.3f (paper %.3f)' % (
        name, supN, best[0][1], best[1][1], tN, supS, best[2][1], best[3][1], tS, rN, a3N, rS, a3S))
    sys.stdout.flush()
json.dump(out, open('Riemann/Deep/STCAudit/elem_check.json', 'w'), indent=1)
