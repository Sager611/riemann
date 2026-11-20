# ceilaudit7 A: zeta check of D(a,b) = sinc(a-b) - sinc(a)sinc(b) on actual zeros (normalised shifts 2*pi*a/log(g/2pi))
import sys, numpy as np, mpmath as mp
mp.mp.dps = 15
z = np.load('/Users/adriansager/playground/riemann/Riemann/Deep/ExploreHelio/results/zeta_zeros.npy')
z = np.sort(z[z > float(sys.argv[1])])[:int(sys.argv[2])]
print(f'{len(z)} zeros in [{z[0]:.1f}, {z[-1]:.1f}]')
pairs = [(0.5, 0.5), (1.0, 1.0), (0.25, -0.25), (0.5, -0.5), (0.9, -0.9), (1.2, -1.2)]
shifts = sorted(set([p[0] for p in pairs] + [p[1] for p in pairs]))
vals = {}
for s in shifts:
    print('shift', s, flush=True)
    vals[s] = np.array([float(mp.siegelz(g + 2*np.pi*s/np.log(g/(2*np.pi)))) for g in z])
Lm = np.mean(np.log(z/(2*np.pi))); Z2 = Lm + 2*0.5772156649
print(f'mean Z^2 = {Z2:.4f} (L + 2*gammaE = {Lm + 2*0.5772:.4f}); zero density check: n/(range*L/2pi) = {len(z)/((z[-1]-z[0])*Lm/(2*np.pi)):.4f}')
for (a, b) in pairs:
    emp = np.mean(vals[a]*vals[b])/Z2
    cue = np.sinc(a-b) - np.sinc(a)*np.sinc(b)
    print(f'D({a:+.2f},{b:+.2f}) zeta={emp:+.4f} CUE-limit={cue:+.4f} ratio={emp/cue if abs(cue)>1e-9 else float("nan"):.3f}')
