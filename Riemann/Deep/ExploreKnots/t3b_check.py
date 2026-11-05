# t3b_check.py — (i) first off-line zeros of the scalar-FE (+) DH-type combination G of the T(2,5) vector, t < TMAX;
# (ii) argument-principle counts of zeros in [0.51, 2.6] x [0.5, T] for G and H (cross-check of the grid search);
# (iii) independent check of reported zeros with mpmath.dirichlet.
import numpy as np, mpmath as mp, sys, json, time
exec(open('t3_zeros.py').read().split("res = {}")[0].replace("TMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 200.0", "TMAX = 60.0"))
T0 = time.time()
TM = float(sys.argv[1]) if len(sys.argv) > 1 else 400.0
def argcount(f, a, s1, s2, t1, t2, h=0.004):
    # boundary of rectangle counterclockwise
    segs = [np.linspace(s1, s2, int((s2-s1)/h)+2) + 1j*t1,
            s2 + 1j*np.linspace(t1, t2, int((t2-t1)/h)+2),
            np.linspace(s2, s1, int((s2-s1)/h)+2) + 1j*t2,
            s1 + 1j*np.linspace(t2, t1, int((t2-t1)/h)+2)]
    tot = 0.0
    for sg in segs:
        v = np.concatenate([L_np(sg[i:i+4000], f, a) for i in range(0, len(sg), 4000)])
        ph = np.unwrap(np.angle(v)); tot += ph[-1]-ph[0]
        if np.max(np.abs(np.diff(ph))) > 1.0: print('  warning: large phase step', np.max(np.abs(np.diff(ph))))
    return tot/(2*np.pi)
for name in ['C2_T25_H_minus', 'B_T25_Im_psi20', 'D_T34', 'E_T27']:
    f, a = fam[name]
    print(name, 'arg-principle count in [0.51,2.6]x[0.5,60]:', round(argcount(f, a, 0.51, 2.6, 0.5, 60.0), 3), flush=True)
# (i) G_plus to TM
f, a = fam['C_T25_G_plus']
sig = np.arange(0.51, 1.60, 0.02); tt = np.arange(0.5, TM, 0.04)
Z = np.empty((len(sig), len(tt)))
for i, sg in enumerate(sig):
    for j0 in range(0, len(tt), 4000):
        Z[i, j0:j0+4000] = np.abs(L_np(sg + 1j*tt[j0:j0+4000], f, a))
cand = []
for i in range(1, len(sig)-1):
    row = Z[i]
    m = (row[1:-1] < row[:-2]) & (row[1:-1] < row[2:]) & (row[1:-1] < Z[i-1, 1:-1]) & (row[1:-1] < Z[i+1, 1:-1])
    for j in np.nonzero(m)[0] + 1:
        if Z[i, j] < 0.2: cand.append((sig[i], tt[j]))
zs = []
for (sg, t) in cand:
    try: z = mp.findroot(lambda s: L_mp(s, f, a), mp.mpc(sg, t), tol=1e-20, maxsteps=50)
    except Exception: continue
    if abs(L_mp(z, f, a)) < 1e-15 and z.real > 0.5 + 1e-8 and 0 < z.imag < TM:
        zz = (round(float(z.real), 10), round(float(z.imag), 8))
        if all(abs(zz[0]-w[0]) + abs(zz[1]-w[1]) > 1e-6 for w in zs): zs.append(zz)
zs.sort(key=lambda w: w[1])
print(f'G_plus: {len(zs)} zeros with Re s > 1/2 for t < {TM:.0f}; first:', zs[:8], f'[{time.time()-T0:.0f}s]', flush=True)
if zs:
    t1 = zs[0][1]
    print('G_plus arg-principle count in [0.51,2.6]x[0.5,%.1f]:' % (t1+1), round(argcount(f, a, 0.51, 2.6, 0.5, t1+1.0), 3))
# (iii) independent check with mp.dirichlet
mp.mp.dps = 40
chk = [('C2_T25_H_minus', (1.1952706486, 5.53931464)), ('B_T25_Im_psi20', (0.8385444513, 8.292146)), ('D_T34', (0.9863068133, 6.98229592)), ('E_T27', (0.8201734848, 9.76668958))]
if zs: chk.append(('C_T25_G_plus', zs[0]))
for name, (x, y) in chk:
    f, a = fam[name]
    co = [mp.mpf(float(np.real(a[n % f]))) for n in range(f)]
    g = lambda s: mp.dirichlet(s, co)
    z = mp.findroot(g, mp.mpc(x, y), tol=1e-30)
    print(f'{name}: zero {mp.nstr(z, 15)}  |value| {mp.nstr(abs(g(z)), 3)}')
