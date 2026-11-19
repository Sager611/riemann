# JerbyAudit: independent (non-Arb) sampled check of the zero-free Gram blocks and zero counts on [g_(n-3), g_(n+3)].
# Own Riemann-Siegel (main sum + C0 term) with phases reduced in high precision at a centre t_c and expanded to second order in
# delta = t - t_c; grid step 1e-4; minimum of |Z| on the block refined with mpmath siegelz. Not a proof (the Arb cover is); a cross-check.
import sys, time, numpy as np
import mpmath as mp
mp.mp.dps = 30
PI = mp.pi
def gram(n):
    x0 = 2*PI*mp.exp(1 + mp.lambertw((8*n + 1)/(8*mp.e)).real)
    return mp.findroot(lambda t: mp.siegeltheta(t) - n*PI, x0)
def setup(tc):
    a = mp.sqrt(tc/(2*PI)); N = int(mp.floor(a)); th = mp.siegeltheta(tc); thp = mp.diff(mp.siegeltheta, tc)
    k = np.arange(1, N + 1, dtype=np.float64)
    c = np.array([float(mp.fmod(th - tc*mp.log(j), 2*PI)) for j in range(1, N + 1)])
    return dict(N=N, k=k, c=c, w=1/np.sqrt(k), fr=float(thp) - np.log(k), th2=float(1/(2*tc)), tc=tc)
def Zgrid(S, d):   # d: numpy array of offsets
    out = np.empty_like(d)
    for i in range(0, len(d), 200):
        dd = d[i:i+200, None]
        ph = S['c'][None, :] + dd*S['fr'][None, :] + 0.5*S['th2']*dd*dd
        out[i:i+200] = 2*np.sum(S['w'][None, :]*np.cos(ph), axis=1)
    t = float(S['tc']) + d; a = np.sqrt(t/(2*np.pi)); N = S['N']; p = a - N
    assert np.all((p > 0) & (p < 1)), 'N changes inside the window'
    den = np.cos(2*np.pi*p); den = np.where(np.abs(den) < 1e-9, 1e-9, den)
    C0 = np.cos(2*np.pi*(p*p - p - 1/16))/den
    return out + (-1)**(N - 1)*a**-0.5*C0
for n in [int(x) for x in sys.argv[1:]]:
    t0 = time.time(); gs = {j: gram(j) for j in range(n - 3, n + 4)}
    tc = gs[n]; S = setup(tc)
    lo, hi = float(gs[n - 3] - tc), float(gs[n + 3] - tc); h = 1e-4
    d = np.arange(lo, hi + h, h); Z = Zgrid(S, d)
    # spot check against mpmath
    chk = max(abs(Z[i] - float(mp.siegelz(tc + mp.mpf(d[i])))) for i in range(0, len(d), len(d)//7))
    sc = np.nonzero(np.sign(Z[1:]) != np.sign(Z[:-1]))[0]
    zeros = d[sc] + h/2
    cnt = []
    for j in range(n - 3, n + 3):
        a_, b_ = float(gs[j] - tc), float(gs[j + 1] - tc); cnt.append(int(np.sum((zeros >= a_) & (zeros < b_))))
    blo, bhi = float(gs[n - 1] - tc), float(gs[n + 1] - tc); m = (d >= blo) & (d <= bhi)
    Zb = np.abs(Z[m]); ib = int(np.argmin(Zb)); dmin = d[m][ib]
    where = 'first grid point inside the block (|Z| rises inward from g_(n-1))' if ib == 0 else ('last grid point inside the block (|Z| rises inward from g_(n+1))' if ib == len(Zb) - 1 else 'interior grid point')
    zm = mp.siegelz(tc + mp.mpf(dmin))
    endmin = min(abs(mp.siegelz(gs[n - 1])), abs(mp.siegelz(gs[n + 1])))
    print(f'n={n}: grid {len(d)} pts, step {h}, max |Z_grid - siegelz| at 8 spot points = {chk:.1e}; zeros per Gram interval [g_j, g_j+1), j = n-3..n+2: {cnt} (total {sum(cnt)} vs 6 expected on average)')
    print(f'   block [g_(n-1), g_(n+1)]: sign changes {int(np.sum((zeros >= blo) & (zeros < bhi)))}; grid min |Z| = {abs(Z[m][ib]):.6f} (siegelz there {mp.nstr(abs(zm), 8)}) at g_n{dmin:+.5f}, the {where}; |Z| at the block ends >= {mp.nstr(endmin, 6)}; sign {"+" if zm > 0 else "-"}  [{time.time() - t0:.1f} s]', flush=True)
