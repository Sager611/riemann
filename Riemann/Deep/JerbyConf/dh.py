# dh.py -- jerbyconf: recompute the two Davenport-Heilbronn numerics of Jerby 2610.21952 (r_* ~ 0.72346, l. 778-783; mu_DH(g_44) ~ 1.45, l. 840-846)
# D(s) = sum c_k k^-s, c = (1, kappa, -kappa, -1, 0) periodic mod 5 (= ((1-i kappa)/2) L(s,chi) + ((1+i kappa)/2) L(s,chibar), chi(2) = i).
# theta_DH(t) = Im log Gamma(3/4 + i t/2) - (t/2) log(pi/5); Z_DH(t) = e^{i theta_DH} D(1/2+it) (real); xi_DH(1/2+it) = |G| Z_DH, G = (pi/5)^{-s/2} Gamma((1+s)/2).
# Section (by analogy with the paper's Z_N(t;a), l. 448-460): S_N(t;r) = cos theta_DH(t) + r sum_{k=2}^{N} c_k k^-1/2 cos(theta_DH(t) - t log k).
import mpmath as mp, sys, json
mp.mp.dps = 30
s5 = mp.sqrt(5); kap = (mp.sqrt(10 - 2*s5) - 2)/(s5 - 1)
c = [mp.mpf(0), mp.mpf(1), kap, -kap, mp.mpf(-1)]
def theta(t): return mp.im(mp.loggamma(mp.mpf(3)/4 + 1j*t/2)) - t/2*mp.log(mp.pi/5)
def D(s): return 5**(-s)*sum(c[a]*mp.zeta(s, mp.mpf(a)/5) for a in range(1, 5))
def Zdh(t): return mp.re(mp.exp(1j*theta(t))*D(mp.mpf(1)/2 + 1j*t))
def Zdh_im(t): return mp.im(mp.exp(1j*theta(t))*D(mp.mpf(1)/2 + 1j*t))
def xi(t):
    s = mp.mpf(1)/2 + 1j*t
    return (mp.pi/5)**(-s/2)*mp.gamma((1 + s)/2)*D(s)
def gram(n, t0):
    return mp.findroot(lambda t: theta(t) - n*mp.pi, t0)
out = {}
g44_data = mp.mpf('85.5848588208579388842977101588649287134')
g44 = gram(44, g44_data)
out['g44'] = mp.nstr(g44, 25); out['g44_minus_dataset'] = mp.nstr(g44 - g44_data, 5)
out['theta(g2_dataset)/pi'] = mp.nstr(theta(mp.mpf('10.4953565915728235541830418434035855369'))/mp.pi, 20)
Z = Zdh(g44); Zp = mp.diff(Zdh, g44); out['Z_DH(g44)'] = mp.nstr(Z, 15); out['Im check'] = mp.nstr(Zdh_im(g44), 5)
out['Z_DH prime(g44)'] = mp.nstr(Zp, 15); out['mu_normalised'] = mp.nstr(abs(Zp/Z), 12)
X = xi(g44); Xp = mp.diff(lambda t: mp.re(xi(t)), g44)
out['xi(g44)'] = mp.nstr(X, 15); out['xi prime(g44)'] = mp.nstr(Xp, 15); out['mu_xi'] = mp.nstr(abs(Xp/mp.re(X)), 12)
out['dlog|G|/dt'] = mp.nstr(mp.diff(lambda t: mp.log(abs((mp.pi/5)**(-(mp.mpf(1)/2 + 1j*t)/2)*mp.gamma((mp.mpf(3)/2 + 1j*t)/2))), g44), 12)
for k in (43, 45):
    gk = gram(k, g44 + (k - 44)*1.49); out[f'g{k}'] = mp.nstr(gk, 20); out[f'Z_DH(g{k})'] = mp.nstr(Zdh(gk), 12)
# sections
def S(t, r, N, d=0):
    th = theta(t)
    f = lambda tt: mp.cos(theta(tt)) + r*sum(c[k % 5]/mp.sqrt(k)*mp.cos(theta(tt) - tt*mp.log(k)) for k in range(2, N + 1))
    return f(t) if d == 0 else mp.diff(f, t, d)
mp.mp.dps = 20
sec = []
for N in range(2, int(sys.argv[1]) if len(sys.argv) > 1 else 31):
    v = S(g44, 1, N); vp = S(g44, 1, N, 1)
    Gl = mp.diff(lambda t: mp.log(abs(mp.gamma((mp.mpf(3)/2 + 1j*t)/2))), g44)
    # continued extremum g(r): S'(g;r) = 0 from g(0) = g44; Delta(r) = S(g(r); r)
    rs = [mp.mpf(i)/200 for i in range(201)]; g = g44; Dl = []
    ok = True
    for r in rs:
        try: g = mp.findroot(lambda t: S(t, r, N, 1), g)
        except Exception: ok = False; break
        Dl.append((r, g, S(g, r, N)))
    rstar = None
    for (r0, g0, d0), (r1, g1, d1) in zip(Dl, Dl[1:]):
        if d0*d1 < 0:
            gg = [g0]
            def Dr(rr):
                gg[0] = mp.findroot(lambda t: S(t, rr, N, 1), gg[0]); return S(gg[0], rr, N)
            rstar = mp.findroot(Dr, (r0, r1), solver='anderson'); break
    sec.append({'N': N, 'S(g44;1)': mp.nstr(v, 10), 'mu_sec': mp.nstr(abs(vp/v), 8), 'mu_sec_unnorm': mp.nstr(abs(vp/v + Gl), 8),
                'r_star': (mp.nstr(rstar, 8) if rstar is not None else None), 'track_ok': ok, 'Delta(1)': mp.nstr(Dl[-1][2], 8) if Dl else None})
    print(sec[-1]); sys.stdout.flush()
out['sections'] = sec
json.dump(out, open(sys.argv[2] if len(sys.argv) > 2 else 'dh_out.json', 'w'), indent=1)
for k, v in out.items():
    if k != 'sections': print(k, v)
