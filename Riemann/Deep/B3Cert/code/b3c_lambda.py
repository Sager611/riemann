# b3cert: (1) explicit constants of Lemma B for a concrete cutoff; (2) independent Lambda-route evaluation of T_sigma(n) (the proof's own
# coordinates: inner m(rho), outer van der Corput range I_sigma) compared with the FFT/xi-route value of b3c_fals.py; (3) the proof's
# intermediate bounds |m|, Var m, |J| on N_sigma tested numerically.
import sys, math, numpy as np
sys.path.insert(0, 'Riemann/Deep/B3Cert/code')
from b3c_fals import setup, Tpm, apm
def S(t):
    t = np.clip(t, 1e-12, 1 - 1e-12); a = np.exp(-1 / t); b = np.exp(-1 / (1 - t)); return a / (a + b)
def chit(u): return np.where((u > 0.5) & (u < 2), S((u - 0.5) * 6) * S((2 - u) * 2), 0.0)
ra, rb = 1 / (4 * math.sqrt(2)), math.sqrt(2)
def alpha_fn(u, T, X0, d, sig):   # chunked (memory cap)
    a = np.empty(len(u), complex)
    for i in range(0, len(u), 4000):
        ap, am = apm(X0 * np.sqrt(u[i:i+4000]), d); a[i:i+4000] = ap if sig > 0 else am
    return chit(u) * (T * u / (2 * np.pi)) ** (-d / 2) * u ** -0.25 * a
def norms(T, X0, d, sig):
    s = np.linspace(2 ** -0.5, 2 ** 0.5, 200001); ds = s[1] - s[0]
    bt = 2 * s * alpha_fn(s ** 2, T, X0, d, sig); out = []
    g = bt
    for k in range(1, 4):
        g = np.gradient(g, ds); out.append(np.sum(np.abs(g)) * ds)
    return out   # B1, B2, B3 = L1 norms of beta~', '', '''
def B2_ibp(T, X0, d, sig):
    u = np.linspace(0.5, 2, 150001); du = u[1] - u[0]; al = alpha_fn(u, T, X0, d, sig); best = 0
    for rho in np.concatenate([np.linspace(-60, sig * ra, 300) if sig > 0 else np.linspace(-60, -rb, 300), np.linspace(sig * rb, 60, 300) if sig > 0 else np.linspace(-ra, 60, 300)]):
        if sig > 0 and ra < rho < rb: continue
        if sig < 0 and -rb < rho < -ra: continue
        ph1 = sig / (2 * np.sqrt(u)) - rho
        g = np.gradient(np.gradient(al / ph1, du) / ph1, du)
        best = max(best, (1 + abs(rho)) ** 2 * np.sum(np.abs(g)) * du)
    return best
def m_of_rho(rho, X0, beta_t, s, ds, sig):
    s0 = sig / (2 * rho); return np.sum(beta_t * np.exp(-1j * X0 * rho * (s - s0) ** 2)) * ds
if __name__ == '__main__':
    T = float(sys.argv[1]); th = float(sys.argv[2]); al, be = complex(sys.argv[3]), complex(sys.argv[4])
    H = T ** th; K = T / H; d = al + be
    w, Phi, N = setup(T, H, be)
    for r in [float(x) for x in sys.argv[5].split(',')]:
        c_r = 2 * math.sqrt(2 * math.pi) * math.sqrt(r); X0 = K * c_r; npq = r * T / H ** 2
        tp, tm = Tpm(T, H, npq, al, be, w, Phi)
        for sig, tx in [(1, tp), (-1, tm)]:
            B1, B2, B3 = norms(T, X0, d, sig); B2i = B2_ibp(T, X0, d, sig)
            M0 = 8 * (2 * ra) ** -0.5 * B1; M1 = (rb - ra) * 8 * (2 * ra) ** -1.5 * (B1 + B2 / ra + B3 / (2 * X0 * ra))
            C1 = (8 * (8 * math.sqrt(2)) ** 0.5 * (3 * M0 + M1) + 2 * B2i / X0) / (2 * math.pi * math.sqrt(2))
            # Lambda route on I_sigma
            ds = 0.04 / X0; s = np.arange(2 ** -0.5, 2 ** 0.5, ds); bt = 2 * s * alpha_fn(s ** 2, T, X0, d, sig)
            rg = sig * np.linspace(ra, rb, 801); mv = np.array([m_of_rho(x, X0, bt, s, ds, sig) for x in rg])
            dm = np.abs(np.gradient(mv, rg)); varm = np.sum(dm[:-1] * np.abs(np.diff(rg)))
            dL = 0.03 / K; Lam = np.arange(min(sig * ra, sig * rb) * c_r, max(sig * ra, sig * rb) * c_r, dL)
            o = np.argsort(rg); mi = np.interp(Lam / c_r, rg[o], mv[o].real) + 1j * np.interp(Lam / c_r, rg[o], mv[o].imag)
            y = 0.5 * H * np.log1p(Lam / H); f = 0.5 * np.exp(-y ** 2 - (1 + 2 * be) * y / H)
            psi = K * H * np.log1p(Lam / H) + K * c_r ** 2 / (4 * Lam)
            JS = np.sum(f * mi * np.exp(1j * psi)) * dL
            kap = 1j * sig * math.sqrt(math.pi / 2) * np.exp(-1j * sig * math.pi / 4)
            tl = kap / math.sqrt(math.pi) * (T / (2 * math.pi)) * X0 ** -0.5 * JS
            # |J| on N_sigma (sample) vs B2 X0^-2 (1+|rho|)^-2
            u = np.linspace(0.5, 2, int(40 * X0 * 3) + 1); duu = u[1] - u[0]; alu = alpha_fn(u, T, X0, d, sig); jr = 0
            for rho in ([-3, -1, 0.0, 0.1, 1.6, 3.0] if sig > 0 else [-3, -1.6, -0.1, 0.0, 1, 3]):
                Jv = np.sum(alu * np.exp(1j * X0 * (sig * np.sqrt(u) - rho * u))) * duu; jr = max(jr, abs(Jv) * X0 ** 2 * (1 + abs(rho)) ** 2)
            sc = H * X0 ** -0.5
            print(f"T={T:.0e} th={th} r={r} X0={X0:.4g} sigma={sig:+d}: |T_sig| xi-route {abs(tx)/sc:.5f}, Lambda-route {abs(tl)/sc:.5f}, rel.diff {abs(tl-tx)/abs(tx):.1e}; "
                  f"sup|m|X0^1/2={np.abs(mv).max()*X0**0.5:.3f} (M0={M0:.1f}), Var m X0^1/2={varm*X0**0.5:.3f} (M1={M1:.1f}), max|J|X0^2(1+|rho|)^2 on N={jr:.2e} (B2={B2i:.3g}); "
                  f"explicit C1={C1:.3g} vs observed {abs(tx)/sc:.3f}", flush=True)
