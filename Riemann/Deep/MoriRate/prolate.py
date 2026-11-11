# prolate.py (morirate audit): prolate spheroidal functions psi_n (n = 0, 4) at bandwidth c by the Legendre (Bouwkamp) tridiagonal
# matrix, refined by Rayleigh-quotient/inverse iteration in mpmath; checks inputs of Mori 24047156 (paper III) and its companion [4], [5]:
#  d_n = 1 - lambda_n(c) against 6621 c^{9/2} e^{-2c} ([4]) and Fuchs 1210 c^{9/2} e^{-2c}; A_n^2 = psi_n(1)^2 <= 4 c d_n/lambda_n ([5, L3.3]);
#  I_0^2 >= pi l lambda/c (paper III L3.2) i.e. beta_0^2 >= pi l/(2c), l = 395/399; inward growth |psi_n(1-x)| <= A_n exp(2c sqrt x) (L3.1).
# lambda_n = c beta_0^2/(pi psi_n(0)^2) from int e^{icxt} psi_n(t) dt = mu_n psi_n(x), |mu_n|^2 = 2 pi lambda_n/c (x = 0).
import sys, math, json
import mpmath as mp

def run(c_mult_pi):
    c = mp.pi*c_mult_pi
    cf = float(c)
    digits = int(-(math.log10(1210) + 4.5*math.log10(cf) - 2*cf/math.log(10))) + 60
    mp.mp.dps = digits
    c = mp.pi*c_mult_pi
    kmax = int(1.6*cf + 2.5*digits) + 20
    if kmax % 2: kmax += 1
    ks = list(range(0, kmax+1, 2))
    m = len(ks)
    D = [k*(k+1) + c**2*mp.mpf(2*k*(k+1) - 1)/((2*k+3)*(2*k-1)) for k in ks]
    E = [c**2*mp.mpf((k+1)*(k+2))/((2*k+3)*mp.sqrt(mp.mpf((2*k+1)*(2*k+5)))) for k in ks[:-1]]
    out = {'c_over_pi': c_mult_pi, 'digits': digits, 'kmax': kmax}
    import numpy as np
    from scipy.linalg import eigh_tridiagonal
    # only the 3 lowest eigenvectors (all of them at c = 1000 pi took 644 MB: memory cap 300 MB)
    if c_mult_pi > 400: raise SystemExit('c/pi > 400 exceeds the memory/time budget')
    w, V = eigh_tridiagonal(np.array([float(x) for x in D]), np.array([float(x) for x in E]), select='i', select_range=(0, 2))
    for n, idx in ((0, 0), (4, 2)):
        chi = mp.mpf(w[idx]); v = [mp.mpf(x) for x in V[:, idx]]
        for it in range(8):
            # solve (T - chi) y = v (Thomas algorithm), normalise; Rayleigh quotient update
            a = [D[i] - chi for i in range(m)]
            cp = [mp.mpf(0)]*m; dp = [mp.mpf(0)]*m
            cp[0] = E[0]/a[0]; dp[0] = v[0]/a[0]
            for i in range(1, m):
                den = a[i] - E[i-1]*cp[i-1]
                cp[i] = E[i]/den if i < m-1 else mp.mpf(0)
                dp[i] = (v[i] - E[i-1]*dp[i-1])/den
            y = [mp.mpf(0)]*m; y[-1] = dp[-1]
            for i in range(m-2, -1, -1): y[i] = dp[i] - cp[i]*y[i+1]
            nr = mp.sqrt(mp.fsum(t*t for t in y)); v = [t/nr for t in y]
            Tv = [D[i]*v[i] + (E[i]*v[i+1] if i < m-1 else 0) + (E[i-1]*v[i-1] if i > 0 else 0) for i in range(m)]
            chi_new = mp.fsum(v[i]*Tv[i] for i in range(m))
            res = mp.sqrt(mp.fsum((Tv[i] - chi_new*v[i])**2 for i in range(m)))
            chi = chi_new
            if res < mp.mpf(10)**(-digits+10): break
        if v[0] < 0: v = [-t for t in v]
        # psi(x) = sum beta_k sqrt(k+1/2) P_k(x)
        P0 = []
        p = mp.mpf(1)
        for j, k in enumerate(ks):
            P0.append(p)
            p = -p*mp.mpf(k+1)/(k+2)
        psi0 = mp.fsum(v[j]*mp.sqrt(ks[j] + mp.mpf(1)/2)*P0[j] for j in range(m))
        psi1 = mp.fsum(v[j]*mp.sqrt(ks[j] + mp.mpf(1)/2) for j in range(m))
        lam = c*v[0]**2/(mp.pi*psi0**2)
        d = 1 - lam
        rec = {'n': n, 'chi_over_c2': float(chi/c**2), 'res': mp.nstr(res, 3), 'd': mp.nstr(d, 8),
               'd_over_6621bound': float(d/(6621*c**mp.mpf(4.5)*mp.e**(-2*c))),
               'd_over_fuchs': float(d/(1210*c**mp.mpf(4.5)*mp.e**(-2*c))),
               'A2_over_4cd': float(psi1**2/(4*c*d/lam)),
               'beta0sq_over_pil2c': float(v[0]**2/(mp.pi*mp.mpf(395)/399/(2*c))) if n == 0 else None}
        # inward growth check on a grid of x in (0, 1]
        worst = 0.0
        for xg in [mp.mpf(10)**(-e) for e in range(1, 9)] + [mp.mpf(j)/20 for j in range(1, 21)]:
            u = 1 - xg
            # Legendre recurrence for P_k(u), k = 0..kmax
            pkm1, pk = mp.mpf(1), u
            vals = {0: mp.mpf(1)}
            for k in range(1, kmax+1):
                if k == 1: vals[1] = u; continue
                pkp = ((2*k-1)*u*vals[k-1] - (k-1)*vals[k-2])/k
                vals[k] = pkp
            ps = mp.fsum(v[j]*mp.sqrt(ks[j] + mp.mpf(1)/2)*vals[ks[j]] for j in range(m))
            ratio = abs(ps)/(abs(psi1)*mp.e**(2*c*mp.sqrt(xg)))
            worst = max(worst, float(ratio))
        rec['inward_worst_ratio'] = worst
        out['n%d' % n] = rec
    print(json.dumps(out), flush=True)

if __name__ == '__main__':
    for cm in [float(a) for a in sys.argv[1:]]:
        run(cm)
