#!/usr/bin/env python3
"""Finite model for the signed class inequality H(sigma) (BeyondB2 REPORT sec. 3).

Rows: all nonprincipal Dirichlet characters chi mod a prime q (proxy for the sextic rows chi_u at U = q).
Exact L-values for all chi at once: L(s,chi_j) = q^{-s} sum_k e(jk/(q-1)) zeta(s, g^k/q)  (one FFT).
Per row:  Lnum = L(1/2-theta_0, chi),  Q = P^{-1/2} sum_{P<p<=2P} conj chi(p) w(p/P)  (P = q^0.21),
          c = 1/L(0.75, eta conj chi)  (eta = chi_7, generic target; eta = 1 also run),  T = Lnum Q c.
Classes of size K: random; top |Q|; top |c| (smallest |L(0.75)|, 'zero-like'); top |c||Q| (contact-like); top |T|;
adversarial = quarter-arc (phase-aligned) subset of the contact-like class.
r = |sum_B T| / sum_B |T|;  sigma_eff = -ln r / ln K  (H(sigma) asks r <= K^{-sigma}; random phases give 1/2).
Run: uv run --no-project --with numpy python finite_model.py   (single thread, ~15 s, ~230 MB at q = 786433)"""
import numpy as np, math, time, resource, sys
t0 = time.time(); rng = np.random.default_rng(20261105)
BERN = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510]
def hurwitz(s, x, N=12):
    tot = np.zeros_like(x)
    for n in range(N): tot += (n + x) ** (-s)
    y = N + x; tot += y ** (1 - s) / (s - 1) + 0.5 * y ** (-s)
    fact = 1.0; poch = s; ypow = y ** (-s - 1)
    for k, b in enumerate(BERN, start=1):
        fact *= (2 * k - 1) * (2 * k) if k > 1 else 2
        tot += b / fact * poch * ypow; poch *= (s + 2 * k - 1) * (s + 2 * k); ypow = ypow / (y * y)
    return tot
def proot(q):
    f, m, p = [], q - 1, 2
    while p * p <= m:
        if m % p == 0:
            f.append(p)
            while m % p == 0: m //= p
        p += 1
    if m > 1: f.append(m)
    return next(g for g in range(2, q) if all(pow(g, (q - 1) // p, q) != 1 for p in f))
def isp(n): return n > 1 and all(n % i for i in range(2, int(n ** 0.5) + 1))
def rows(q, th=0.04736, sa=0.75, eta=7, zP=0.21):
    g = proot(q); n = q - 1; ak = np.empty(n, dtype=np.int64); a = 1
    for k in range(n): ak[k] = a; a = a * g % q
    x = ak / q
    Lw = q ** (-(0.5 - th)) * n * np.fft.ifft(hurwitz(0.5 - th, x)); Ls = q ** (-sa) * n * np.fft.ifft(hurwitz(sa, x))
    j = np.arange(n); c = 1.0 / Ls[(eta - j) % n]
    P = q ** zP; ka = np.empty(q, dtype=np.int64); ka[ak] = np.arange(n); arr = np.zeros(n)
    for p in range(int(P) + 1, int(2 * P) + 1):
        if isp(p) and p != q: arr[ka[p]] += math.exp(-1.0 / (p / P - 0.999))
    Q = n * np.fft.ifft(arr / math.sqrt(P))[(-j) % n]
    msk = (j != 0) & (((eta - j) % n) != 0)
    return Lw[msk], Q[msk], c[msk]
for q in [65537, 163841, 786433]:
    for eta in [7, 0]:
        Lw, Q, c = rows(q, eta=eta); T = Lw * Q * c
        rr = lambda idx: abs(T[idx].sum()) / np.abs(T[idx]).sum()
        print("q=%d eta=%s  |family mean T|/E|T| = %.3f" % (q, "chi_7" if eta else "1", abs(T.mean()) / np.abs(T).mean()))
        for nm, key in [("top |T|", np.abs(T)), ("top |c||Q| (contact-like)", np.abs(c) * np.abs(Q)), ("top |c| (zero-like)", np.abs(c)),
                        ("top |Q|", np.abs(Q)), ("random (median of 100)", None)]:
            line = []; Ks = []; lr = []
            for R in [0.4, 0.5, 0.6266, 0.7, 0.8]:
                K = int(round(len(T) ** R))
                if key is None: r = float(np.median([rr(rng.choice(len(T), K, replace=False)) for _ in range(100)]))
                else: r = rr(np.argsort(-key)[:K])
                Ks.append(math.log(K)); lr.append(math.log(r)); line.append("%.3f" % r)
            sl = -np.polyfit(Ks, lr, 1)[0]
            print("   %-27s r at K=q^(.4,.5,.6266,.7,.8): %s   fitted sigma %.3f" % (nm, " ".join(line), sl))
        K = int(round(len(T) ** 0.6266)); o = np.argsort(-(np.abs(c) * np.abs(Q)))[:K]; ph = np.angle(T[o])
        best = max((o[np.abs(np.angle(np.exp(1j * (ph - p0)))) < math.pi / 4] for p0 in np.linspace(-math.pi, math.pi, 73)), key=len)
        print("   adversarial quarter-arc subset of the contact-like class: %d of %d rows, r = %.3f" % (len(best), K, rr(best)))
        sys.stdout.flush()
print("time %.1fs, maxrss %.0f MB" % (time.time() - t0, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20))
