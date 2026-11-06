# t2_cauer.py (explore22): Cauer-II RC ladder of the 'zero network' Z(p) = sum_rho 1/(w_rho + p), w_rho = z_rho^2,
# built ONLY from Taylor data at the centre s = 1/2 (no zeros used), for xi, mv (xi with one moved pair) and dh.
# Z is a passive RC impedance (Stieltjes) iff every w_rho > 0 iff RH for that function (Grommer 1914).
# Its Stieltjes continued fraction c1/(1+c2 p/(1+c3 p/(...))) has c_{2n} = p1_n/p0_n, c_{2n+1} = p0_{n+1}/p1_n,
# p0_n, p1_n the LU pivots of the Hankel matrices [m_{i+j}], [m_{i+j+1}], m_k = sigma_{k+1} = sum_rho w_rho^{-(k+1)}.
# Usage: t2_cauer.py <xi|mv|dh> <dps> <nmax> [hfac]
import sys, time, json
import mpmath as mp

kind = sys.argv[1]; MVD = '0.15'
if kind.startswith('mv:'): kind, MVD = 'mv', kind[3:]
dps = int(sys.argv[2]); nmax = int(sys.argv[3]); hfac = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
K = 2 * nmax + 3          # number of power sums needed
mp.mp.dps = dps + 20
t0 = time.time()

def moments_xi(K):
    # M_{2n} = int_0^inf Phi(u) u^{2n} du, Phi(u) = sum_m (2 pi^2 m^4 e^{9u} - 3 pi m^2 e^{5u}) exp(-pi m^2 e^{4u})
    # G(z) = int_0^inf Phi cos(zu) du = (1/8) Xi(z/2): zeros at z = 2 gamma.  Strip |Im u| < pi/8, use d = pi/16.
    nats = mp.mpf(mp.mp.dps) * mp.log(10) + 40
    h = (2 * mp.pi * (mp.pi / 16) / nats) * hfac
    U = mp.mpf(0.25) * mp.log((nats + 9 * 2) / mp.pi + 2 * K * 1.0) + mp.mpf('0.15')
    J = int(U / h) + 1
    M = [mp.mpf(0)] * K
    for j in range(J + 1):
        u = j * h
        e4 = mp.exp(4 * u); e5 = mp.exp(5 * u); e9 = mp.exp(9 * u)
        phi = mp.mpf(0); m = 1
        while True:
            ex = mp.exp(-mp.pi * m * m * e4)
            term = (2 * mp.pi ** 2 * m ** 4 * e9 - 3 * mp.pi * m * m * e5) * ex
            phi += term
            if ex < mp.mpf(10) ** (-(mp.mp.dps + 30)) / (1 + e9 * m ** 4): break
            m += 1
        wgt = h / 2 if j == 0 else h
        u2 = u * u; p = mp.mpf(1)
        for n in range(K):
            if n > 0 and j == 0: break
            M[n] += wgt * phi * p
            p *= u2
    # F(w) = G(2 sqrt w) = sum_n (-1)^n 4^n w^n M_{2n}/(2n)!  -> zeros at w = gamma^2
    return [(-1) ** n * mp.mpf(4) ** n * M[n] / mp.factorial(2 * n) for n in range(K)], J, h

def moments_dh(K):
    # Lambda(1/2+iz) = 2 int_0^inf Psi(v) cos(zv) dv, Psi(v) = 2 e^{3v/2} theta(e^{2v}), theta(x) = sum a_n n e^{-pi n^2 x/5}
    # Psi even (functional equation). Strip |Im v| < pi/4, use d = pi/8.
    kap = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
    A = {1: mp.mpf(1), 2: kap, 3: -kap, 4: mp.mpf(-1), 0: mp.mpf(0)}
    nats = mp.mpf(mp.mp.dps) * mp.log(10) + 40
    h = (2 * mp.pi * (mp.pi / 8) / nats) * hfac
    U = mp.mpf(0.5) * mp.log(5 * (nats + 2 * K * 2.0) / mp.pi) + mp.mpf('0.2')
    J = int(U / h) + 1
    M = [mp.mpf(0)] * K
    for j in range(J + 1):
        v = j * h
        x = mp.exp(2 * v)
        th = mp.mpf(0); n = 1
        while True:
            ex = mp.exp(-mp.pi * n * n * x / 5)
            th += A[n % 5] * n * ex
            if ex < mp.mpf(10) ** (-(mp.mp.dps + 30)) / n: break
            n += 1
        psi = 2 * mp.exp(1.5 * v) * th
        wgt = h / 2 if j == 0 else h
        v2 = v * v; p = mp.mpf(1)
        for k in range(K):
            if k > 0 and j == 0: break
            M[k] += wgt * psi * p
            p *= v2
    return [(-1) ** n * 2 * M[n] / mp.factorial(2 * n) for n in range(K)], J, h

def power_sums(f, K):
    g = [c / f[0] for c in f]
    s = [mp.mpf(0)] * (K)
    for n in range(1, K):
        s[n] = -n * g[n] - mp.fsum(s[j] * g[n - j] for j in range(1, n))
    return s   # s[k] = sum_rho w_rho^{-k}

if kind in ('xi', 'mv'):
    f, J, h = moments_xi(K + 1)
else:
    f, J, h = moments_dh(K + 1)
tq = time.time() - t0
sig = power_sums(f, K + 1)
if kind == 'mv':
    with mp.workdps(mp.mp.dps):
        ga = mp.zetazero(6).imag; gb = mp.zetazero(7).imag
        gbar = (ga + gb) / 2; dl = mp.mpf(MVD)
        za = mp.mpc(gbar, -dl); zb = mp.mpc(gbar, dl)
        for k in range(1, K + 1):
            sig[k] = sig[k] - ga ** (-2 * k) - gb ** (-2 * k) + mp.re(za ** (-2 * k) + zb ** (-2 * k))

m = [sig[k + 1] for k in range(2 * nmax + 2)]
def lu_pivots(shift):
    A = [[m[i + j + shift] for j in range(nmax)] for i in range(nmax)]
    piv = []
    for k in range(nmax):
        p = A[k][k]; piv.append(p)
        for i in range(k + 1, nmax):
            fac = A[i][k] / p
            Ai = A[i]; Ak = A[k]
            for j in range(k + 1, nmax):
                Ai[j] -= fac * Ak[j]
    return piv
p0 = lu_pivots(0); p1 = lu_pivots(1)
c = [m[0]]
for n in range(1, nmax + 1):
    c.append(p1[n - 1] / p0[n - 1])          # c_{2n}
    if n < nmax: c.append(p0[n] / p1[n - 1])  # c_{2n+1}
first_neg = next((k + 1 for k, ck in enumerate(c) if ck < 0), None)
out = {'kind': kind, 'mv_delta': MVD if kind == 'mv' else None, 'dps': dps, 'nmax': nmax, 'hfac': hfac, 'nodes': J, 'h': float(h), 't_quad_s': round(tq, 1),
       't_total_s': round(time.time() - t0, 1), 'sigma1': mp.nstr(sig[1], 25), 'sigma2': mp.nstr(sig[2], 20),
       'first_negative_c_index': first_neg,
       'c_sign': ''.join('+' if ck > 0 else '-' for ck in c),
       'log10_c': [round(float(mp.log10(abs(ck))), 3) for ck in c],
       'p0_digits': [mp.nstr(p, 12) for p in p0[:]], 'p1_digits': [mp.nstr(p, 12) for p in p1[:]]}
print(json.dumps(out), flush=True)
