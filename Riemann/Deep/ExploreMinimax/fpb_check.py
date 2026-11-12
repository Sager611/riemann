#!/usr/bin/env python3
"""explore36 / fpb_check: exact global Le Cam pair F_{p,b}(s) = zeta(s) (1 + b p^{-s} + p^{1-2s}).
A(s) = 1 + b p^{-s} + p^{1-2s} = (1 - al p^{-s})(1 - be p^{-s}), al + be = -b, al be = p.  p^s A(s) is symmetric
under s -> 1 - s, so F has zeta's gamma factor, conductor p^2, root number +1, Euler product, a_n = 1 and
Lambda_F(n) = Lambda(n) for n < p.  Zeros of A: p^{s-1/2} = w, w + 1/w = -b/sqrt(p).
Checks: (1) real parts of A's zeros; (2) Poisson identity sum_{A zeros} g^(gamma) = log p [2 g(0) + sum_{k!=0}
(w^k + w^-k) g(k log p)] (zero side vs prime side) for a C^6 test g, support below and above log p;
(3) Lambda_F(p^2) sign; (4) two-bump Weil form Q_A(f) < 0 iff |b| > 2 sqrt p on a window just above p."""
import math, cmath, functools
import numpy as np
print = functools.partial(print, flush=True)

def zeros_A(p, b, mmax):
    disc = cmath.sqrt((b / math.sqrt(p)) ** 2 - 4)
    ws = [(-b / math.sqrt(p) + disc) / 2, (-b / math.sqrt(p) - disc) / 2]
    out = []
    lp = math.log(p)
    for w in ws:
        lw = cmath.log(w)
        for m in range(-mmax, mmax + 1):
            s = 0.5 + (lw + 2j * math.pi * m) / lp
            out.append(s)
    return ws, out

def g_fun(y, Lg):
    u = np.clip(np.abs(y) / Lg, 0, 1)
    return np.where(np.abs(y) < Lg, (1 - u ** 2) ** 7, 0.0)

_X, _W = np.polynomial.legendre.leggauss(2000)

def ghat_vec(zs, Lg):
    # g^(z) = int g(y) e^{i z y} dy for complex z (entire); Gauss-Legendre on [-Lg, Lg], vectorised in chunks
    y = Lg * _X; w = Lg * _W * g_fun(Lg * _X, Lg)
    tot = 0.0 + 0.0j
    zs = np.asarray(zs)
    for i in range(0, len(zs), 500):
        tot += np.sum(np.exp(1j * np.outer(zs[i:i + 500], y)) @ w)
    return tot

def main():
    for p, b in [(7, 5.0), (7, 7.0), (11, 8.0), (11, 11.5)]:
        ws, zs = zeros_A(p, b, 0)
        re = sorted(set(round(z.real, 6) for z in zs))
        al_be_sq = b * b - 2 * p
        lam_p2 = math.log(p) * (1 - al_be_sq)
        on = abs(b) <= 2 * math.sqrt(p)
        print(f"p={p} b={b}: |b| {'<=' if on else '>'} 2sqrt(p)={2*math.sqrt(p):.3f}; |w|={[round(abs(w),4) for w in ws]}; "
              f"Re(zeros of A) = {re}; Lambda_F(p^2) = {lam_p2:.3f}")
        lp = math.log(p)
        for Lg in [0.9 * lp, 1.5 * lp]:
            _, zs = zeros_A(p, b, 250)
            zero_side = ghat_vec([(s - 0.5) / 1j for s in zs], Lg)
            prime_side = 0.0
            for w in ws:
                for k in range(-3, 4):
                    prime_side += lp * (w ** k) * g_fun(np.array([k * lp]), Lg)[0]
            # sum over both roots of w^{k}: equals log p * sum_k (w1^k + w2^k) g(k log p), k = 0 term 2 g(0)
            print(f"   support {Lg/lp:.1f} log p: zero side {zero_side.real:+.6f}{zero_side.imag:+.1e}i  prime side {prime_side.real:+.6f}")
        # two-bump test on window c = p e^{eps}: f = bump at -lp/2 and +lp/2 (width eps), Q_A(f)/||f||^2
        for eps in [0.05]:
            # g = f*f~ for two equal narrow bumps: g(0) = 2B, g(+-lp) = B (B = bump energy); only k = 0, +-1 inside support
            B = 1.0
            wsum = (ws[0] + 1 / ws[0]).real          # = -b/sqrt(p)
            QA = lp * (2 * (2 * B) + wsum * (B + B))  # log p [2 g(0) + (w + 1/w)(g(lp) + g(-lp))] with g(0)=2B
            print(f"   two-bump window c = p e^{eps}: Q_A/||f||^2 = {QA/(2*B):+.4f} log-units (negative iff |b| > 2 sqrt p)")

if __name__ == '__main__':
    main()
