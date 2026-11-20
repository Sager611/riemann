# li_universal.py (explore56): is lambda_n > 0 forced by the change-of-measure structure alone?
# Class: positive symmetric laws. V = log Y; under the Esscher tilt Q_{1/2} (density ∝ Y^{1/2}) V has an even law g.
# Then F(s) = E_P[Y^s] = M(s - 1/2)/M(1/2), M(w) = E_g[e^{wV}], satisfies F(0)=F(1)=1, F(s)=F(1-s) (BPY: F = 2 xi for zeta).
# Li: lambda_n = (1/(n-1)!) d^n/ds^n [s^{n-1} log F(s)] at s=1 = sum_m C(n,m) kt_m/(m-1)!, kt_m = (log M)^{(m)}(1/2)
# (cumulants of V under Q = Y*P, the size-biased / Girsanov-tilted law). kt_1 = KL(Q|P) >= 0, kt_2 = Var_Q(V) >= 0 (BPY 2001).
# Family: g = p*(d_{a}+d_{-a})/2 + (1-p)*(d_{b}+d_{-b})/2, M(w) = p cosh(a w) + (1-p) cosh(b w).
import mpmath as mp, itertools, sys
mp.mp.dps = 30
def cumulants(p, a, b, m):
    M = lambda w: p*mp.cosh(a*w) + (1-p)*mp.cosh(b*w)
    return [mp.diff(lambda w: mp.log(M(w)), mp.mpf(1)/2, k) for k in range(1, m+1)]
def lam(n, kt):
    return sum(mp.binomial(n, m)*kt[m-1]/mp.factorial(m-1) for m in range(1, n+1))
def zeta_check():
    # same formula on zeta: kt_m = (log 2xi)^{(m)}(1); expect lambda_1..3 = 0.0230957, 0.0923457, 0.2076389
    xi = lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
    kt = [mp.diff(lambda s: mp.log(2*xi(s)), 1, k, method='quad', radius=0.25) for k in range(1, 6)]  # Cauchy circle: avoids the removable point s=1
    return [lam(n, kt) for n in range(1, 6)], kt
if __name__ == '__main__':
    lz, ktz = zeta_check()
    print('zeta: lambda_1..5 =', [mp.nstr(x, 10) for x in lz])
    print('zeta: Esscher cumulants kt_1..5 =', [mp.nstr(x, 8) for x in ktz])
    best = {}
    for p in [0.5, 0.8, 0.9, 0.95, 0.99, 0.999]:
        for a in [0.5, 1, 2, 4, 6, 8]:
            for b in [1.5, 3, 6, 9, 12, 18, 24]:
                if b <= a: continue
                kt = cumulants(mp.mpf(p), mp.mpf(a), mp.mpf(b), 10)
                for n in range(1, 11):
                    v = lam(n, kt)
                    if n not in best or v < best[n][0]: best[n] = (v, p, a, b)
    for n in sorted(best):
        v, p, a, b = best[n]
        print('n=%d min lambda_n over grid = %s at p=%g a=%g b=%g' % (n, mp.nstr(v, 8), p, a, b))
