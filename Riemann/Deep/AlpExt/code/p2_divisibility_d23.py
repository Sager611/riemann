# alpext P2: inert-prime divisibility (remark item 2; alpcore Test A; AlpToy item 2) for the CM curve with j = j((1+sqrt-23)/2),
# CM by O_F, F = Q(sqrt-23), defined over K3 = Q(alpha), alpha^3 = alpha + 1 (O_K3 = Z[alpha], disc -23), the real cubic
# subfield of the Hilbert class field H. Model y^2 = 4x^3 - g2 x - g3, g2 = 12 j(j-1728), g3 = 8 j(j-1728)^2,
# Delta = 1728^2 j^2 (j-1728)^3. D^a(wp - z^-2)(0) = a! c_k, a = 2k-2 (DLMF 23.9.7 recursion), exact in K3.
# Valuations at every prime P | p of K3 (degree 1, 2 or 3; unramified for p != 23) via Hensel-lifted factors of x^3-x-1.
# Primes with v_P(Delta) = 6s > 0 (j0 = 0 or 1728 mod P) use the quadratic twist by lam, v_P(lam) = -s (good reduction).
# Prediction: p inert in F (Frob = transposition: P1 of degree 1 and P2 of degree 2, both supersingular):
#   v_P(a! c_k) = a p/(p^2-1) + O(log a) (= a/p + a/p^3 + ...); p split in F (ordinary): v_P = O(log a).
import sys, mpmath as mp
from flint import fmpq_poly, fmpz_poly, fmpq, fmpz
from math import log
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 150
mp.mp.dps = 60
H = [mp.mpf(12771880859375), mp.mpf(-5151296875), mp.mpf(3491750), mp.mpf(1)]
jr = mp.polyroots(H[::-1], maxsteps=200, extraprec=200)
ar = mp.polyroots([1, 0, -1, -1], maxsteps=200, extraprec=200)
areal = [a for a in ar if abs(mp.im(a)) < 1e-40][0]; acx = [a for a in ar if mp.im(a) > 1e-40][0]
jreal = [j for j in jr if abs(mp.im(j)) < 1e-30][0]; jcx = [j for j in jr if abs(mp.im(j)) > 1e-30]
sol = None
for jc in jcx:
    M = mp.matrix([[1, areal, areal ** 2], [1, acx, acx ** 2], [1, mp.conj(acx), mp.conj(acx) ** 2]])
    v = mp.lu_solve(M, mp.matrix([jreal, jc, mp.conj(jc)]))
    co = [mp.re(v[i]) for i in range(3)]
    if all(abs(c - mp.nint(c)) < 1e-20 for c in co) and all(abs(mp.im(v[i])) < 1e-20 for i in range(3)):
        sol = [int(mp.nint(c)) for c in co]
print('j0 = %d + %d*alpha + %d*alpha^2' % tuple(sol))
fK = fmpq_poly([-1, -1, 0, 1])
def red(u): return u % fK
J = fmpq_poly(sol)
J1728 = J - 1728
g2 = red(12 * J * J1728); g3 = red(8 * J * red(J1728 * J1728))
c = {2: g2 / 20, 3: g3 / 28}
for k in range(4, KMAX + 1):
    s = fmpq_poly([0])
    for m in range(2, k - 1): s += c[m] * c[k - m]
    c[k] = red(s) * fmpq(3, (2 * k + 1) * (k - 3))
def vp_int(n, p):
    n = abs(int(n))
    if n == 0: return 10 ** 9
    v = 0
    while n % p == 0: n //= p; v += 1
    return v
def hensel_root(r, p, Mexp):
    pm = p ** Mexp; R = r
    for _ in range(12):
        fv = (R ** 3 - R - 1) % pm; dv = (3 * R * R - 1) % pm
        R = (R - fv * pow(dv, -1, pm)) % pm
    assert (R ** 3 - R - 1) % pm == 0
    return R
def primes_over(p, Mexp):
    pm = p ** Mexp
    roots = [r for r in range(p) if (r ** 3 - r - 1) % p == 0]
    out = []
    if len(roots) == 3:
        for r in roots: out.append((1, [-hensel_root(r, p, Mexp) % pm, 1]))
    elif len(roots) == 1:
        R = hensel_root(roots[0], p, Mexp)
        out.append((1, [-R % pm, 1]))
        out.append((2, [(R * R - 1) % pm, R % pm, 1]))      # x^3-x-1 = (x-R)(x^2+Rx+R^2-1)
    else:
        out.append((3, [-1 % pm, -1 % pm, 0, 1]))
    return out, pm
def vP(u, h, pm, p):  # valuation of the K3-element u (fmpq_poly) at the prime given by the monic factor h mod p^M
    den = int(u.denom()); num = [int(x) for x in u.numer().coeffs()]
    # reduce num mod h over Z/p^M (h monic)
    num = [x % pm for x in num] + [0] * 3
    d = len(h) - 1
    for i in range(len(num) - 1, d - 1, -1):
        q = num[i]
        if q:
            for t in range(d + 1): num[i - d + t] = (num[i - d + t] - q * h[t]) % pm
    vnum = min(vp_int(x, p) if x else 10 ** 9 for x in num[:d])
    return vnum - vp_int(den, p)
def vfact(a, p):
    v = 0; q = p
    while q <= a: v += a // q; q *= p
    return v
def leg(a, p): r = pow(a % p, (p - 1) // 2, p); return -1 if r == p - 1 else r
Mexp = 600
print(f'KMAX={KMAX} (a <= {2*KMAX-2}); precision p^{Mexp}')
for p in (5, 7, 11, 17, 19, 37, 43, 13, 29, 31, 59, 101):
    pr, pm = primes_over(p, Mexp)
    kind = 'inert in F' if leg(-23, p) == -1 else 'split in F'
    for (deg, h) in pr:
        vD = vP(red(J * J * red(J1728 * J1728 * J1728)), h, pm, p)
        sh = 0
        if vD != 0:
            # quadratic twist by lam with v_P(lam) = -vD/6 (still CM by O_F, j = j0, over K3): c_k -> lam^k c_k
            v2 = vP(g2, h, pm, p); v3 = vP(g3, h, pm, p)
            if vD % 6 or v2 - 2 * (vD // 6) < 0 or v3 - 3 * (vD // 6) < 0:
                print(f'p={p} deg={deg}: v_P(Delta) = {vD}, no good twist over K3, skipped'); continue
            sh = vD // 6
            print(f'p={p} deg={deg}: v_P(Delta)={vD}, v(g2)={v2}, v(g3)={v3}: twisted model, v(c_k) shifted by -{sh}k')
        rows = []
        for k in range(2, KMAX + 1):
            a = 2 * k - 2
            V = vfact(a, p) + vP(c[k], h, pm, p) - sh * k
            rows.append((a, V))
        if kind == 'inert in F':
            dev = [(V - a * p / (p * p - 1)) / max(1.0, log(a) / log(p)) for a, V in rows]
        else:
            dev = [V / max(1.0, log(a) / log(p)) for a, V in rows]
        a_last, V_last = rows[-1]
        print(f'p={p:3d} ({kind}, P of degree {deg}): a={a_last}: v={V_last}, a p/(p^2-1)={a_last*p/(p*p-1):.1f}, a/(p-1)={a_last/(p-1):.1f}; '
              f'(v - pred)/log_p a in [{min(dev):+.2f}, {max(dev):+.2f}] (pred = a p/(p^2-1) if inert else 0)')
