# Numerical checks of the analytic identities quoted in Baluyot's thesis, ch. 3 (mpmath).
# (1) Estermann functional equation with shifts (used in Lemma 3.1.3)
# (2) Lemma 3.1.3 (= Conrey 1989, Lemma 5) at a test point
# (3) Lemma 3.1.7 (= Conrey 1989, Lemma 3) at a test point
# (4) Lemma 3.1.2 bound |D(0,2sigma-1,0,H/K)| <= 2K^2 on a grid
import sys
from math import gcd
import mpmath as mp
mp.mp.dps = 25
SEC = set(sys.argv[1].split(',')) if len(sys.argv) > 1 else {'1', '2', '3', '4'}
e = lambda x: mp.exp(2j * mp.pi * x)


def D(s, al, be, H, K):
    # sum_{m,n} m^{-s-al} n^{-s-be} e(mnH/K) via Hurwitz zeta (analytic continuation)
    return sum(e(mp.mpf(a * b * H) / K) * mp.power(K, -2 * s - al - be)
               * mp.zeta(s + al, mp.mpf(a) / K) * mp.zeta(s + be, mp.mpf(b) / K)
               for a in range(1, K + 1) for b in range(1, K + 1))


def inv(H, K):
    return pow(H, -1, K)

# (1) Estermann: D(s,al,be,H/K) = (1/pi)(K/2pi)^{1-2s-al-be} G(1-s-al)G(1-s-be)
#     {D(1-s,-al,-be,Hb/K) cos(pi(al-be)/2) - D(1-s,-al,-be,-Hb/K) cos(pi(2s+al+be)/2)}
for (s, al, be, H, K) in ([] if '1' not in SEC else [(mp.mpc(0.3, 2.1), 0.2, 0.0, 2, 5), (mp.mpc(-0.4, 0.7), 0.13, -0.07, 3, 7), (mp.mpc(0.5, 9.0), 0.3, 0.1, 1, 4)]):
    Hb = inv(H, K)
    lhs = D(s, al, be, H, K)
    rhs = (1 / mp.pi) * mp.power(K / (2 * mp.pi), 1 - 2 * s - al - be) * mp.gamma(1 - s - al) * mp.gamma(1 - s - be) * (
        D(1 - s, -al, -be, Hb, K) * mp.cos(mp.pi * (al - be) / 2) - D(1 - s, -al, -be, (-Hb) % K, K) * mp.cos(mp.pi * (2 * s + al + be) / 2))
    print("(1) Estermann FE s=%s al=%s be=%s H/K=%d/%d: |lhs-rhs|/|lhs| = %.2e" % (mp.nstr(s, 4), al, be, H, K, abs(lhs - rhs) / abs(lhs)))

# (2) Lemma 3.1.3: S(w,al,be,H/K) = sum m^{-al} n^{-be} e(mnH/K) e(mnw), Im w > 0
def S_direct(w, al, be, H, K, L=400):
    tot = mp.mpc(0)
    for m in range(1, L):
        for n in range(1, L // m + 1):
            tot += mp.power(m, -al) * mp.power(n, -be) * e(mp.mpf(m * n * H) / K) * e(m * n * w)
    return tot


def S_lemma(w, al, be, H, K, c=1.5):
    Hb = inv(H, K)
    z = -2j * mp.pi * w
    t1 = mp.power(K, -1 + al - be) * mp.zeta(1 + be - al) * mp.gamma(1 - al) * mp.power(z, al - 1)
    t2 = mp.power(K, -1 + be - al) * mp.zeta(1 + al - be) * mp.gamma(1 - be) * mp.power(z, be - 1)
    t3 = D(0, al, be, H, K)

    def f(t):
        s = c + 1j * t
        return (mp.power(z, s - 1) * mp.gamma(1 - s) * mp.power(K / (2 * mp.pi), 2 * s - 1 - al - be)
                * mp.gamma(s - al) * mp.gamma(s - be)
                * (D(s, -al, -be, Hb, K) * mp.cos(mp.pi * (al - be) / 2)
                   - D(s, -al, -be, (-Hb) % K, K) * mp.cos(mp.pi * (2 - 2 * s + al + be) / 2)))
    I = mp.quad(f, [-40, -15, -5, 0, 5, 15, 40]) * 1j  # ds = i dt
    t4 = I / (2 * mp.pi ** 2 * 1j)
    return t1 + t2 + t3 + t4

for (w, al, be, H, K) in ([] if '2' not in SEC else [(mp.mpc(0.1, 0.25), 0.2, 0.0, 2, 5), (mp.mpc(-0.2, 0.3), 0.35, 0.0, 1, 3)]):
    a = S_direct(w, al, be, H, K)
    b = S_lemma(w, al, be, H, K)
    print("(2) Lemma 3.1.3 w=%s al=%s H/K=%d/%d: direct=%s lemma=%s rel.diff=%.2e" % (mp.nstr(w, 3), al, H, K, mp.nstr(a, 10), mp.nstr(b, 10), abs(a - b) / abs(a)))


# (3) Lemma 3.1.7: (1/(i Delta sqrt(pi))) int_(c) exp((s-s0)^2/Delta^2) chi(1+beta-s) omega^{-s} ds
#     = omega^{-beta} int_0^inf z^{s0-beta} exp(-Delta^2 (log z)^2/4) (e(-omega z) + e(omega z)) dz/z
def chi(s):
    return mp.power(2, s) * mp.power(mp.pi, s - 1) * mp.sin(mp.pi * s / 2) * mp.gamma(1 - s)

for (Dl, s0, be, om, c) in ([] if '3' not in SEC else [(3.0, mp.mpc(0.6, 1.3), 0.2, 0.7, 1.5), (2.0, mp.mpc(0.55, -0.8), 0.1, 1.9, 1.3)]):
    lhs = mp.quad(lambda t: mp.exp(((c + 1j * t) - s0) ** 2 / Dl ** 2) * chi(1 + be - (c + 1j * t)) * mp.power(om, -(c + 1j * t)),
                  [-30, -10, 0, 10, 30]) * 1j / (1j * Dl * mp.sqrt(mp.pi))
    rhs = mp.power(om, -be) * mp.quad(lambda x: mp.exp(x * (s0 - be)) * mp.exp(-Dl ** 2 * x ** 2 / 4)
                                      * (e(-om * mp.exp(x)) + e(om * mp.exp(x))), [-9, -4, -2, -1, 0] + [k / 8 for k in range(1, 57)])
    print("(3) Lemma 3.1.7 Delta=%s s0=%s beta=%s omega=%s: lhs=%s rhs=%s rel.diff=%.2e" % (Dl, mp.nstr(s0, 3), be, om, mp.nstr(lhs, 10), mp.nstr(rhs, 10), abs(lhs - rhs) / abs(lhs)))

# (4) Lemma 3.1.2: |D(0, 2sigma-1, 0, H/K)| <= 2K^2
worst = 0
for K in (range(1, 13) if '4' in SEC else []):
    for H in range(0 if K == 1 else 1, K):
        if gcd(H, K) != 1 and K > 1:
            continue
        for sg in (0.51, 0.6, 0.75):
            v = abs(D(0, 2 * sg - 1, 0, H, K)) / (2 * K ** 2)
            worst = max(worst, v)
if '4' in SEC: print("(4) Lemma 3.1.2: max |D(0,2sigma-1,0,H/K)|/(2K^2) over K<=12, sigma in {.51,.6,.75}: %.4f" % worst)
