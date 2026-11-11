# T3 (flutter side): Theodorsen's lift-deficiency function is passive (positive real); no zeros of its denominator on the
# physical sheet (argument principle); Re C(k) in (1/2, 1] <=> |H0(k)| < |H1(k)| (Nicholson). Wagner: 1 - phi(s) completely monotone.
import mpmath as mp
mp.mp.dps = 15
import functools, sys
print = functools.partial(print, flush=True)
H = lambda n, k: mp.hankel2(n, k)
C = lambda k: H(1, k)/(H(1, k) + 1j*H(0, k))
Cl = lambda p: mp.besselk(1, p)/(mp.besselk(0, p) + mp.besselk(1, p))    # Laplace form, p = i k
mx = (0, None); okRe = True
for j in range(1, 81):
    k = mp.mpf(j)/4
    c = C(k); r = abs(H(0, k))/abs(H(1, k))
    if not (0.5 < c.real <= 1 and c.imag <= 0 and r < 1): okRe = False
    d = abs(c - Cl(1j*k))
    if d > mx[0]: mx = (d, k)
print(f"k in (0, 20]: 1/2 < Re C <= 1, Im C <= 0, |H0|<|H1| all hold: {okRe}; max |C(k) - K1/(K0+K1)(ik)| = {mp.nstr(mx[0], 3)}")
print("C(k) samples:", [(float(k), mp.nstr(C(k), 6)) for k in [0.05, 0.2, 1, 5, 20]])
# argument principle for D(p) = K0(p) + K1(p) on the physical sheet |arg p| < pi: keyhole contour, r = 1e-6 .. 60
def winding(R=30, r=1e-4, eps=1e-9, n=300):
    D = lambda p: mp.besselk(0, p) + mp.besselk(1, p)
    pts = []
    for j in range(n + 1): pts.append(R*mp.expj(-mp.pi + eps + (2*mp.pi - 2*eps)*j/n))          # big circle ccw
    for j in range(n + 1): pts.append(mp.expj(mp.pi - eps)*(R - (R - r)*j/n))                    # upper side of cut inwards
    for j in range(n + 1): pts.append(r*mp.expj(mp.pi - eps - (2*mp.pi - 2*eps)*j/n))           # small circle cw
    for j in range(n + 1): pts.append(mp.expj(-mp.pi + eps)*(r + (R - r)*j/n))                  # lower side outwards
    tot = 0
    vals = [D(p) for p in pts]
    for a, b in zip(vals[:-1], vals[1:]): tot += mp.im(mp.log(b/a))
    return tot/(2*mp.pi)
print(f"zeros of K0 + K1 inside the keyhole (physical sheet, 1e-4 < |p| < 30): {mp.nstr(winding(), 6)}")
