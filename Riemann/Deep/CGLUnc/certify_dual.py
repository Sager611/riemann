"""Arb certificate: an S-periodic configuration X (density nu = p/S) with energy e(X) < nu - H_W for a window W of the
16-harmonic family. Then no TEN-type certificate (local Gram-defect bound c*r - B*span <= Dcirc valid for all
configurations, combined with AF's chain) proves kappa > nu for W.
  H_W  = 3/2 - cot(t0)/sqrt2 - sum_j c_j^2 (1/2 - 1/(4 pi^2 j^2)) / (2 sin^2 t0),  t0 = 1/sqrt2
  Khat(a) = (psi*psi)(a) / (2 sin^2 t0),  psi = cos(sqrt2 s) + sum_j c_j cos(2 pi j s) on [-1/2, 1/2]
  e(X) = S^-2 [Khat(0) p^2 + 2 sum_{m=1}^{S-1} Khat(m/S) |sum_x exp(2 pi i m x/S)|^2] - p/S   (exact Poisson sum)
Usage: certify_dual.py dual.json [window: TEN | json-file with c]"""
import sys, json
from fractions import Fraction
from flint import arb, acb, ctx

ctx.prec = 200
CAMN = [14540341, -17051850, 6940151, 5482678, -200220, 3939293, -5615406, 3367926, -4987106,
        5017749, -4274153, 5029336, -4508275, 744709, 541564, -371255]


def window(arg):
    if arg == "TEN":
        return [Fraction(1)] + [Fraction(v, 10 ** 9) for v in CAMN]
    d = json.load(open(arg))
    return [Fraction(1)] + [Fraction(round(v * 10 ** 12), 10 ** 12) for v in d["c"][1:]]


def A(q):
    return arb(q.numerator) / arb(q.denominator)


def main():
    dual = json.load(open(sys.argv[1]))
    cw = window(sys.argv[2] if len(sys.argv) > 2 else "TEN")
    c = [A(q) for q in cw]
    pi = arb.pi()
    s2 = arb(2).sqrt()
    t0 = 1 / s2
    om = [s2] + [2 * pi * j for j in range(1, 17)]
    sin2 = t0.sin() ** 2
    H = arb(3) / 2 - (t0.cos() / t0.sin()) / s2
    for j in range(1, 17):
        H -= c[j] ** 2 * (arb(1) / 2 - 1 / (4 * pi ** 2 * j ** 2)) / (2 * sin2)

    def I(a, b, Om, ph):  # int_a^b cos(Om s + ph) ds
        if Om == 0:
            return (b - a) * ph.cos()
        return ((Om * b + ph).sin() - (Om * a + ph).sin()) / Om

    def conv(al):  # (psi*psi)(al), 0 <= al < 1
        lo, hi = arb(-1) / 2, arb(1) / 2 - al
        tot = arb(0)
        for j in range(17):
            for k in range(17):
                if c[j] == 0 or c[k] == 0:
                    continue
                d = 0 if j == k else om[j] - om[k]
                v = I(lo, hi, d, -om[k] * al) + I(lo, hi, om[j] + om[k], om[k] * al)
                tot += c[j] * c[k] * v / 2
        return tot

    S = int(dual["S"]); xs = [Fraction(round(v * 10 ** 9), 10 ** 9) for v in dual["x"]]
    p = len(xs)
    kh0 = conv(arb(0)) / (2 * sin2)
    tot = kh0 * p * p
    for m in range(1, S):
        kh = conv(arb(m) / S) / (2 * sin2)
        mu = acb(0)
        for x in xs:
            mu += (2 * pi * acb(0, 1) * m * A(x) / S).exp()
        tot += 2 * kh * abs(mu) ** 2
    e = tot / S ** 2 - arb(p) / S
    nu = arb(p) / S
    slack = nu - H - e
    print("S =", S, " p =", p, " nu =", Fraction(p, S), "=", nu.str(12))
    print("H_W   =", H.str(25))
    print("e(X)  =", e.str(25))
    print("slack = nu - H_W - e(X) =", slack.str(15), " certified > 0:", slack > 0)


main()
