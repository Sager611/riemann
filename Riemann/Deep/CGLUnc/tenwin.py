"""TEN-class window data: psi, H_W = 2 - R(psi), kernel k_W, Khat = (psi*psi)/(int psi)^2 on [-1,1].
Window family (arena TEN): psi(s) = cos(sqrt2 s) + sum_{j=1}^{16} c_j cos(2 pi j s) on [-1/2,1/2]."""
import numpy as np

CAMN = [14540341, -17051850, 6940151, 5482678, -200220, 3939293, -5615406, 3367926, -4987106,
        5017749, -4274153, 5029336, -4508275, 744709, 541564, -371255]
TEN_C = np.array([1.0] + [c / 1e9 for c in CAMN])
MT_C = np.array([1.0] + [0.0] * 16)
OM = np.array([np.sqrt(2.0)] + [2 * np.pi * j for j in range(1, 17)])

_xg, _wg = np.polynomial.legendre.leggauss(400)


def psi(s, c, om=OM):
    s = np.asarray(s, dtype=float)
    return np.cos(np.multiply.outer(s, om)) @ c


def quad(f, a, b):
    x = 0.5 * (b - a) * _xg + 0.5 * (a + b)
    return 0.5 * (b - a) * np.sum(_wg * f(x))


def int_psi(c, om=OM):
    return quad(lambda s: psi(s, c, om), -0.5, 0.5)


def R_of(c, om=OM):
    I1 = int_psi(c, om)
    I2 = quad(lambda s: psi(s, c, om) ** 2, -0.5, 0.5)
    # double integral of |u - v| psi(u) psi(v): 2 * int_0^1 w (psi*psi)(w) dw
    w = 0.5 * _xg + 0.5
    ww = 0.5 * _wg
    ac = np.array([autoc(t, c, om) for t in w])
    I3 = 2 * np.sum(ww * w * ac)
    return (I2 + I3) / I1 ** 2


def autoc(t, c, om=OM):
    t = abs(t)
    if t >= 1:
        return 0.0
    return quad(lambda s: psi(s, c, om) * psi(s + t, c, om), -0.5, 0.5 - t)


def H_of(c, om=OM):
    return 2 - R_of(c, om)


def Khat(alpha, c, om=OM):
    I1 = int_psi(c, om)
    return np.array([autoc(a, c, om) for a in np.atleast_1d(alpha)]) / I1 ** 2


def k_of(x, c, om=OM):
    I1 = int_psi(c, om)
    x = np.atleast_1d(np.asarray(x, dtype=float))
    return np.array([quad(lambda s: psi(s, c, om) * np.cos(2 * np.pi * xx * s), -0.5, 0.5) for xx in x]) / I1


if __name__ == "__main__":
    for name, c in [("MT", MT_C), ("TEN", TEN_C)]:
        print(name, "H =", repr(H_of(c)), " int psi =", int_psi(c), " Khat(0) =", Khat(0.0, c)[0],
              " min psi =", psi(np.linspace(-0.5, 0.5, 2001), c).min())
    print("H_MT exact 3/2 - cot(1/sqrt2)/sqrt2 =", 1.5 - 1 / np.tan(1 / np.sqrt(2)) / np.sqrt(2))
