# E: N6 R393 Polya-Wiman barrier: do x-derivatives heal an off-line zero pair of a genus-0 function?
# Model: G(x) = prod_k (1 + x/k^2) = sinh(pi sqrt x)/(pi sqrt x) (zeros X_k = -k^2, z-spacing 1 with x = -z^2);
# replace the zeros z = K, K+1 by the pair z = K + 1/2 +- i*dz (x-zeros X, conj X). Track the zero of F^(n) that starts at X.
# Heuristic (derivative = imaginary translation of the local factor by n*dx/pi, dx = local x-spacing 2K): healed once n >~ pi*dz.
import mpmath as mp, math, sys
mp.mp.dps = 40
K = int(sys.argv[1]) if len(sys.argv) > 1 else 200
def G(x):
    s = mp.sqrt(x)
    return mp.sinh(mp.pi*s)/(mp.pi*s) if s != 0 else mp.mpf(1)
XK, XK1 = -mp.mpf(K)**2, -mp.mpf(K+1)**2
def H(x):  # entire: G without the zeros at -K^2, -(K+1)^2
    return G(x)/((1 - x/XK)*(1 - x/XK1))
def Hder(x0, mmax, r, npts=256):
    # Cauchy integrals: H^(m)(x0) = m!/(2 pi i) \oint H(z)/(z-x0)^{m+1} dz on |z - x0| = r
    vals = [H(x0 + r*mp.expjpi(2*mp.mpf(j)/npts)) for j in range(npts)]
    out = []
    for m in range(mmax+1):
        s = mp.mpf(0)
        for j in range(npts):
            w = mp.expjpi(-2*mp.mpf(j)*m/npts)
            s += vals[j]*w
        out.append(mp.factorial(m)*s/(npts*r**m))
    return out
def Fn(x, n, dz):
    X = -(mp.mpf(K) + mp.mpf(1)/2 + 1j*dz)**2
    # P(x) = (1 - x/X)(1 - x/conj X)
    a = 1/X; b = 1/mp.conj(X)
    P = (1 - x*a)*(1 - x*b); P1 = -(a + b) + 2*x*a*b; P2 = 2*a*b
    r = mp.mpf(K)  # radius: < distance to nearest zeros of H (about 2K from X_K)
    h = Hder(x, n+1, r)
    f = P*h[n] + n*P1*h[n-1] + n*(n-1)/2*P2*h[n-2] if n >= 2 else (P*h[0] if n == 0 else P*h[1] + P1*h[0])
    g = P*h[n+1] + (n+1)*P1*h[n] + (n+1)*n/2*P2*h[n-1] if n >= 1 else P*h[1] + P1*h[0]
    return f, g
for dz in [mp.mpf('0.3'), mp.mpf('1.0'), mp.mpf('3.0')]:
    X = -(mp.mpf(K) + mp.mpf(1)/2 + 1j*dz)**2
    x = X
    line = []
    for n in range(0, 25):
        for it in range(30):
            f, g = Fn(x, n, dz)
            step = f/g
            x = x - step
            if abs(step) < mp.mpf(10)**(-15)*abs(x): break
        zim = mp.im(mp.sqrt(-x))  # back to z = sqrt(-x)
        line.append(f"{n}:{float(abs(mp.im(x))):.3g}/{float(abs(zim)):.3g}")
        if abs(mp.im(x)) < 1e-12*abs(x): break
    print(f"K={K} dz={float(dz)} (pred heal n ~ pi*dz = {math.pi*float(dz):.1f}): n:|Im x|/|Im z| ->", " ".join(line), flush=True)
