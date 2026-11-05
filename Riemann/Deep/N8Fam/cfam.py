# C_fam(f) for dilated windows under several family form-factor models; plain and certified bounds.
# Run: OMP_NUM_THREADS=1 uv run --no-project --with numpy python3 cfam.py
# Window f: even probability density on I=[-1/2,1/2]; f_lam(u) = f(u/lam)/lam on [-lam/2, lam/2].
# Energy constant C(f_lam) = int f_lam^2 + iint phi(|s-t|) f_lam(s) f_lam(t)
#                          = (1/lam) int f^2 + int_{-1}^{1} phi(lam |d|) R(d) dd,  R = f*f (autocorrelation).
# phi = form factor without the delta: zeta (BGST) phi(e)=e on [0,1].
import numpy as np
A2 = np.array([122290, -45012, 16241, 1887, 2711, 3731, -5307, 2982, -5153, 5358]) * 1e-6  # paper's second window
X_CERT = 8722e-6          # x = 4361/500000
B_CERT = 7 * 627e-6       # b = 7 alpha = 4389e-6
C0P = 6639737 / 5000000   # paper's rational cap for C(f)
BQ = 0.87338196661320752  # zero-free half-plane exponent (Bebel / Liu B_2), TREE.md node N7
GL_X, GL_W = np.polynomial.legendre.leggauss(160)

def gl(a, b):
    x = 0.5 * (b - a) * GL_X + 0.5 * (b + a)
    return x, 0.5 * (b - a) * GL_W

def cos_window(a):
    a = np.asarray(a, float)
    def f(u):
        u = np.asarray(u, float)
        return 1 + sum(aj * np.cos(2 * np.pi * (j + 1) * u) for j, aj in enumerate(a))
    return f

def mt_window():
    k = np.sqrt(2.0)
    norm = 2 * np.sin(k / 2) / k
    return lambda u: np.cos(k * np.asarray(u, float)) / norm

def flat():
    return lambda u: np.ones_like(np.asarray(u, float))

def R_of(f, d):
    # autocorrelation at d in [0,1]
    out = np.empty_like(d)
    for i, di in enumerate(d):
        s, w = gl(-0.5, 0.5 - di)
        out[i] = np.sum(w * f(s) * f(s + di))
    return out

def phi_factory(model, c=1.263, B=BQ):
    if model == 'zeta':
        return (lambda e: e), []
    if model == 'gue':
        return (lambda e: np.minimum(e, 1.0)), [1.0]
    if model == 'A':      # exact |alpha| up to 1, large-sieve bound c|alpha| on (1,2)
        return (lambda e: np.where(e <= 1, e, c * e)), [1.0]
    if model == 'AQ':     # exact up to 1, GUE value 1 on (1,1/B) (quasi-GRH), c|alpha| beyond
        return (lambda e: np.where(e <= 1, e, np.where(e < 1 / B, 1.0, c * e))), [1.0, 1 / B]
    if model == 'Q':      # exact up to 1, GUE value 1 on (1,1/B); window must have lam <= 1/B
        return (lambda e: np.where(e <= 1, e, 1.0)), [1.0]
    raise ValueError(model)

def C_model(f, lam, model, **kw):
    phi, kinks = phi_factory(model, **kw)
    s, w = gl(-0.5, 0.5)
    I0 = np.sum(w * f(s) ** 2)
    cuts = sorted(set([0.0, 1.0] + [k / lam for k in kinks if 0 < k / lam < 1]))
    I1 = 0.0
    for a, b in zip(cuts[:-1], cuts[1:]):
        d, wd = gl(a, b)
        I1 += 2 * np.sum(wd * phi(lam * d) * R_of(f, d))
    return I0 / lam + I1

def report(name, f, model, lam_max, **kw):
    lams = np.linspace(0.5, lam_max, 301)
    C = np.array([C_model(f, l, model, **kw) for l in lams])
    plain = 2 - C
    cert = (2 - C - lams * B_CERT) / (1 - X_CERT)
    i, j = int(np.argmax(plain)), int(np.argmax(cert))
    print('%-6s %-4s %-18s plain max %.6f at lam %.4f | scaled cert max %.6f at lam %.4f (plain there %.6f)'
          % (name, model, kw, plain[i], lams[i], cert[j], lams[j], plain[j]))
    return lams, plain, cert

if __name__ == '__main__':
    f2 = cos_window(A2)
    print('check zeta C, paper window: %.12f (paper 1.32794727284)' % C_model(f2, 1.0, 'zeta'))
    print('check zeta C, flat: %.12f (4/3)' % C_model(flat(), 1.0, 'zeta'))
    print('check zeta C, MT: %.12f (1/2+cot(1/sqrt2)/sqrt2 = %.12f)' % (C_model(mt_window(), 1.0, 'zeta'), 0.5 + 1 / np.sqrt(2) / np.tan(1 / np.sqrt(2))))
    print('check gue lam=2 flat: %.12f (13/12)' % C_model(flat(), 2.0, 'gue'))
    print('paper cert at lam=1: s = %.7f (paper 0.6735382 with cap C0P: %.7f)' % ((2 - C_model(f2, 1.0, 'zeta') - B_CERT) / (1 - X_CERT), (2 - C0P - B_CERT) / (1 - X_CERT)))
    for name, f in [('paper', f2), ('flat', flat()), ('MT', mt_window())]:
        report(name, f, 'gue', 1.9999)
        report(name, f, 'A', 1.9999)
        report(name, f, 'A', 1.9999, c=1.0)
        report(name, f, 'AQ', 1.9999)
        report(name, f, 'Q', 1 / BQ)
        report(name, f, 'zeta', 1.0)
    for name, f in [('flat', flat()), ('MT', mt_window()), ('paper', f2)]:
        print('%s model A at lam=1.5: plain %.6f' % (name, 2 - C_model(f, 1.5, 'A')))
