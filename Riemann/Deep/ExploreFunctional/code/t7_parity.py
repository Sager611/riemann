# t7_parity.py (explore39): reflection identity Q(Eg) - Q(Rg) = -4 int int_{[0,a]^2} j(x+y) g(x) g(y); with g the positive odd
# ground state (odd sector Markov for c <= 2): lambda_ev <= lambda_odd - 2 I(c), I = int int J(x+y) g g / |g|^2_{[0,a]}.
import sys, numpy as np, math, json
sys.path.insert(0, '.scratch/explore39'); from xform import Form
rho = 1.324717957244746
for c in [1.2, rho, 1.5, rho**2, 1.9, 2.0]:
    F = Form(c); n = 20; a = F.a
    Qo, _ = F.Q(n, 'odd'); eo, Vo = np.linalg.eigh(Qo); v = Vo[:, 0]
    Qe, _ = F.Q(n, 'even'); ee = np.linalg.eigvalsh(Qe)
    xs, ws = F.gauss(0.0, a, 300); g = F.B(xs, n, 'odd') @ v
    if np.median(g) < 0: g = -g
    X, Y = np.meshgrid(xs, xs, indexing='ij'); Jxy = F.jc(X + Y)
    I = float(ws @ (Jxy*np.outer(g, g)) @ ws)/float(ws @ (g*g))
    print(json.dumps({'c': round(c, 6), 'lam_odd': round(float(eo[0]), 7), 'lam_even': round(float(ee[0]), 7), 'I(c)': round(I, 6),
                      'bound lam_odd-2I': round(float(eo[0]) - 2*I, 7), 'g>0 on (0,a]': bool(np.all(g > 0))}), flush=True)
