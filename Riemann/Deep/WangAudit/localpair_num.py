# Rough numerical sanity check of thm:local-pair / (eq:unweighted) normalisation with actual zeros (delta = 0).
# Expect only rough agreement at T ~ 10^3 (errors O(H) vs main HL, i.e. O(1/L) ~ 15-25%).
import sys, numpy as np
G = np.array([float(l.split()[1]) for l in open(sys.argv[1])])
def Cl(l): return l/2 + 1/(np.sqrt(2)*np.tan(l/np.sqrt(2)))
def K(xi, l):  # Fourier transform of f_l(u) = cos(sqrt2 u)/(sqrt2 sin(l/sqrt2)) on [-l/2, l/2], real xi
    a = l/2; b = np.sqrt(2); c = 2*np.pi*xi
    def sinc_part(w):
        w = np.where(np.abs(w) < 1e-12, 1e-12, w)
        return np.sin(w*a)/w
    return (sinc_part(b-c) + sinc_part(b+c))/(np.sqrt(2)*np.sin(l/np.sqrt(2)))
for (T, th) in ((1000.0, 0.9), (1200.0, 0.75), (600.0, 0.95)):
    H = T**th
    if T + H > G[-1]: print('skip', T, th); continue
    g = G[(G > T) & (G <= T+H)]; N = len(g)
    d = g[:, None] - g[None, :]
    w = 4/(4 + d**2)
    for Lname, L in (('logT', np.log(T)), ('log(T/2pi)', np.log((T+H/2)/(2*np.pi)))):
        print(f'T={T:.0f} theta={th} H={H:.1f} N(I)={N} HL/2pi={H*L/(2*np.pi):.1f} [{Lname}]')
        for lf in (0.5, 0.8, 0.95):
            l = lf*th
            k2 = K(d*L/(2*np.pi), l)**2
            SK = k2.sum(); W = (k2*w).sum()
            print(f'   lambda={l:.3f}: S_K/N={SK/N:.4f}  C_lambda={Cl(l):.4f}  2-S_K/N={2-SK/N:.4f}  W_I(Q)/(HL/2pi)={W/(H*L/(2*np.pi)):.4f}')
        al = np.array([0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0, 1.2, 1.5])
        Fv = [ (np.cos(d*L*a)*w).sum()/(H*L/(2*np.pi)) for a in al ]
        pred = al + L*np.exp(-2*L*al)
        print('   F(alpha)=2pi F_I(T^a)/(HL):', ' '.join(f'{a:.2f}:{f:.3f}/{p:.3f}' for a, f, p in zip(al, Fv, pred)))
