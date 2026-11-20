# ceilaudit9 E62: untwisted inputs vs explore62's dilution witness (zeta data; Riemann-Siegel from ExploreMalliavin/zdil.py).
# (a) macroscopic K(alpha) = <F^2(v)F^2(v+x)>/<F^2>^2, alpha fixed in t (x = alpha*theta'/pi): witness (Gaussian or energy-window bulk,
#     independent placement) has K >= ~1 there; (b) share of one-point moments <|F|^m>, <|F'|^m> carried by the bulk {|F| < b}, b = L^(1/4)/log L
#     (Thm 3' bulk Phi(Y) has F' ~ 0, so it misses the bulk share of <|F'|^m>); (c) repair: zeta's own humps, signs reassigned in blocks of M humps,
#     made zero-free by Fhat = sigma*sqrt(F^2 + d^2 F'^2): zeros per unit, TV distance of the joint (F, F') law, F'' energy, bulk energy.
import sys, numpy as np
sys.path.insert(0, 'Riemann/Deep/ExploreMalliavin')
from zdil import th, thp, zrs, vgrid
PI = np.pi
def tv2(a1, b1, a2, b2, ex, ey):
    H1, _, _ = np.histogram2d(a1, b1, [ex, ey]); H2, _, _ = np.histogram2d(a2, b2, [ex, ey])
    return 0.5*np.abs(H1/H1.sum() - H2/H2.sum()).sum()
def run(t0, W, S=32, M=50):
    v, t = vgrid(t0, W, 1/S); Z = zrs(t); F = Z/np.sqrt(2*thp(t)); F /= np.sqrt((F*F).mean())
    L = np.log(t0/(2*PI)); h = 1/S; Fp = np.gradient(F, h); Fpp = np.gradient(Fp, h)
    zc = np.nonzero(np.sign(F[:-1])*np.sign(F[1:]) < 0)[0]
    print(f"t0={t0:.0e} L={L:.2f} points={len(F)} zeros={len(zc)} per unit {len(zc)/W:.4f} M4={(F**4).mean():.3f} <F'^2>={(Fp**2).mean():.3f} (pi^2/3={PI**2/3:.3f}) <F''^2>={(Fpp**2).mean():.2f} (pi^4/5={PI**4/5:.2f})", flush=True)
    # (a)
    F2 = F*F; al = np.arange(0.2, 12.01, 0.2); Ks = []
    for a in al:
        k = int(round(a*thp(np.array([t0]))[0]/PI*S)); Ks.append((F2[:-k]*F2[k:]).mean())
    Ks = np.array(Ks); i = np.argmin(Ks)
    print("  (a) K(alpha), alpha in t: " + " ".join(f"{a:.1f}:{K:.3f}" for a, K in zip(al[::5], Ks[::5])) + f"; min K = {Ks[i]:.3f} at alpha = {al[i]:.1f}; fraction of alpha in [0.2,12] with K<1: {(Ks<1).mean():.2f}", flush=True)
    # (b)
    b = L**0.25/np.log(L)
    for bb, nm in ((b, 'b=L^(1/4)/logL'), (0.5, 'b=0.5'), (0.2, 'b=0.2')):
        B = np.abs(F) < bb
        sh = lambda X, m: (np.abs(X[B])**m).sum()/(np.abs(X)**m).sum()
        print(f"  (b) {nm}={bb:.3f}: |bulk| {B.mean():.3f}; bulk share of <|F|^m> m=.25,.5,1,2: " + ",".join(f"{sh(F,m):.3f}" for m in (.25,.5,1,2))
              + "; of <|F'|^m>: " + ",".join(f"{sh(Fp,m):.3f}" for m in (.25,.5,1,2)), flush=True)
    # (c) humps between consecutive sign changes; sigma constant on blocks of M humps
    edges = np.concatenate([[0], zc+1, [len(F)]]); hump = np.repeat(np.arange(len(edges)-1), np.diff(edges))
    sig = np.where((hump//M) % 2 == 0, 1.0, -1.0)
    ex = np.linspace(-2, 2, 17); q = np.quantile(np.abs(Fp), 0.99); ey = np.linspace(-q, q, 17)
    tvself = tv2(F, Fp, -F, -Fp, ex, ey)   # empirical asymmetry of zeta's own joint law (noise floor)
    print(f"  (c) TV((F,F'),(-F,-F')) noise floor {tvself:.4f}; small humps (H<b) {np.mean([np.abs(F[edges[k]:edges[k+1]]).max() < b for k in range(len(edges)-1)]):.3f} of humps", flush=True)
    for d in (0.0, 0.05, 0.2):
        Fh = sig*np.sqrt(F*F + (d*Fp)**2); Fhp = np.gradient(Fh, h); Fhpp = np.gradient(Fhp, h)
        nzh = int((np.sign(Fh[:-1])*np.sign(Fh[1:]) < 0).sum())
        # joint law (F,F') vs (Fh,Fh'); zeta's law symmetrized is the target (law of (F,F') is ~symmetric)
        tv = tv2(F, Fp, Fh, Fhp, ex, ey); tvs = tv2(np.concatenate([F,-F]), np.concatenate([Fp,-Fp]), Fh, Fhp, ex, ey)
        mom = ",".join(f"{(np.abs(Fhp)**m).mean()/(np.abs(Fp)**m).mean():.3f}" for m in (.25,.5,1,2))
        B = np.abs(F) < b
        print(f"  (c) d={d:.2f}: zeros/unit {nzh/W:.4f} (zeta {len(zc)/W:.4f}); TV(F,F') vs zeta {tv:.4f} (vs symmetrized {tvs:.4f}); <|Fh'|^m>/<|F'|^m> m=.25,.5,1,2: {mom}; "
              f"<Fh''^2>/<F''^2> {(Fhpp**2).mean()/(Fpp**2).mean():.3f}, on bulk |F|<b: {(Fhpp[B]**2).sum()/(Fpp**2).sum():.3f} (zeta's own {(Fpp[B]**2).sum()/(Fpp**2).sum():.3f})", flush=True)
if __name__ == '__main__':
    W = float(sys.argv[1])
    for t0 in [float(x) for x in sys.argv[2].split(',')]: run(t0, W)
