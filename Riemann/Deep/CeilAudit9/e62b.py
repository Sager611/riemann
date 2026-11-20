# ceilaudit9 E62b: kink budget of the hump-rearranged witness. Small humps (max|F| < b) take the sign of the preceding large hump;
# a kink (touch, no crossing) sits at every boundary where the witness keeps the sign but F changes it. Smoothing a kink of slope s over width delta
# costs (3 pi/8) s^2/delta of F'' energy and changes the joint law on ~2 delta per kink, so the witness is admissible in the limit iff
# Omega = (1/W) sum_kinks F'(gamma)^2 = o(<F''^2>) while #kept zeros = o(W).
import sys, numpy as np
sys.path.insert(0, 'Riemann/Deep/ExploreMalliavin')
from zdil import thp, zrs, vgrid
PI = np.pi
W = float(sys.argv[1])
for t0 in [float(x) for x in sys.argv[2].split(',')]:
    S = 32; v, t = vgrid(t0, W, 1/S); Z = zrs(t); F = Z/np.sqrt(2*thp(t)); F /= np.sqrt((F*F).mean())
    L = np.log(t0/(2*PI)); h = 1/S; Fp = np.gradient(F, h); Fpp = np.gradient(Fp, h); E2 = (Fpp**2).mean()
    zc = np.nonzero(np.sign(F[:-1])*np.sign(F[1:]) < 0)[0]; edges = np.concatenate([[0], zc+1, [len(F)]])
    H = np.array([np.abs(F[edges[k]:edges[k+1]]).max() for k in range(len(edges)-1)])
    sl = np.abs(Fp[zc])   # slope at each internal boundary (zero) between hump k and k+1
    for b in (L**0.25/np.log(L), 0.5, 0.3):
        big = H >= b; sgn = np.sign([F[edges[k]] if edges[k] < len(F) else 1 for k in range(len(H))])
        tau = np.empty(len(H)); cur = sgn[0]
        for k in range(len(H)):
            if big[k]: cur = sgn[k]
            tau[k] = cur
        keep = tau[:-1] != tau[1:]                  # witness crosses here (a true zero of F kept)
        kink = ~keep                                # F crosses, witness touches
        Om_all = (sl**2).sum()/W; Om = (sl[kink]**2).sum()/W
        print(f"t0={t0:.0e} L={L:.2f} b={b:.3f}: humps {len(H)}, large {big.mean():.3f}; witness zeros/unit {keep.sum()/W:.4f} (zeta {len(zc)/W:.4f}); "
              f"Omega_kinks/<F''^2> {Om/E2:.4f} (all {len(zc)} zeros rearranged: {Om_all/E2:.4f}); kinks {kink.mean():.3f} of zeros; "
              f"cost at delta=0.1: {(3*PI/8)*Om/0.1/E2:.3f} of <F''^2>", flush=True)
