# Check of Formalism.HN.SlopeIneq (HarderNarasimhan/Statement.lean:91-100; proved: slope_inequality_general,
# HarderNarasimhanProofs/SlopeIneq.lean:507) on two explicit lattice maps. Conventions (LatticeSlopes/Statement.lean:15-25):
# covol = sqrt(det Gram) of the Z-lattice, deg = -log covol, slope = deg / Z-rank; muMax over nonzero R-submodules,
# muMin over quotients by saturated proper R-submodules. R = Z (Ex. 1) and R = Z[i] (Ex. 2: hermitian, real metric Re h).
import itertools, numpy as np
from math import log, sqrt

def zlat_stats(G, box=6):
    """R = Z, rank 2 or 3: return deg, rk, muMax, muMin by enumeration (saturated sublattices via primitive vectors)."""
    n = G.shape[0]; detG = np.linalg.det(G); deg = -0.5 * log(detG)
    vecs = [np.array(v) for v in itertools.product(range(-box, box + 1), repeat=n) if any(v)]
    norm = lambda v: sqrt(v @ G @ v)
    lam1 = min(norm(v) for v in vecs)
    slopes_sub = [deg / n, -log(lam1)]           # rank-n and rank-1 subs
    quot = [deg / n, deg + log(lam1) if n == 2 else (deg + log(lam1)) / (n - 1)]  # quotient by Zv (rank n-1)
    if n == 3:
        Gi = np.linalg.inv(G); detv = lambda v: sqrt(v @ Gi @ v)   # dual norm
        lamd = min(detv(v) for v in vecs)
        # saturated rank-2 N = M cap v^perp (v primitive dual): covol(N) = covol(M) * ||v||_dual
        slopes_sub.append((deg - log(lamd)) / 2)
        quot.append(deg - (deg - log(lamd)))     # quotient by rank-2 N: rank 1, deg L - deg N = -log ||v||_dual... = log(1/lamd)
    return deg, n, max(slopes_sub), min(quot)

def gauss_stats(H, box=3):
    """R = Z[i], rank 2 over Z[i] (Z-rank 4). Real metric Re h. Saturated Z[i]-submodules of Z[i]-rank 1: Z[i]v."""
    detH = np.linalg.det(H).real; deg = -log(detH)       # covol of the Z-lattice = det H (det of real Gram = det(H)^2)
    rng = range(-box, box + 1)
    best = None
    for a, b, c, d in itertools.product(rng, rng, rng, rng):
        if a == b == c == d == 0: continue
        v = np.array([a + 1j * b, c + 1j * d]); nv = sqrt((v.conj() @ H @ v).real)
        best = nv if best is None or nv < best else best
    lam1 = best
    # deg(Z[i]v) = -log(covol) = -log(||v||^2) (Z[i] has covolume 1 in C); slope over Z-rank 2 = -log ||v||
    muMax = max(deg / 4, -log(lam1)); muMin = min(deg / 4, (deg + 2 * log(lam1)) / 2)
    return deg, 4, muMax, muMin

def opnorm(P, GL, GM):  # sup ||P x||_M / ||x||_L
    w = np.linalg.eigvals(np.linalg.solve(GL, P.conj().T @ GM @ P))
    return sqrt(max(w.real))

print("Ex. 1 (R = Z): L = Z^2, Gram A = [[2,1],[1,3]]; M = Z^3, Gram B = [[1,0,0],[0,2,1],[0,1,5]]; phi = P (3x2 integer, injective)")
A = np.array([[2., 1.], [1., 3.]]); B = np.array([[1., 0, 0], [0, 2., 1.], [0, 1., 5.]])
for P in (np.array([[1., 0], [1., 1.], [0, 1.]]), np.array([[3., 1.], [0, 2.], [1., -1.]])):
    C = opnorm(P, A, B)
    dL, rL, MxL, MnL = zlat_stats(A); dM, rM, MxM, MnM = zlat_stats(B)
    dI, rI, MxI, MnI = zlat_stats(P.T @ B @ P)
    print(f" P={P.tolist()}: C={C:.4f} logC={log(C):.4f} | L: deg={dL:.4f} muMax={MxL:.4f} muMin={MnL:.4f} | M: muMax={MxM:.4f} | im: muMin={MnI:.4f}")
    print(f"   (b) muMax L <= muMax M + log C : {MxL:.4f} <= {MxM + log(C):.4f}  {MxL <= MxM + log(C) + 1e-12}")
    print(f"   (a) muMin L <= muMax M + log C : {MnL:.4f} <= {MxM + log(C):.4f}  {MnL <= MxM + log(C) + 1e-12};  muMin L - log C <= muMin im : {MnL - log(C):.4f} <= {MnI:.4f}  {MnL - log(C) <= MnI + 1e-12}")
    # divisibility: 2P lands in 2M; as a map into M the archimedean bound doubles, as a map into the sublattice 2M it is unchanged
    C2 = opnorm(2 * P, A, B); _, _, Mx2M, _ = zlat_stats(4 * B)
    print(f"   divisible map 2P: into M needs log C' = {log(C2):.4f}; into 2M (slope muMax(2M) = {Mx2M:.4f} = muMax M - log 2) same bound {Mx2M + log(C2):.4f}")

print("Ex. 2 (R = Z[i], hermitian): L = Z[i]^2, H_L = [[2,1+i],[1-i,3]]; M = Z[i]^2, H_M = [[1,0],[0,2]]")
HL = np.array([[2, 1 + 1j], [1 - 1j, 3]]); HM = np.array([[1, 0], [0, 2]], dtype=complex)
for P in (np.array([[1 + 1j, 1], [0, 2 - 1j]]), np.array([[1, 1j], [1, -1j]])):
    C = opnorm(P, HL, HM)
    dL, rL, MxL, MnL = gauss_stats(HL); dM, rM, MxM, MnM = gauss_stats(HM); dI, rI, MxI, MnI = gauss_stats(P.conj().T @ HM @ P)
    print(f" P={P.tolist()}: C={C:.4f} logC={log(C):.4f} | L: deg={dL:.4f} muMax={MxL:.4f} muMin={MnL:.4f} | M: muMax={MxM:.4f} | im: muMin={MnI:.4f}")
    print(f"   (b) {MxL:.4f} <= {MxM + log(C):.4f} {MxL <= MxM + log(C) + 1e-12}; (a) {MnL:.4f} <= {MxM + log(C):.4f} {MnL <= MxM + log(C) + 1e-12}; image: {MnL - log(C):.4f} <= {MnI:.4f} {MnL - log(C) <= MnI + 1e-12}")
# Sharpness: for an isometric embedding of a semistable L the inequality (b) is an equality
I2 = np.eye(2); P = np.array([[1., 0], [0, 1.], [0, 0]]); B3 = np.diag([1., 1., 1.])
print("Sharpness: L = Z^2 standard into Z^3 standard isometrically: muMax L =", zlat_stats(I2)[2], " muMax M + log C =", zlat_stats(B3)[2] + log(opnorm(P, I2, B3)))
