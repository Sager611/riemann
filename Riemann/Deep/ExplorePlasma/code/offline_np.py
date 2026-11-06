# explore20: xi-like kernel with the first two zeros merged into an off-line quadruple gamma_m +- i delta (numpy, float64)
# Xi_off(z) = Xi(z) N(z^2) / ((z^2 - g1^2)(z^2 - g2^2)),  N(T) = (T - gm^2 + d^2)^2 + 4 gm^2 d^2
#          = Xi + A Xi/(z^2-g1^2) + B Xi/(z^2-g2^2),  A = N(g1^2)/(g1^2-g2^2), B = N(g2^2)/(g2^2-g1^2)
# kernel (Xi = 2 int_R Phi e^{izu}):  Phi_off = Phi + A psi1 + B psi2,
#   psi_j(u) = -(1/g) int_u^inf sin(g(w-u)) Phi(w) dw = -(1/g)[cos(gu) Is(u) - sin(gu) Ic(u)],  Is/Ic = int_u^inf sin/cos(gw) Phi(w) dw
#   psi_j'(u) = cos(gu) Ic(u) + sin(gu) Is(u),  psi_j'' = -Phi - g^2 psi_j   (psi even, psi'(0) = Xi(g)/4 = 0)
import numpy as np
G1, G2 = 14.134725141734693790, 21.022039638771554993
def Phi_np(w, k=0, N=10):
    w = np.asarray(w, dtype=float); x = np.exp(2*w); s = np.zeros_like(w)
    for n in range(1, N+1):
        a = np.pi*n*n
        for c, al in ((2*a*a, 4.5), (-3*a, 2.5)):
            e = c*np.exp(al*w - a*x)
            s += e if k == 0 else ((al - 2*a*x)*e if k == 1 else ((al - 2*a*x)**2 - 4*a*x)*e)
    return s
def grid(U=4.0, h=5e-4, m=8):
    xg, wg = np.polynomial.legendre.leggauss(m)
    edges = np.arange(0.0, U + h/2, h)
    lo, hi = edges[:-1], edges[1:]
    nodes = (lo[:, None] + hi[:, None])/2 + (hi - lo)[:, None]/2*xg[None, :]
    weights = (hi - lo)[:, None]/2*wg[None, :]
    return edges, nodes, weights
def psis(g, edges, nodes, weights):
    P = Phi_np(nodes)
    Is_int = (np.sin(g*nodes)*P*weights).sum(1); Ic_int = (np.cos(g*nodes)*P*weights).sum(1)
    Is = np.concatenate([np.cumsum(Is_int[::-1])[::-1], [0.0]]); Ic = np.concatenate([np.cumsum(Ic_int[::-1])[::-1], [0.0]])
    u = edges
    psi = -(np.cos(g*u)*Is - np.sin(g*u)*Ic)/g
    dpsi = np.cos(g*u)*Ic + np.sin(g*u)*Is
    d2psi = -Phi_np(u) - g*g*psi
    return psi, dpsi, d2psi
def coeffs(delta, g1=G1, g2=G2):
    gm = (g1 + g2)/2; d = float(delta)
    N = lambda T: (T - gm**2 + d**2)**2 + 4*gm**2*d**2
    return gm, d, N(g1**2)/(g1**2 - g2**2), N(g2**2)/(g2**2 - g1**2), N
def make(delta, U=4.0, h=5e-4):
    edges, nodes, weights = grid(U, h)
    gm, d, A, B, N = coeffs(delta)
    p1, dp1, d2p1 = psis(G1, edges, nodes, weights); p2, dp2, d2p2 = psis(G2, edges, nodes, weights)
    F0 = Phi_np(edges) + A*p1 + B*p2
    F1 = Phi_np(edges, 1) + A*dp1 + B*dp2
    F2 = Phi_np(edges, 2) + A*d2p1 + B*d2p2
    return edges, F0, F1, F2, (gm, d, A, B, N, dp1[0], dp2[0])
