# N3: numerical check of the multivariate explicit formula behind Lemma Z (j = 1, 2), real zeta zeros (mpmath),
# t near 1000. LHS_j = sum over n_1..n_j of prod Λ(n_i)φ(log n_i) f(Σ log n_i/L) (n_1⋯n_j)^{-1/2-it};
# RHS_j = (-1)^j Σ_{ρ_1..ρ_j} ∫ f(Σu/L) ∏ φ(u_i) e^{(ρ_i-1/2-it)u_i} du (zeros with |γ - t| ≤ Γ).
import numpy as np, mpmath as mp, sys
L, th, eta, b = np.log(1000.0), 0.55, 0.15, 0.75
def bump01(x):
    x = np.clip(x, 0, 1); g = lambda y: np.where(y > 0, np.exp(-1/np.maximum(y, 1e-300)), 0.0)
    return g(x)/(g(x)+g(1-x))
rho = lambda v: (1-bump01((v-th)/(b-th)))*bump01(v/eta-1)
phi = lambda u: rho(u/L)/u
f = lambda x: np.where(x < th, np.exp(1-1/np.maximum(1-(x/th)**2, 1e-300)), 0.0)
# von Mangoldt on n < e^{bL}
N = int(np.exp(b*L))+2
Lam = np.zeros(N)
for p in range(2, N):
    if all(p % q for q in range(2, int(p**0.5)+1)):
        pk = p
        while pk < N: Lam[pk] = np.log(p); pk *= p
ns = np.arange(N); w = np.where(Lam > 0, Lam*phi(np.log(np.maximum(ns,2))), 0.0); idx = np.nonzero(w)[0]
zs = [float(mp.im(mp.zetazero(n))) for n in range(int(sys.argv[1]), int(sys.argv[2]))]
gam = np.array(zs)
u = np.linspace(eta*L*0.999, b*L*1.001, 3001); du = u[1]-u[0]
for t in (1000.0, 1003.7):
    Z = np.exp(1j*np.outer(u, gam - t)).sum(axis=1)                    # Σ_ρ e^{(ρ-1/2-it)u}, RH numerically
    lhs1 = sum(w[n]*f(np.log(n)/L)*n**(-0.5-1j*t) for n in idx)
    rhs1 = -np.sum(f(u/L)*phi(u)*Z)*du
    lhs2 = sum(w[m]*w[n]*f(np.log(m*n)/L)*(m*n)**(-0.5-1j*t) for m in idx for n in idx if m*n < np.exp(th*L))
    g = phi(u)*Z; S = u[:,None]+u[None,:]
    rhs2 = np.sum(f(S/L)*g[:,None]*g[None,:])*du*du
    print(f"t={t}: zeros={len(gam)} j=1 LHS={lhs1:.5f} RHS={rhs1:.5f} |diff|={abs(lhs1-rhs1):.2e}; j=2 LHS={lhs2:.5f} RHS={rhs2:.5f} |diff|={abs(lhs2-rhs2):.2e}", flush=True)

