# N1: Lemma Z mechanism. I_j(a) = ∫ f(Σv) ∏ (ρ_η(v_i)/v_i) e^{a v_i} dv (j copies of one planted zero of depth a/L at t = γ)
# versus the lossy route |I_1|^j (DPC's |e_r| << (1+|S_1|)^r). Units: v = u/L, X^d = e^{θ a}.
import numpy as np
th, eta, b = 0.55, 0.1, 0.62
def bump01(x):  # smooth step 0 (x<=0) -> 1 (x>=1)
    x = np.clip(x, 0, 1); g = lambda y: np.where(y > 0, np.exp(-1/np.maximum(y, 1e-300)), 0.0)
    return g(x)/(g(x)+g(1-x))
def rho(v): return (1-bump01((v-th)/(b-th)))*bump01(v/eta-1)
def f(x): return np.where(x < th, np.exp(1-1/np.maximum(1-(x/th)**2, 1e-300)), 0.0)
def I(j, a, n):
    v = np.linspace(eta, b, n); dv = v[1]-v[0]; g = rho(v)/v*np.exp(a*v)
    if j == 1: return np.sum(f(v)*g)*dv
    if j == 2:
        S = v[:,None]+v[None,:]; return np.sum(f(S)*g[:,None]*g[None,:])*dv*dv
    S = v[:,None,None]+v[None,:,None]+v[None,None,:]
    return np.sum(f(S)*g[:,None,None]*g[None,:,None]*g[None,None,:])*dv**3
print("a   e^{θa}     I1        I2        I3      I1^2      I1^3   I2/e^{θa} I3/e^{θa}")
for a in [0, 2, 5, 10, 20, 40]:
    i1, i2, i3 = I(1,a,40001), I(2,a,1601), I(3,a,161)
    E = np.exp(th*a)
    print(f"{a:3d} {E:9.3g} {i1:9.3g} {i2:9.3g} {i3:9.3g} {i1**2:9.3g} {i1**3:9.3g} {i2/E:9.3g} {i3/E:9.3g}")

