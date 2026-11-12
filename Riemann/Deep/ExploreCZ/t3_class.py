
import numpy as np, sys
sys.path.insert(0,'.')
# exponent bookkeeping: R <= N^{LV}; T = N^tau, V = N^sigma, |a_n|<=1, G_k <= N^{k+o(1)}
def LV_class(s,tau,K=12):   # positivity class = k-th power MVT (local L^{2k} masses only)
    return min(max(tau,k)+k*(1-2*s) for k in range(1,K+1))
def LV_hux(s,tau): return max(2-2*s, tau+4-6*s)
def LV_gm(s,tau):  return max(2-2*s, 18/5-4*s, tau+12/5-4*s)
taus=[1.2,1.5,2.0,2.5,3.0]; sigs=[0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]
print("LV exponents: cell = class | best(Hux,GM,class)  [* = phase needed: best < class]")
print("sigma \\ tau "+"".join(f"{t:>16}" for t in taus))
for s in sigs:
    row=f"{s:10.2f}  "
    for t in taus:
        c=LV_class(s,t); b=min(c,LV_hux(s,t),LV_gm(s,t))
        row+=f"{c:7.3f}|{b:6.3f}{'*' if b<c-1e-9 else ' '} "
    print(row)
# density: Guth-Maynard sec.13 / Ingham detector reduction (copied from ExploreKakeya t1_energy_lp.py)
def density(s, LV, nus=np.linspace(0.01,0.5,491)):
    worst=2*(1-s)
    for nu in nus:
        vals=[]; k=1
        while k*nu<=3.0:
            x=k*nu; tau=1/x
            lv = 2-2*s if tau<=1 else LV(s,tau)
            vals.append(x*lv); k+=1
        worst=max(worst,min(vals))
    return worst/(1-s)
print()
print("density A(sigma): class (positivity) | Ingham 3/(2-s) | class+Huxley | Huxley 3/(3s-1) | GM 15/(3+5s) on [.7,.8]")
for s in [0.55,0.6,0.65,0.7,0.72,0.75,0.78,0.8,0.85,0.9]:
    a=density(s,LV_class); ah=density(s,lambda s_,t_: min(LV_class(s_,t_),LV_hux(s_,t_)))
    gm=15/(3+5*s) if 0.7<=s<=0.8 else float('nan')
    print(f"  {s:.2f}: {a:.4f} | {3/(2-s):.4f} | {ah:.4f} | {3/(3*s-1):.4f} | {gm:.4f}")
# profile admissibility (Fejer spikes) -- explicit constants
lam=np.log(2.0)
u=np.linspace(-4000,4000,1600001); du=u[1]-u[0]
F=np.where(np.abs(u)<1e-12,1.0,(np.sin(lam*u/2)/(lam*u/2+1e-300))**2)
IF={k:float(np.sum(F**k)*du) for k in (1,2,3,4)}
print()
print("Fejer F(t)=(sin(t log2/2)/(t log2/2))^2, Fhat supp [-log2,log2]: int F^k =",{k:round(v,4) for k,v in IF.items()}, "(int F = 2pi/log2 =",round(2*np.pi/lam,4),")")
def overlap(L,k):  # max_t (sum_r F(t-rL))^k / F-mass factor: bound via sup of the comb
    t=np.linspace(0,L,2001); s=np.zeros_like(t)
    for r in range(-2000,2001): 
        x=t-r*L; s+=np.where(np.abs(x)<1e-12,1.0,(np.sin(lam*x/2)/(lam*x/2+1e-300))**2)
    # mass per period of comb^k divided by int F^k
    return float(np.trapezoid(s**k,t)/IF[k])
for L in (20.0,50.0,200.0):
    print(f"  comb spacing L={L}: (period mass of comb^k)/int F^k = "+", ".join(f"k{k}:{overlap(L,k):.4f}" for k in (1,2,3)))
