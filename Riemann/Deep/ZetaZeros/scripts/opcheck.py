
# opcheck.py -- numerical check of the operator layer (AF/Lam construction as used in the paper).
# zeros on the line -> real z; off-line pairs {rho, 1-conj(rho)} -> {z, conj z};
# A = sum_z nu_z F_z <F_{conj z}, .>,  F_z(u) = sqrt(f(u)) e^{-2 pi i u z}.
# Nonzero spectra via the small matrices (cyclicity): A ~ N (K(w - z))_{z,w}; P, Q likewise.
# Checks on random configurations: tr A = N, tr A^2 = sum nu nu K(z-w)^2 (complex K),
# A Hermitian (spectrum real), each Q has <= p positive eigenvalues, and Lemma 2.1:
# tr A^2 >= 2c tr A - 2c tr P + tr phi_c(P) - c^2 p.
import json, numpy as np
rng = np.random.default_rng(7)
d = json.load(open('src/anc/certificate/inputs.json')); a = np.array(d['window_coefficients'],float)/1e6
Mq = 801; u = np.linspace(-0.5,0.5,Mq); wq = np.full(Mq, u[1]-u[0]); wq[0]=wq[-1]=wq[0]/2
f = 1+sum(a[j]*np.cos(2*np.pi*(j+1)*u) for j in range(6))
def Kmat(z, w):  # K(w_j - z_i) by quadrature, complex arguments allowed
    D = w[None,:] - z[:,None]
    return np.einsum('k,ijk->ij', wq*f, np.exp(-2j*np.pi*u[None,None,:]*D[:,:,None]))
def phi(l,c): return np.where(l<=c, l*l, 2*c*l-c*c)
def spec(z, nu):  # nonzero spectrum of sum nu_z F_z <F_{conj z}, .> = eig of diag(nu) (<F_{conj z_i}, F_{z_j}>) = diag(nu) K(z_j - z_i)
    if len(z)==0: return np.zeros(0)
    return np.linalg.eigvals(np.diag(nu) @ Kmat(z, z))
worst = np.inf; maxim = 0
for trial in range(400):
    c = float(rng.choice([2.0, 3.35])); cap = 2 if c > 2 else 1
    Z=[]; NU=[]; KIND=[]
    y = np.cumsum(0.3 + rng.exponential(0.9, size=rng.integers(3,10)))
    for yy in y:
        if rng.random() < 0.2:
            bb = float(rng.choice([1e-3, 0.05, 0.3, 1.0, 2.5]))*rng.random(); nu = int(rng.choice([1,1,1,2]))
            Z += [yy+1j*bb, yy-1j*bb]; NU += [nu,nu]; KIND += ['off','off']
        else:
            Z.append(yy+0j); NU.append(int(rng.choice([1,1,1,2,2,3]))); KIND.append('on')
    Z=np.array(Z); NU=np.array(NU,float); KIND=np.array(KIND)
    N = NU.sum()
    lamA = spec(Z, NU); trA = lamA.sum(); trA2 = (lamA**2).sum()
    maxim = max(maxim, np.abs(lamA.imag).max()/max(1,np.abs(lamA).max()))
    pred2 = np.sum(NU[:,None]*NU[None,:]*Kmat(Z,Z)**2)
    inP = (KIND=='on') & (NU<=cap)
    lamP = spec(Z[inP], NU[inP]).real; lamQ = spec(Z[~inP], NU[~inP])
    p = int(((KIND=='on') & (NU>cap)).sum() + (KIND=='off').sum()//2)
    npos = int((lamQ.real > 1e-9*max(1,np.abs(lamQ).max() if lamQ.size else 1)).sum())
    lower = 2*c*trA.real - 2*c*lamP.sum() + phi(np.clip(lamP,0,None),c).sum() - c*c*p
    worst = min(worst, trA2.real - lower)
    assert abs(trA - N) < 1e-8*N, (trA, N)
    assert abs(trA2 - pred2) < 1e-7*max(1,abs(pred2)), (trA2, pred2)
    assert npos <= p, (npos, p, lamQ)
    assert trA2.real >= lower - 1e-7*max(1,abs(trA2)), (trA2, lower)
print('400 random configurations (on-line mult 1-3, off-line pairs with Im z up to 2.5, mult 1-2): all checks passed')
print('max relative imaginary part of spec(A) = %.2e (A Hermitian => 0)' % maxim)
print('min over trials of trA^2 - Lemma 2.1 lower bound = %.4g (>= 0 required)' % worst)

