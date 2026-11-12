import numpy as np, math, sys
d = np.load(sys.argv[1]); xs = d['xs']; L = float(d['L']); R = float(d['R']); T0 = float(d['T0']); dt = float(d['dt'])
logV = d['logV'].astype(float); logP = d['logP'].astype(float); J = logV.shape[1]; T1 = T0 + dt*(J-1)
th = lambda t: t/2*np.log(t/(2*np.pi)) - t/2 - np.pi/8 + 1/(48*t) + 7/(5760*t**3)
Nz = (th(T1) - th(T0))/np.pi
tt = T0 + dt*np.arange(J); Z = np.real(np.exp(1j*th(tt))*d['Fhalf'])
nsc = int(np.sum(np.sign(Z[1:]) != np.sign(Z[:-1])))
kn = 2*np.pi*Nz/((T1-T0)*L)
def winding(Vline, Vb, Vt):
    # counterclockwise: bottom (t=T0, sig 1/2->3), right (sig=3), top (t=T1, 3->1/2), left (sig=1/2, T1->T0)
    ab = np.unwrap(np.angle(Vb)); at = np.unwrap(np.angle(Vt[::-1])); al = np.unwrap(np.angle(Vline[::-1]))
    right = np.angle(Vt[-1]/Vb[-1])     # sigma=3 edge, |V-1|<1 there
    tot = (ab[-1]-ab[0]) + right + (at[-1]-at[0]) + (al[-1]-al[0])
    jump = np.max(np.abs(np.diff(np.unwrap(np.angle(Vline)))))
    return tot/(2*np.pi), jump
nV, jV = winding(d['Vhalf'], d['Vb'], d['Vt'])
nP, jP = winding(d['Phalf'], d['Pb'], d['Pt'])
print(f"T0={T0} T1={T1:.1f} L={L:.4f} a={float(d['a']):.4f} R={R} N_zeta(RvM main)={Nz:.1f} signchanges Z={nsc} kn={kn:.4f}")
print(f"zeros right of 1/2 in [T0,T1]: V {nV:.3f} (max arg step {jV:.2f}), Psi {nP:.3f} (step {jP:.2f}); fraction nu_V={nV/Nz:.4f} nu_Psi={nP/Nz:.4f}")
g = np.euler_gamma
print(" x      sigma   m_f     h_f     J_f    var_f   J/var | m_V     h_V     J_V  | m_P     h_P    J_P  | h_f/(|x|kn)")
rows = {}
for i, x in enumerate(xs):
    lf = logV[i] + logP[i]
    m_f = lf.mean(); mx = lf.max(); h_f = 0.5*(np.log(np.mean(np.exp(2*(lf-mx)))) + 2*mx)
    m_V = logV[i].mean(); mv = logV[i].max(); h_V = 0.5*(np.log(np.mean(np.exp(2*(logV[i]-mv)))) + 2*mv)
    m_P = logP[i].mean(); mp = logP[i].max(); h_P = 0.5*(np.log(np.mean(np.exp(2*(logP[i]-mp)))) + 2*mp)
    v = lf.var(); rows[x] = (m_f, h_f, m_V, h_V, m_P, h_P)
    lp = h_f/(abs(x)*kn) if x < 0 else float('nan')
    print(f"{x:6.3f} {0.5+x/L:6.4f} {m_f:7.4f} {h_f:7.4f} {h_f-m_f:6.4f} {v:6.4f} {(h_f-m_f)/v:5.3f} | {m_V:7.4f} {h_V:7.4f} {h_V-m_V:6.4f} | {m_P:7.4f} {h_P:7.4f} {h_P-m_P:6.4f} | {lp:.4f}")
iR = int(np.argmin(abs(xs + R))); i0 = int(np.argmin(abs(xs)))
m_f = lambda i: logV[i].mean() + logP[i].mean()
hR = rows[xs[iR]][1]; mR = rows[xs[iR]][0]; m0 = rows[xs[i0]][0]; nf = (nV + nP)/Nz; nuV = nV/Nz
A = (hR - mR)/(R*kn); B = (mR - m0)/(R*kn) - nf; C = m0/(R*kn); D = nf - nuV
print(f"Levinson finite-T bound h_f(-R)/(R kn) = {hR/(R*kn):.4f}; actual nu_V(0) = {nuV:.4f}; nu_f(0) = {nf:.4f}")
print(f"decomposition: A concavity {A:.4f} + B weights (a<beta<1/2, beta>1/2 overweight) {B:.4f} + C m_f(0) {C:.4f} + D mollifier zeros {D:.4f} = {A+B+C+D:.4f}")
print(f"LP optimum over grid lines x<0: {min(rows[x][1]/(abs(x)*kn) for x in xs if x < 0):.4f}; m at last line (sigma={0.5+xs[-1]/L:.2f}) = {rows[xs[-1]][0]:.5f}")
print(f"kappa-equivalents 1-2nu: Levinson bound {1-2*hR/(R*kn):.4f}, truth with V {1-2*nuV:.4f}; Gaussian Jensen gap gamma/2 = {g/2:.4f}")
