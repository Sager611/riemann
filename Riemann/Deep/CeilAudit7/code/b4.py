# ceilaudit7 B: (1) CUE one-point 4th-order ratio R = <Z^2 Z_x^2>/<Z^4> (normalised units, 1 zero per unit) via Toeplitz;
# (2) LP: block mixtures of band-limited Gaussians, energy-weighted spectrum exactly flat on [-pi,pi] (second order = class S of explore60),
#     matching <X^4>/<X^2>^2 = M4, <X^2X'^2>/<X^4> = R, <|X|>/<X^2>^(1/2) = M1; minimise the zero density (Rice per component, length-weighted).
import sys, numpy as np
from scipy.optimize import linprog
def F(N, phi):
    q = np.exp(-1j*phi); f = {2: q, -2: np.conj(q), 1: -2*q-2, -1: -2*np.conj(q)-2, 0: 4+2*np.cos(phi)}
    M = N; T = np.zeros((M, M), complex); idx = np.arange(M)
    for m, v in f.items():
        j = idx[(idx-m >= 0) & (idx-m < M)]; T[j, j-m] = v
    s, ld = np.linalg.slogdet(T); return ld, s
mode = sys.argv[1]
if mode == 'ratio':
    for N in (100, 200, 400, 800):
        l0, _ = F(N, 0.0); ex = np.log((N+1)*(N+2)**2*(N+3)/12)
        out = []
        for h in (0.02, 0.01):
            lh, _ = F(N, 2*np.pi*h/N); out.append(-2*(np.expm1(lh-l0))/h**2/4)
        R = (4*out[1]-out[0])/3
        print(f'N={N} logF(0)-exact={l0-ex:+.2e} R(h=.02)={out[0]:.6f} R(h=.01)={out[1]:.6f} Richardson R={R:.6f} (Gaussian pi^2/9={np.pi**2/9:.6f}); M4_CUE={(N+2)**2*(N+3)/(12*(N+1)):.2f}')
elif mode == 'lp':
    M4, R, M1 = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    J = int(sys.argv[5]); edges = np.linspace(0, np.pi, J+1)
    amps = np.exp(np.linspace(np.log(1e-3), np.log(1e3), 61))
    cols = []; meta = []
    for j in range(J):
        w1, w2 = edges[j], edges[j+1]; l2 = (w2**3-w1**3)/(3*(w2-w1)); rho = np.sqrt(l2)/np.pi
        for A in amps:
            band = np.zeros(J); band[j] = A/(1.0/J)  # energy density per band, normalised so flat = 1 in every band
            cols.append(np.concatenate([[1.0], band*(1.0/J)*J, [3*A**2, A**2*l2, np.sqrt(2/np.pi)*np.sqrt(A)]])); meta.append((j, A, rho))
    V = np.array(cols).T  # rows: length, J band energies (target 1 each -> sum p A = 1 per band... see below), M4, M22, M1
    # band energy row j: sum_k p_k A_k [band_k = j] must equal 1/J (flat spectrum: each band carries 1/J of the energy)
    Aeq = []; beq = []
    for j in range(J):
        Aeq.append([meta[k][1] if meta[k][0] == j else 0.0 for k in range(len(meta))]); beq.append(1.0/J)
    Aeq.append([3*meta[k][1]**2 for k in range(len(meta))]); beq.append(M4)
    l2s = [((edges[m[0]+1]**3-edges[m[0]]**3)/(3*(edges[m[0]+1]-edges[m[0]]))) for m in meta]
    Aeq.append([meta[k][1]**2*l2s[k] for k in range(len(meta))]); beq.append(R*M4/3*3/3)  # <X^2X'^2>/<X^2>^2 = R*M4
    beq[-1] = R*M4
    Aeq.append([np.sqrt(2/np.pi)*np.sqrt(meta[k][1]) for k in range(len(meta))]); beq.append(M1)
    Aub = [[1.0]*len(meta)]; bub = [1.0]
    c = [m[2] for m in meta]
    res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=np.array(Aeq), b_eq=np.array(beq), bounds=[(0, None)]*len(meta), method='highs')
    print(f'M4={M4} R={R} M1={M1} J={J}: status={res.status} min zero density={res.fun if res.status==0 else None} (1/sqrt3={1/np.sqrt(3):.4f})')
    if res.status == 0:
        x = res.x; used = np.where(x > 1e-12)[0]
        print(f'  live length={x.sum():.4f}; components used:')
        for k in used: print(f'   band {meta[k][0]} A={meta[k][1]:.4g} p={x[k]:.4g} rho={meta[k][2]:.3f}')
        print(f'  check: energy per band={[round(sum(x[k]*meta[k][1] for k in used if meta[k][0]==j),6) for j in range(J)]}; M4={sum(3*x[k]*meta[k][1]**2 for k in used):.4f}; M22/M4={sum(x[k]*meta[k][1]**2*l2s[k] for k in used)/sum(x[k]*meta[k][1]**2 for k in used)*3/3:.4f}')
