# explore61 s11: Theorem C numbers.  d_N <= |chi - h| + sum_j |e_j| dist(P_W1 rho_j, V_N)  (in-W1 power function), h = hybrid optimizer:
#   bound on P:  P_C = 2 D sqrt(Q_pw) + Q_pw,  Q_pw = (sum |e_j| sqrt((S-Om)_jj))^2  >= Q = e'(S-Om)e;
#   full-space power bound (no W1 projection): P_full = 2 D sqrt(F) + F,  F = (sum |e_j| sqrt(S_jj))^2  >= P.
import glob, math, numpy as np
for f in sorted(glob.glob('logs/s1_N*_Q*_K5.npz'), key=lambda f: (int(f.split('_N')[1].split('_')[0]), int(f.split('_Q')[1].split('_')[0]))):
    z = np.load(f); N = int(f.split('_N')[1].split('_')[0]); Qg = int(f.split('_Q')[1].split('_')[0]); L = math.log(N)
    S, Om, e = z['S'], z['Om'], z['e']; P, R = float(z['P']), float(z['R']); Q = P - R; D2 = float(z['d2']) - P; D = math.sqrt(D2)
    Qpw = float(np.sum(np.abs(e) * np.sqrt(np.clip(np.diag(S - Om), 0, None)))) ** 2
    F = float(np.sum(np.abs(e) * np.sqrt(np.diag(S)))) ** 2
    PC = 2 * D * math.sqrt(Qpw) + Qpw; Pf = 2 * D * math.sqrt(F) + F
    print(f'N={N:5d} Q={Qg:2d}: L^2P={L*L*P:.4f} | in-W1: Q_pw/Q={Qpw/Q:.2f} L^3 Q_pw={L**3*Qpw:.3f} => L^2 P_C={L*L*PC:.4f} (P_C/P={PC/P:.2f}) '
          f'| full space: F/P={F/P:.1f} L^2 P_full={L*L*Pf:.3f} (P_full/P={Pf/P:.1f}) | L^1.5 sqrt(Q_pw)={L**1.5*math.sqrt(Qpw):.3f}  L sqrt(F)={L*math.sqrt(F):.3f}')
