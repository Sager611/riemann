# explore61 s12: markdown table of the split for every s1 run (K0 = 5 head), for REPORT.md.
import glob, math, re, numpy as np
from scipy.linalg import eigh
zf = {}
for line in open('logs/s6_optzero.log'):
    m = re.search(r'N=(\d+) Q=(\d+):.*measured L\^3 Q=([0-9.]+)\s+predicted L\^2 sum\|G_e\|\^2=([0-9.]+)', line)
    if m: zf[(int(m.group(1)), int(m.group(2)))] = float(m.group(3)) / float(m.group(4))
rows = []
for f in glob.glob('logs/s1_N*_Q*_K5.npz'):
    z = np.load(f); N = int(f.split('_N')[1].split('_')[0]); Qg = int(f.split('_Q')[1].split('_')[0]); L = math.log(N)
    S, Om, e = z['S'], z['Om'], z['e']; P, R = float(z['P']), float(z['R']); Q = P - R; D2 = float(z['d2']) - P
    cosB = R / math.sqrt((D2 - R) * Q); lam = eigh(S - Om, Om, eigvals_only=True).max()
    Qpw = float(np.sum(np.abs(e) * np.sqrt(np.clip(np.diag(S - Om), 0, None)))) ** 2
    F = float(np.sum(np.abs(e) * np.sqrt(np.diag(S)))) ** 2
    rows.append((N, Qg, f'| {N} | {Qg} | {L*L*P:.4f} | {L*L*R:.4f} | {R/P:.3f} | {L**3*Q:.4f} | {L*Q/R:.3f} | {cosB:.3f} | {P/D2:.4f} ≤ {Q/R:.4f} | '
                 f'{zf.get((N, Qg), float("nan")):.3f} | {L*lam:.2f} | {Qpw/Q:.2f} | {F/P:.1f} | {L*np.abs(e).sum():.2f} |'))
print('| N | grid Q | L²P | L²R | R/P | L³Q | c_e = LQ/R | cos | P/D² ≤ Q/R (Lemma B) | L·Q / zero sum | L·λmax(S−Ω, Ω) | Q_pw/Q | F/P | L‖e‖₁ |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for r in sorted(rows): print(r[2])
