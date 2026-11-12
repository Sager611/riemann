
import re, math, sys
import numpy as np
rows = []
for line in open(sys.argv[1]):
    m = re.match(r'(-?\d+) (\d+) (\d+) (-?\d+) \[(.*)\] (\d+)', line.strip())
    if not m: continue
    D, h, L, feq = int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4))
    vals = [float(x.replace(' E', 'e').replace('E','e')) for x in m.group(5).split(',')]
    rows.append((abs(D), h, feq, vals))
rows.sort()
labels = ['(0,0)','(0,1.5)','(0,4)','(.05,0)','(.05,1.5)','(.05,4)','(.15,0)','(.15,1.5)','(.15,4)']
print('n D =', len(rows), ' |D| range', rows[0][0], rows[-1][0], ' worst feq bits', max(r[2] for r in rows))
wins = [(12,25),(25,50),(50,100),(100,200),(200,400),(400,800)]
W = []
for a,b in wins:
    sel = [r for r in rows if a < r[0] <= b]
    if not sel: continue
    H = sum(r[1] for r in sel)
    S = [sum(r[3][k] for r in sel) for k in range(9)]
    W.append((a,b,len(sel),H,S))
print('window  #D  sum_h   S/sum_h for (delta,t)=', labels)
for a,b,n,H,S in W:
    print(f'({a},{b}] {n:3d} {H:5d} ' + ' '.join(f'{s/H:7.2f}' for s in S))
print('exponent of S(X) in X between consecutive windows (S ~ X^e); family size exponent for comparison:')
for i in range(1,len(W)):
    a0,b0,_,H0,S0 = W[i-1]; a1,b1,_,H1,S1 = W[i]
    r = math.log(b1/b0)
    print(f'({a0},{b0}]->({a1},{b1}]: sum_h {math.log(H1/H0)/r:5.2f} | ' + ' '.join(f'{math.log(S1[k]/S0[k])/r:5.2f}' for k in range(9)))
# per-D normalized R_D = S_D/h ; regress log R_D on log|D| and on log log|D|
X = np.array([r[0] for r in rows if r[0] > 50]); 
for k in [1,2,4,5,6,7,8]:
    R = np.array([r[3][k]/r[1] for r in rows if r[0] > 50])
    ok = R > 0
    A = np.vstack([np.log(X[ok]), np.ones(ok.sum())]).T
    e, c = np.linalg.lstsq(A, np.log(R[ok]), rcond=None)[0]
    # fit R = a + b log X
    B = np.vstack([np.log(X[ok]), np.ones(ok.sum())]).T
    b2, a2 = np.linalg.lstsq(B, R[ok], rcond=None)[0]
    print(f'{labels[k]:9s} per-D R_D=S_D/h: power fit exponent {e:5.2f}; linear-in-log fit R ~ {a2:6.2f} + {b2:5.2f} log|D|; max R_D/|D|^0.5 = {np.max(R/np.sqrt(X)):.3f}')
