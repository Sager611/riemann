
import re, math
import numpy as np
def load(fn, idx):
    rows=[]
    for line in open(fn):
        m = re.match(r'(-?\d+) (\d+) (\d+) (-?\d+) \[(.*)\] (\d+)', line.strip())
        if not m: continue
        vals=[float(x.replace(' E','e').replace('E','e')) for x in m.group(5).split(',')]
        rows.append((abs(int(m.group(1))), int(m.group(2)), int(m.group(4)), [vals[i] for i in idx]))
    return rows
rows = load('msq.out',[0,1,3,4,6,7]) + load('msq2.out',[0,1,2,3,4,5])
rows.sort()
labels=['(0,0)','(0,1.5)','(.05,0)','(.05,1.5)','(.15,0)','(.15,1.5)']
print('D count', len(rows), 'max |D|', rows[-1][0], 'worst feq bits', max(r[2] for r in rows))
wins=[(37.5,75),(75,150),(150,300),(300,600)]
W=[]
for a,b in wins:
    sel=[r for r in rows if a<r[0]<=b]; H=sum(r[1] for r in sel); S=[sum(r[3][k] for r in sel) for k in range(6)]
    W.append((a,b,len(sel),H,S))
print('window  #D sum_h | S/sum_h for', labels)
for a,b,n,H,S in W: print(f'({a:g},{b:g}] {n:3d} {H:5d} | '+' '.join(f'{s/H:7.2f}' for s in S))
print('window exponents d log S / d log X (family size sum_h first):')
for i in range(1,len(W)):
    a0,b0,_,H0,S0=W[i-1]; a1,b1,_,H1,S1=W[i]; r=math.log(b1/b0)
    print(f'({a0:g},{b0:g}]->({a1:g},{b1:g}]: sum_h {math.log(H1/H0)/r:5.2f} | '+' '.join(f'{math.log(S1[k]/S0[k])/r:5.2f}' for k in range(6)))
# overall fit of S over the windows
Xs=np.array([b for a,b,_,_,_ in W]); Hs=np.array([H for _,_,_,H,_ in W])
print('least-squares exponent over the 4 windows: sum_h', round(np.polyfit(np.log(Xs),np.log(Hs),1)[0],2), '| S:', ' '.join(str(round(np.polyfit(np.log(Xs),np.log([w[4][k] for w in W]),1)[0],2)) for k in range(6)))
print('ratio S/sum_h fitted as a + b log X vs c X^{1/2}: residual comparison')
for k in range(6):
    y=np.array([w[4][k]/w[3] for w in W]); L=np.log(Xs)
    A=np.vstack([L,np.ones_like(L)]).T; coef,res,_,_=np.linalg.lstsq(A,y,rcond=None); r1=np.sum((A@coef-y)**2)
    c=np.sum(y*np.sqrt(Xs))/np.sum(Xs); r2=np.sum((c*np.sqrt(Xs)-y)**2)
    print(f'  {labels[k]:9s} log-model a+b logX: a={coef[1]:6.2f} b={coef[0]:5.2f} rss={r1:7.3f} | sqrt-model c={c:5.3f} rss={r2:7.3f}')
