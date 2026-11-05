# ba_ratio.py: exact B/A at interior critical roots vs the split-contour claims of Thm 6.1-6.3/Lemma 6.4
# B/A = 2 x_m J''_{d,k}(x_m) / (d J_{d,k}(x_m))   (from Thm 3.2: A ~ J_{d,k}(x_m), B = z f_k'(z) ~ J'_{d-1,k+1}(x_m))
import sys, json, math, mpmath as mp
mp.mp.dps=60
g=[mp.mpf(s) for s in json.load(open(sys.argv[1]))]
R=json.load(open(sys.argv[2])); S=json.load(open(sys.argv[3]))
def J(d,k,x,der=0):
    c=[mp.binomial(d,j)*g[k+j] for j in range(d+1)]
    return mp.polyval(list(reversed(c)) if False else c[::-1], x) if der==0 else mp.diff(lambda t: mp.polyval(c[::-1],t), x, der)
dstar=0.09756
out=[]
for key,val in R.items():
    d,k=map(int,key.split(','))
    for m in range(1,d-2):
        z=mp.mpf(val['z'][m-1]); x=-z*z/2
        c=[mp.binomial(d,j)*g[k+j] for j in range(d+1)]
        P=lambda t: mp.polyval(c[::-1],t)
        Jx=P(x); J2=mp.diff(P,x,2)
        BA=2*x*J2/(d*Jx)
        s=S.get(f"{d},{k},{m}")
        if s is None: continue
        sig,tau=s['sigma'],s['tau']; Phia=s['Phia']; thu=math.atan2(tau,sig)
        Lpos=(sig*sig+tau*tau)*math.sin(Phia+dstar+thu)/math.sin(Phia+dstar-thu)
        out.append((d,k,m,float(BA),Lpos,float(BA)/Lpos))
print('anchor (9,0): m, B/A, Lambda_pos(actual saddle), ratio')
for o in out:
    if o[0]==9 and o[1]==0: print(o[2], round(o[3],5), round(o[4],5), round(o[5],4))
mn=min(out,key=lambda o:o[5]); mx=max(out,key=lambda o:o[5]); print('min (B/A)/Lambda_pos over grid:', mn); print('max (B/A)/Lambda_pos over grid:', mx)
viol=[o for o in out if o[5] < 1-0.853]
print('count with B/A < (1-0.853)*Lambda_pos :', len(viol), viol[:8], 'PASS' if not viol else 'FAIL')
viol2=[o for o in out if o[3] < 0.72231-0.49688 and o[0]==9 and o[1]==0]
print('(9,0) roots with B/A < 0.2254 (anchor claim):', viol2)
