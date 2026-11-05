# roots_floor.py: roots of J^{d-1,k+1}, rescaled frequencies omega_m, check Thm 5.1/5.2 floor
# Run: uv run --no-project --with mpmath python3 roots_floor.py GAMJSON DMAX KMAX OUTJSON
import sys, json, mpmath as mp
mp.mp.dps=60
g=[mp.mpf(s) for s in json.load(open(sys.argv[1]))]
DMAX=int(sys.argv[2]); KMAX=int(sys.argv[3]); OUT=sys.argv[4]
mu=lambda n: g[n-1]/g[n]
delta=lambda n: 1-mu(n)/mu(n+1)
def OmegaLS(d,k):
    return mp.sqrt((2*d+k-mp.mpf(1)/2)*mu(k+2)/(1+(d-2)*mp.sqrt(delta(k+2))))
res={}
fails=[]; minfloor=(1e9,None)
for d in range(4,DMAX+1):
    for k in range(0,KMAX+1):
        if d+k<9 or d+k+1>len(g)-1: continue
        N=d-1
        coeffs=[mp.binomial(N,j)*g[k+1+j] for j in range(N+1)]  # ascending
        r=mp.polyroots(coeffs[::-1], maxsteps=400, extraprec=200)
        xs=sorted([mp.re(x) for x in r], reverse=True)  # closest to 0 first
        imax=max(abs(mp.im(x)) for x in r)
        M=4*d+2*k-1
        z=[mp.sqrt(-2*x) for x in xs]
        om=[z_*mp.sqrt(M)/2 for z_ in z]
        OL=OmegaLS(d,k)
        ok = om[0] >= OL and OL >= mp.mpf('12.2757')
        if not ok: fails.append((d,k,float(om[0]),float(OL)))
        if om[0]<minfloor[0]: minfloor=(float(om[0]),(d,k))
        res[f"{d},{k}"]={'z':[mp.nstr(v,30) for v in z],'omega':[mp.nstr(v,20) for v in om],'OmegaLS':mp.nstr(OL,15),'maximag':mp.nstr(imax,5)}
json.dump(res,open(OUT,'w'))
print('pairs',len(res),'floor fails',fails[:10])
print('min omega_1 over grid',minfloor)
print('Thm 5.1-5.2 floor omega_1>=Omega_LS>=12.2757 on grid:', 'PASS' if not fails else 'FAIL')
for key in ['9,0','10,0','15,0','20,0','25,0','30,0','4,5','5,4','8,1']:
    if key in res:
        o=res[key]['omega']; print(key,'omega_1=',o[0][:10],'omega_{d-3}=',o[-3][:10] if len(o)>=3 else None,'omega_{d-1}=',o[-1][:10],'OmegaLS=',res[key]['OmegaLS'][:10])
