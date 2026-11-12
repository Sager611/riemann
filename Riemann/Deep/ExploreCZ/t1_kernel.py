
import numpy as np, time, sys
t0=time.time()
def S_on(N, u, chunk=4000):
    n=np.arange(N+1,2*N+1,dtype=np.float64); ln=np.log(n)
    out=np.empty(len(u),dtype=np.complex128)
    for i in range(0,len(u),chunk):
        uu=u[i:i+chunk]
        out[i:i+chunk]=np.exp(1j*np.outer(uu,ln)).sum(axis=1)
    return out
for N in (300,1000):
    U=2*np.pi*N
    u=np.arange(1.0, 12.5*U, 0.25)
    S=np.abs(S_on(N,u))
    cz=u<U*0.9
    print(f"N={N}: CZ range u in [1,0.9*2piN]: max |S|*u/N = {np.max(S[cz]*u[cz]/N):.3f}; max|S| on [N,0.9*2piN] = {np.max(S[(u>N)&cz]):.2f} (sqrt N={np.sqrt(N):.1f})")
    rows=[]
    for m in range(1,13):
        b=(u>=U*(m-0.5))&(u<U*(m+0.5))
        mx=np.max(S[b]); um=u[b][np.argmax(S[b])]
        rows.append((m,mx/np.sqrt(N),mx*um/N, np.sqrt(np.mean(S[b]**2))/np.sqrt(N)))
    print("  band m (u~2piNm): max|S|/sqrtN, CZ-size ratio max|S|*u/N, rms|S|/sqrtN")
    for m,a,b_,c in rows: print(f"   m={m:2d}  {a:6.3f}  {b_:9.1f}  {c:6.3f}")
print("time",round(time.time()-t0,1))
