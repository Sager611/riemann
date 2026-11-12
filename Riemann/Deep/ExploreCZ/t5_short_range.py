
import numpy as np
rng=np.random.default_rng(7)
N=1000; n=np.arange(N+1,2*N+1); ln=np.log(n.astype(float))
def gram(W):
    E=np.exp(1j*np.outer(W,ln))          # R x N
    return E@E.conj().T                   # G_rs = sum_n n^{i(t_r-t_s)}
for (R,Lint,kind) in [(400,0.9*2*np.pi*N,'random'),(400,0.9*2*np.pi*N,'AP spacing ~14'),(900,0.9*2*np.pi*N,'random dense')]:
    if kind.startswith('AP'): W=np.arange(R)*(Lint/R)
    else:
        W=np.sort(rng.uniform(0,Lint,R*3)); keep=[W[0]]
        for x in W[1:]:
            if x-keep[-1]>=1.0: keep.append(x)
        W=np.array(keep[:R])
    Gm=gram(W); lam=np.linalg.eigvalsh(Gm).max()
    off=np.abs(Gm-np.diag(np.diag(Gm))).sum(axis=1).max()
    print(f"{kind:16s} R={len(W)} in [0,{Lint:.0f}]: lambda_max/N={lam/N:.3f}  Schur/triangle (N+max row sum)/N={(N+off)/N:.2f}  MV-type (|I|/1+N)/N... large sieve (Lint+N)/N={(Lint+N)/N:.2f}")
