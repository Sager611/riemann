
import json, numpy as np, sys
def load(fn):
    d=json.load(open(fn)); h=d['h']; eta=np.array(d['eta']); K=len(eta); k0=d['k0']
    Kx = 3*K
    mu = np.full(Kx+1, h); mu[0]=1.0; mu[1:K+1]=h*(1+eta)
    return d,h,k0,K,mu
for fn in sys.argv[1:]:
    d,h,k0,K,mu=load(fn)
    print(fn.split('/')[-1],'lam',d['lam'],'theta',round(d['theta'],5),'h',round(h,5),'k0',k0,'K',K)
    big = np.argsort(-mu[1:K+1])[:8]+1
    print('  largest atoms (u, mu):', [(round(k*h,3), round(mu[k],3)) for k in big])
    # triangle: mu_p+mu_q - mu_{p+q} <= 1, mu_p + mu_q - mu_{|p-q|} <= 1
    best1=(-9,None); best2=(-9,None)
    for p in range(1,K+1):
        q=np.arange(p,K+1)
        v1 = mu[p]+mu[q]-mu[p+q]-1
        i=np.argmax(v1)
        if v1[i]>best1[0]: best1=(v1[i],(p,q[i]))
        v2 = mu[p]+mu[q]-mu[q-p]-1
        v2[q==p] = -9
        i=np.argmax(v2)
        if v2[i]>best2[0]: best2=(v2[i],(p,q[i]))
    print('  max mu_p+mu_q-mu_{p+q}-1 = %.4f at (p,q)=%s u=(%.3f,%.3f)'%(best1[0],best1[1],best1[1][0]*h,best1[1][1]*h))
    print('  max mu_p+mu_q-mu_{q-p}-1 = %.4f at (p,q)=%s u=(%.3f,%.3f)'%(best2[0],best2[1],best2[1][0]*h,best2[1][1]*h))
    # packing: windows of k0 consecutive sites p..p+k0-1 (p>=1): sum <= 1
    cs=np.concatenate([[0],np.cumsum(mu[1:])])
    w=np.array([cs[p+k0-1]-cs[p-1] for p in range(1,K+1)])
    i=np.argmax(w); print('  max Palm mass in a window of length < lam: %.4f at u in [%.3f,%.3f]'%(w[i],(i+1)*h,(i+k0)*h))
    # count of violated triangle pairs
    nv=0
    for p in range(1,K+1):
        q=np.arange(p,K+1); nv+=np.sum(mu[p]+mu[q]-mu[p+q]-1>1e-9)
    print('  # violated (p,q) sum-type triangles:', nv)
