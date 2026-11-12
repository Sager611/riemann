
import json, numpy as np, sys
def yamada_check(fn, Lmax=12.0, dL=0.0005):
    d=json.load(open(fn)); h=d['h']; eta=np.array(d['eta']); K=len(eta)
    w = h*(1+eta)              # atom weights at kh, k=1..K (beyond K: weight h)
    ks = np.arange(1,K+1)*h
    Ls = np.arange(dL, Lmax, dL)
    # E N_L^2 = L + 2 sum_k w_k (L-kh)_+ ; for L < K h only finite part matters
    EN2 = Ls + 2*np.array([np.sum(w*np.clip(L-ks,0,None)) for L in Ls])
    var = EN2 - Ls**2
    th = Ls-np.floor(Ls); yam = th*(1-th)
    slack = var - yam
    i = np.argmin(slack)
    print(fn.split('/')[-1], 'lam',d['lam'],'theta',d['theta'],'h',round(h,5),'X',d['X'])
    print('  min Var-theta(1-theta) = %.5f at L=%.4f (Var=%.5f, Yamada=%.5f)'%(slack[i],Ls[i],var[i],yam[i]))
    neg = Ls[slack< -1e-9]
    if len(neg): print('  VIOLATED on L in [%.4f, %.4f], count %d grid pts; segments:'%(neg.min(),neg.max(),len(neg)))
    # print segments
    segs=[]; 
    if len(neg):
        start=neg[0]; prev=neg[0]
        for x in neg[1:]:
            if x-prev>1.5*dL: segs.append((start,prev)); start=x
            prev=x
        segs.append((start,prev))
        for s in segs[:20]:
            j0=np.argmin(np.where((Ls>=s[0]-1e-12)&(Ls<=s[1]+1e-12),slack,9))
            print('    [%.4f,%.4f] min slack %.5f at L=%.4f'%(s[0],s[1],slack[j0],Ls[j0]))
    # also print Var at some L
    for L in [0.5,1.0,1.25,1.5,1.75,2.0,2.5,3.0,4.0,5.0,8.0]:
        j=int(round(L/dL))-1; print('   L=%.2f Var=%.4f  Yamada=%.4f  GUE~'%(Ls[j],var[j],yam[j]))
for fn in sys.argv[1:]: yamada_check(fn)
