
# chamlu.py: spot check of lem:chamlu (max over sigma in [-1/64,1] of 1/|zeta(sigma+iT)|) and of the sigma<-1/64 bound.
import time
from flint import acb, ctx, arb
ctx.prec=96
paper={1e6+1:0.894198297,1893193.5:0.591342108,1e7+1:0.551087906,1e8+1:0.579764046,1e9:0.669256578,1e10+1:0.536165863}
for T,v in paper.items():
    t=time.time(); n=1200; best=(1e9,None)
    for k in range(n+1):
        s=-1/64+(1+1/64)*k/n
        z=abs(acb(s,T).zeta()); zm=float(z.mid())
        if zm<best[0]: best=(zm,s)
    # refine around best
    s0=best[1]; h=(1+1/64)/n
    for k in range(201):
        s=max(-1/64,min(1,s0-h+2*h*k/200)); zm=float(abs(acb(s,T).zeta()).mid())
        if zm<best[0]: best=(zm,s)
    # sigma < -1/64 samples
    neg=max(float((1/abs(acb(s,T).zeta())).mid()) for s in [-0.02,-0.1,-0.5,-1,-2,-5])
    print('T=%-12g max 1/|zeta| on grid = %.9f at sigma=%.5f ; paper %.9f ; ratio %.6f ; max over sample sigma<-1/64: %.4f ; %.1fs'%(T,1/best[0],best[1],v,(1/best[0])/v,neg,time.time()-t),flush=True)
