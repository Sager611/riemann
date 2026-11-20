
import numpy as np, math, sys
sys.argv=["x","0"]
exec(open("b1_cue.py").read().split("designs=[(60,8)")[0])
for (K,m) in [(60,8),(120,16),(190,16)]:
    v,B,(lo,hi)=design(K,m,us)
    uu=np.linspace(lo,hi,20001); du=uu[1]-uu[0]
    f=np.exp(Y0*uu/2)*v(uu)
    env=[]
    for k in [0,1,2,5,10,20,50]:
        dx=k*A
        env.append(abs(np.sum(f*np.exp(1j*uu*dx/2))*du))
    e0=env[0]
    print(K,m,"kernel envelope |K(kP)|/|K(0)| for k=0,1,2,5,10,20,50:",[f"{e/e0:.2e}" for e in env])
    # periodic row consistency: sum over k of K(kP) at frequency u*: coherent count
