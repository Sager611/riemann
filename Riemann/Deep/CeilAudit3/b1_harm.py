
import numpy as np, math, sys
sys.argv=["x","0"]
exec(open("b1_cue.py").read().split("designs=[(60,8)")[0])
for (K,m) in [(60,8),(120,16),(190,16)]:
    v,B,(lo,hi)=design(K,m,us)
    A=[2*pr*math.cosh(k*pr*Y0)*float(v(np.array([k*us]))[0]) for k in (1,2,3,4)]
    # D-mean part bound: sup over j>=K of j!(2/x)^j-type terms is ~0; report v(0+) behaviour
    print(K,m,"Band",round(B,2),"A1..A4:",[f"{a:.3e}" for a in A],"ratio(A1)",round(A[0]/B,4),"|A2+A3+A4|/A1",f"{sum(abs(a) for a in A[1:])/abs(A[0]):.2e}")
