# Exact reproduction of PR137 Anchor/ChainData.lean recursion (176 stages), fractions.
from fractions import Fraction as F
import math, sys
STEP=F(1,4000); TOP=F(311,400); LEFT=F(1467,2000); DELTA=F(1,25000000); PAD=F(1,10**7)
def fH(R,r): return r/R
def fu(R,r): return (1-r/R)/6
def fY(R,z,k,r): return r/R*z+(1/k-2)*fu(R,r)
def fS(R,z,k,r): return (1-r)/2+fY(R,z,k,r)
def fs(R,z,k,r): return fS(R,z,k,r)-fu(R,r)
def fwb(R,z,a0,a1): return (R/2+z-a0)/(F(1,2)-a1)
def fthr(R,z,a0,a1):
    wb=fwb(R,z,a0,a1); return (5*wb/4-F(5,6))/(3*wb/2-(1+R)/2-z)
def fF0(R,z,k,a0,a1): return F(5,6)-(F(5,4)-3*k/2)*fwb(R,z,a0,a1)-k*(1+R)/2-k*z
def fg(R,z,k,a0,a1,r,nu): return (r+fs(R,z,k,r)-fH(R,r)/2-fH(R,r)*a0-fu(R,r)/nu)/(F(1,2)-a1)
def fw(R,z,k,a0,a1,r,nu): return fg(R,z,k,a0,a1,r,nu)/fH(R,r)+F(1,10**6)
def fsig(R,z,k,a0,a1,r):
    w=fw(R,z,k,a0,a1,r,F(3,4)); return (1-w)/2+a0-a1*w
def fline(a0,a1,r): return a0-a1*r
def qfloor(q): return q.numerator//q.denominator
def qceil(q): return -qfloor(-q)
def gates(R,z,k,a0,a1,r):
    S=fS(R,z,k,r); u=fu(R,r); H=fH(R,r); s0=fs(R,z,k,r); Y=fY(R,z,k,r)
    wlo=fw(R,z,k,a0,a1,r,k); whi=fw(R,z,k,a0,a1,r,F(3,4))
    g=[H*fF0(R,z,k,a0,a1)-F(3,100000), F(137,1000)-S, F(14,100)-S/H, F(9,10)-r-S-4*u, F(91,100)-r-S-4*u,
       F(5,6)+F(21,250)-r-S-3*u-F(1,100), 3*H-2*r-8*s0-F(1,5), r-F(1,2)-Y+u-F(1,5),
       r/2-s0-F(3,2)*(F(5,6)-k)*(u/(k*k))/(F(1,2)-a1)-F(1,10), H+u-H/6-F(1,15)-F(3,4)*r-s0-F(3,100),
       wlo-F(2,3), F(79,100)-whi, 1-F(3,4)*whi-F(3,4)*max(F(0),fsig(R,z,k,a0,a1,r))-F(1,6),
       F(79,100)-r-(2+F(1,10**11))*u]
    return g
class St:
    def __init__(s,R,z,a0,a1):
        s.R,s.z,s.a0,s.a1=R,z,a0,a1
        s.left=R-STEP; s.wb=fwb(R,z,a0,a1); s.thr=fthr(R,z,a0,a1)
        s.k=F(qceil((s.thr+F(1,5000))*10**6),10**6)
        s.F0=fF0(R,z,s.k,a0,a1)
    def Y(s,r): return fY(s.R,s.z,s.k,r)
    def nxt(s):
        y0=s.Y(s.left)-DELTA-2*PAD; y1=s.Y(s.R)-DELTA-2*PAD
        own1=-(y1-y0)/(s.R-s.left); own0=y0+own1*s.left
        ng=F(qfloor((s.Y(s.left)-DELTA)*10**10),10**10)
        s.own0,s.own1=own0,own1
        return St(s.left,ng,own0,own1)
def run(n=176):
    st=[St(TOP,F(0),F(0),F(0))]
    for i in range(n-1): st.append(st[-1].nxt())
    st[-1].nxt()
    return st
if __name__=='__main__':
    st=run()
    g0=F(552429,10**7); g1=F(7227847,10**8)
    worst=None; mins=[None]*14
    for i,s in enumerate(st):
        for r in (s.left,s.R):
            gg=gates(s.R,s.z,s.k,s.a0,s.a1,r)
            for j,v in enumerate(gg):
                if mins[j] is None or v<mins[j][0]: mins[j]=(v,i,float(r))
            m=s.Y(r)-DELTA-fline(g0,g1,r)
            if worst is None or m<worst[0]: worst=(m,i,float(r))
    for i in (0,1,2,50,100,150,174,175):
        s=st[i]; print("stage %3d R=%.6f gain=%.10f k=%.6f thr=%.8f wb=%.6f a0=%.6f a1=%.6f Y(left)=%.10f"%(i,float(s.R),float(s.z),float(s.k),float(s.thr),float(s.wb),float(s.a0),float(s.a1),float(s.Y(s.left))))
    print("left 175 =",st[175].left, st[175].left==LEFT)
    print("line g0-g1 r minorant margin min (Y-DELTA-line):",float(worst[0]),"stage",worst[1],"r",worst[2])
    for j,m in enumerate(mins): print("gate%2d min %.3e at stage %d r=%.5f"%(j+1,float(m[0]),m[1],m[2]))
