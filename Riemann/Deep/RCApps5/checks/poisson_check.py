#!/usr/bin/env python3
# Lemma 6 of RCApps3 (Poisson in T) against brute force, for M with square factors p^2 || M (the pairs n = e m, n' = e m'
# that Lemma P keeps in the w = 1 count), X' = X, several odd squarefree c prime to M.
import math, signal
signal.alarm(600)
def jacobi(a, n):
    a %= n; r = 1
    while a:
        while a%2==0:
            a//=2
            if n%8 in (3,5): r=-r
        a,n = n,a
        if a%4==3 and n%4==3: r=-r
        a %= n
    return r if n==1 else 0
def mob(n):
    r=1;p=2
    while p*p<=n:
        if n%p==0:
            n//=p
            if n%p==0: return 0
            r=-r
        p+=1
    return -r if n>1 else r
def gcd(a,b):
    while b: a,b=b,a%b
    return a
def step(x):  # smooth 0->1 on [0,1]
    if x<=0: return 0.0
    if x>=1: return 1.0
    f=lambda t: math.exp(-1/t) if t>0 else 0.0
    return f(x)/(f(x)+f(1-x))
def Phi(u): # C^infty, supp [1/4,2], =1 on [1/2,1]
    return step((u-0.25)/0.25)*step((2-u)/1.0)
def K(v, n=4000):
    # K(v) = v^{-1} int_0^1 Phi(s/v)(1-s)^{-1/2} ds = v^{-1} int_0^1 2 Phi((1-t^2)/v) dt
    h=1.0/n; s=0.0
    for i in range(n):
        t=(i+0.5)*h; s+=2*Phi((1-t*t)/v)
    return s*h/v
def A(M, y):
    q=4*y*y; cnt=0
    for T in range(q):
        d=T*T-4*M
        if d%(y*y): continue
        if (d//(y*y))%4 in (0,1): cnt+=1
    return cnt/4
def brute(M, X, c):
    tot=0.0
    ymax=int(math.sqrt(4*M/(X/4)))+1
    for y in range(1,ymax+1):
        # 4M - T^2 = |Delta| y^2 in [X/4, 2X] y^2
        lo=4*M-2*X*y*y; hi=4*M-X*y*y/4
        if hi<0: continue
        Tmax=int(math.isqrt(max(0,int(hi))))+1; Tmin=int(math.isqrt(max(0,int(lo))))-1 if lo>0 else 0
        for T in range(max(0,Tmin), Tmax+1):
            d=T*T-4*M
            if d>=0 or d%(y*y): continue
            D=d//(y*y)
            if D%4 not in (0,1): continue
            if (T-y*D)%2: continue
            w=Phi(-D/X)
            if w==0: continue
            mult=(1 if T==0 else 2)*2   # +-T, +-y
            tot+=mult*w*jacobi(D,c)
    return tot
def main_term(M, X, c):
    s=0.0; y=1
    while True:
        v=y*y*X/(4*M)
        if v>8.5: break
        if gcd(y,c)==1:
            k=K(v)
            if k!=0: s+=A(M,y)*k
        y+=1
    return mob(c)/c*X/math.sqrt(M)*s
cases=[(19*19*23*29,  4*10**5, 1),(19*19*23*29, 4*10**5, 3*5),(19*19*23*29*31, 2*10**6, 1),(23*23*29*31*37,3*10**6,1),
       (23*23*29*31*37, 3*10**6, 7),(19*19*23*23*29, 2*10**6, 1),(19*19*23*23*29, 2*10**6, 3),(41*43*47*53, 2*10**6, 1)]
for M,X,c in cases:
    b=brute(M,X,c); m=main_term(M,X,c)
    print(f"M={M} X={X} c={c} X/(c sqrt M)={X/(c*math.sqrt(M)):.1f}  brute={b:.4f}  main={m:.4f}  diff={b-m:.4f}")

