"""d0check.py (ugoaudit): independent re-derivation of D0 (DERIVATION 3.2-3.5) with mpmath interval arithmetic (not flint).
uv run --no-project --with mpmath python3 d0check.py"""
from mpmath import iv, mp
from fractions import Fraction
iv.prec = 200
def primes(m):
    s=[1]*(m+1); s[0]=s[1]=0
    for i in range(2,int(m**0.5)+1):
        if s[i]:
            for j in range(i*i,m+1,i): s[j]=0
    return [i for i in range(m+1) if s[i]]
PR=primes(2000)
def T(k):            # 2^k prod_{i<=k} p_i/(p_i-1), exact
    t=Fraction(2**k)
    for p in PR[:k]: t*=Fraction(p,p-1)
    return t
def Pi(k):
    r=1
    for p in PR[:k]: r*=p
    return r
def g_lo(D,c):       # certified lower bound of g_c(D) = c sqrt D / (log(D+4) max(3, loglog D))
    D=iv.mpf(D); ll=iv.log(iv.log(D))
    mhi = max(iv.mpf(3).b, ll.b)          # upper end of max(3, loglog D)
    den = iv.log(D+4).b * mhi
    return (c*iv.sqrt(D)/den).a
def gt(D,c,t):       # certified g_c(D) > t
    return g_lo(D,c) > iv.mpf(t.numerator)/t.denominator
def Dplus(c,t):      # least integer D >= 16 with certified g_c(D) > t (g increasing on [16, oo))
    if gt(16,c,t): return 16
    hi=32
    while not gt(hi,c,t): hi*=2
    lo=hi//2
    while hi-lo>1:
        mid=(lo+hi)//2
        if gt(mid,c,t): hi=mid
        else: lo=mid
    return hi
def run(cstr):
    c=iv.mpf(cstr)
    rows=[]
    for k in range(1,61):
        t=T(k); d=Dplus(c,t); rows.append((k,d,Pi(k),Pi(k)>=d))
    poss=[r for r in rows if not r[3]]
    D0=max(r[1] for r in poss)
    kw=max(poss,key=lambda r:r[1])[0]
    K=min(k for k in range(15,61) if all(r[3] for r in rows[k-1:]))
    return D0,kw,K,[r[0] for r in poss]
for cs in ('0.11','0.037','0.036'):
    D0,kw,K,poss=run(cs)
    print(f'c={cs}: D0={D0} ({float(D0):.4e}), worst k={kw}, possible k={poss}, base K={K}, D0<=1e18: {D0<=10**18}')
X=iv.log(iv.mpf(Pi(15)))
st=iv.sqrt(iv.mpf(53))*(X/(2*X+1))*(iv.log(X)/iv.log(2*X+1))
print('p_16 =',PR[15],' log Pi_15 in',X,' step constant in',st,' > 2*53/52 =',2*53/52,':',st.a>2*53/52)
print('Pi_15 =',Pi(15))
