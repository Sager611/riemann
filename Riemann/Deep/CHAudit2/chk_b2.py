# chaudit2: exact square-free falsification tests for Sec. 8 of 2511.14736v1 (Lemmas difintsq, sugr, trapeze, oldtrivial,
# Prop. andalas, MOST's |R(x)| <= 1.12543 x^(1/4)), with an exact mu-sieve to N = 1e7 (x <= 5e13).
import math, random, sys
import numpy as np
def P(*a):
    print(" ".join(str(t) for t in a)); sys.stdout.flush()
N = 10**7
mu = np.ones(N + 1, dtype=np.int8); mu[0] = 0
isp = np.ones(N + 1, dtype=bool); isp[:2] = False
for p in range(2, int(N**0.5) + 1):
    if isp[p]:
        isp[p*p::p] = False
primes = np.nonzero(isp)[0]
del isp
for p in primes:
    p = int(p)
    mu[p::p] *= -1
    if p*p <= N:
        mu[p*p::p*p] = 0
del primes
M = np.cumsum(mu, dtype=np.int32)
Q = np.cumsum(mu != 0, dtype=np.int32)
z2 = math.pi**2/6
P("sieve done; M(1e7) =", int(M[N]), " Q(1e7) =", int(Q[N]), "(known: M(1e7) = 1928, Q(1e7) = 6079291)")

def Qx(x):  # exact Q(x) for integer x <= N^2
    r = math.isqrt(x); s = 0
    for lo in range(1, r + 1, 10**6):
        d = np.arange(lo, min(r, lo + 10**6 - 1) + 1, dtype=np.int64)
        s += int(np.sum(mu[lo:lo + d.size].astype(np.int64) * (x // (d*d))))
    return s
def Mv(t):  # M(t) for real t >= 0
    return int(M[int(math.floor(t))]) if t >= 1 else 0

# (0) identity Q(x) = sum_k M(sqrt(x/k)) (Lemma difintsq, first line)
random.seed(3); bad = 0
for _ in range(200):
    x = random.randint(1, 10**6)
    k = np.arange(1, x + 1, dtype=np.float64)
    s = int(np.sum(M[np.floor(np.sqrt(x / k)).astype(np.int64)]))
    if s != int(Q[x]): bad += 1
P("(0) Q(x) == sum_{k<=x} M(sqrt(x/k)) at 200 random x <= 1e6: failures =", bad)

# (i) Lemma oldtrivial |R(x)| <= sqrt x and MOST |R(x)| <= 1.12543 x^(1/4), all real x in [1, 1e7] (left limits included)
b12 = (0, 0); b14 = (0, 0); big = []
for lo in range(1, N + 1, 10**6):
    hi = min(N, lo + 10**6 - 1)
    n = np.arange(lo, hi + 1, dtype=np.float64)
    Rn = Q[lo:hi + 1].astype(np.float64)
    amax = np.maximum(np.abs(Rn - n/z2), np.abs(Rn - (n + 1)/z2))
    r12 = amax/np.sqrt(n); r14 = amax/n**0.25
    i = int(np.argmax(r12)); j = int(np.argmax(r14))
    if r12[i] > b12[0]: b12 = (float(r12[i]), lo + i)
    if r14[j] > b14[0]: b14 = (float(r14[j]), lo + j)
    big += (np.nonzero(r12 > 1)[0] + lo).tolist()[:20]
P("(i) sup_{1<=x<1e7+1} |R(x)|/sqrt(x) <=", b12, "; sup |R(x)|/x^(1/4) <=", b14, "(upper estimates: max of |R| at n and (n+1)^-, divided by the value at n)")
P("    x with |R|/sqrt x > 1 (upper estimate):", big[:20])
# ramando, exactly as stated, N = 2..1e6
NN = np.arange(1, 10**6 + 1, dtype=np.float64)
S2 = np.cumsum(1/np.sqrt(NN)); S4 = np.cumsum(NN**-0.25)
lhs1 = 2*np.sqrt(NN + 0.5) + S2; rhs1 = 4*np.sqrt(NN - 0.5)
lhs2 = 4/3*(NN + 0.5)**0.75 + S4; rhs2 = 8/3*(NN - 0.5)**0.75
P("ramando: max over 2<=N<=1e6 of (2 sqrt(N+1/2) + sum n^-1/2) - 4 sqrt(N-1/2) =", float(np.max((lhs1 - rhs1)[1:])), "(<= 0 needed)")
P("ramando: max over 3<=N<=1e6 of LHS2 - RHS2 =", float(np.max((lhs2 - rhs2)[2:])), "(<= 5/7 = 0.714286 needed); max over 37<=N<=1e6:", float(np.max((lhs2 - rhs2)[36:])), "(<= 0 needed)")
P("ramando at N = 36 (not claimed):", float((lhs2 - rhs2)[35]))
del NN, S2, S4, lhs1, rhs1, lhs2, rhs2
# samples above 1e7 (exact Q via mu-sum), integer x and x^-
mx = (0, 0)
for _ in range(120):
    x = int(10**random.uniform(7, math.log10(5e13)))
    qx = Qx(x); sq = int(Qx(x) - Qx(x - 1))
    for val in (abs(qx - x/z2), abs(qx - sq - x/z2)):
        if val/x**0.25 > mx[0]: mx = (val/x**0.25, x)
P("    120 random x in [1e7, 5e13]: max |R(x)|/x^(1/4) (x and x^-) =", mx)

# (iii) per-k tests of Lemma trapeze (with the edesmo pairs) and Lemma sugr
pairs = [(1, 1), (3/4, 3/2), (2/3, 8/3), (16/25, 114/25), (768/1225, 9458/1225), (18432/29645, 361192/29645), (442368/715715, 14328304/715715)]
def pieces(x, a, b):
    tlo, thi = math.sqrt(x/b), math.sqrt(x/a)
    ns = np.arange(math.floor(tlo) + 1, math.floor(thi) + 1, dtype=np.float64)
    bp = x/ns**2
    bp = bp[(bp > a) & (bp < b)]
    pts = np.sort(np.concatenate(([a], bp, [b])))
    mids = 0.5*(pts[1:] + pts[:-1]); lens = np.diff(pts)
    vals = M[np.floor(np.sqrt(x/mids)).astype(np.int64)].astype(np.float64)
    return lens, vals, tlo, thi
worstT = {p: (0,) for p in pairs}; worstS = (0,); nT = 0
xs = [int(10**random.uniform(8, 13)) for _ in range(14)] + [10**8, 10**10, 10**12, 10**13]
for x in xs:
    kmax = min(4000, x // 2)
    for k in range(1, kmax + 1):
        lens, vals, tlo, thi = pieces(float(x), k - 0.5, k + 0.5)
        mk = Mv(math.sqrt(x/k))
        J = float(np.sum(lens*np.abs(vals - mk))); S = abs(float(np.sum(lens*(vals - mk))))
        cnt = int(Q[math.ceil(thi) - 1] - Q[math.floor(tlo)])
        nT += 1
        if cnt == 0:
            if S > 1e-9: worstS = (float('inf'), x, k)
        elif S/(cnt/2) > worstS[0]: worstS = (S/(cnt/2), x, k, cnt)
        for c1, c2 in pairs:
            rhs = c1/4*math.sqrt(x)*(1/math.sqrt(k - 0.5) - 1/math.sqrt(k + 0.5)) + c2
            if J/rhs > worstT[(c1, c2)][0]: worstT[(c1, c2)] = (J/rhs, x, k)
P("(iii) per-k tests on", nT, "(x, k) pairs, x in [1e8, 1e13], k <= 4000:")
P("    sugr: max |int (M(sqrt(x/u)) - M(sqrt(x/k))) du| / (#sqfree/2) =", worstS)
for p in pairs: P("    trapeze pair", p, ": max LHS/RHS =", worstT[p])

# (iv) Prop. andalas end to end: D(x) = R(x) - [sum_{k<=K} M(sqrt(x/k)) - int_0^{K+1/2} M(sqrt(x/u)) du]
def Dx(x, K):
    U = K + 0.5
    Nn = math.isqrt(int(x // U)) if U <= x else 0
    while (Nn + 1)**2 * U <= x: Nn += 1
    while Nn > 0 and Nn**2 * U > x: Nn -= 1
    s_int = 0; s_frac = 0.0
    for lo in range(1, Nn + 1, 10**6):
        nn = np.arange(lo, min(Nn, lo + 10**6 - 1) + 1, dtype=np.int64)
        m = mu[lo:lo + nn.size].astype(np.int64); n2 = nn*nn
        s_int += int(np.sum(m * (x // n2)))
        s_frac += float(np.sum(m * ((x % n2).astype(np.float64) / n2.astype(np.float64))))
    sumM = sum(Mv(math.sqrt(x/k)) for k in range(1, K + 1))
    return Qx(x) - sumM + U*Mv(math.sqrt(x/U)) - s_int - s_frac
worstA = (0,); nA = 0
for x in xs[:10] + [10**8, 10**10, 10**12]:
    c3 = {}
    for c1, c2 in pairs[3:5]:
        cc3 = 0.5*(1/(2*z2) - c1/4)
        ups0 = (cc3*math.sqrt(x)/c2)**(2/3); Kp0 = math.ceil(ups0 - 0.5)
        y0 = x**0.2/(16*0.570591/c1)**0.8; K0 = max(0, math.ceil(y0 - 0.5))
        trials = [(K0, Kp0)] + [(K, K + random.randint(0, 3*Kp0 + 5)) for K in [random.randint(0, 60) for _ in range(5)]]
        for K, Kp in trials:
            if Kp < K: Kp = K
            D = Dx(x, K)
            E = c1*math.sqrt(x)/4*(1/math.sqrt(K + 0.5) - 1/math.sqrt(Kp + 0.5)) + c2*(Kp - K) + 0.5*int(Q[int(math.floor(math.sqrt(x/(Kp + 0.5))))])
            nA += 1
            if abs(D)/E > worstA[0]: worstA = (abs(D)/E, x, K, Kp, (c1, c2), D, E)
P("(iv) andalas end to end on", nA, "(x, K, K', pair) cases: max |D|/bound =", worstA)
