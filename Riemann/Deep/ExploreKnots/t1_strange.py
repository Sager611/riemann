# t1_strange.py — exact check of strange identities (Kontsevich–Zagier F for the trefoil; Hikami's T(2,5) series)
# as formal power series in t (q = e^{-t}); compares with -1/2 e^{shift t/D} sum_k (-1)^k L(-2k-1,chi) (t/D)^k / k!
from fractions import Fraction as Fr
from math import factorial, comb
M = 14

def mul(a, b):
    c = [Fr(0)]*M
    for i, x in enumerate(a):
        if x == 0: continue
        for j in range(M - i):
            if b[j]: c[i+j] += x*b[j]
    return c

def inv_unit(a):
    assert a[0] != 0
    b = [Fr(0)]*M; b[0] = 1/a[0]
    for n in range(1, M):
        b[n] = -sum(a[k]*b[n-k] for k in range(1, n+1))/a[0]
    return b

def one_minus_q(j):  # 1 - e^{-j t}
    return [Fr(0)] + [-Fr((-j)**k, factorial(k)) for k in range(1, M)]

def ratio(a, b):  # (1-q^a)/(1-q^b), both valuation 1
    A = one_minus_q(a)[1:] + [Fr(0)]
    B = one_minus_q(b)[1:] + [Fr(0)]
    # truncation: shifting loses top coefficient; recompute exactly
    A = [-Fr((-a)**(k+1), factorial(k+1)) for k in range(M)]
    B = [-Fr((-b)**(k+1), factorial(k+1)) for k in range(M)]
    return mul(A, inv_unit(B))

def qpow(e):  # q^e = e^{-e t}
    return [Fr((-e)**k, factorial(k)) for k in range(M)]

ONE = [Fr(1)] + [Fr(0)]*(M-1)
poch = [ONE]
for k in range(1, M+1):
    poch.append(mul(poch[-1], one_minus_q(k)))

def qbinom(n, k):
    r = ONE
    for i in range(1, k+1):
        r = mul(r, ratio(n-k+i, i))
    return r

def add(a, b): return [x+y for x, y in zip(a, b)]

# Bernoulli polynomials
def bernoulli_numbers(n):
    B = [Fr(0)]*(n+1); B[0] = Fr(1)
    for m in range(1, n+1):
        B[m] = -sum(comb(m+1, k)*B[k] for k in range(m))/(m+1)
    return B
BN = bernoulli_numbers(2*M+4)
def Bpoly(n, x):
    return sum(comb(n, k)*BN[k]*x**(n-k) for k in range(n+1))

def Lneg(m, chi, f):  # L(-m, chi) for f-periodic chi
    return -Fr(f)**m/(m+1)*sum(chi(r)*Bpoly(m+1, Fr(r, f)) for r in range(1, f+1))

def theta_side(chi, f, D, shift):
    s = [Fr((-1)**k)*Lneg(2*k+1, chi, f)/(factorial(k)*Fr(D)**k) for k in range(M)]
    e = [Fr(shift, D)**k/factorial(k) for k in range(M)]
    return [Fr(-1, 2)*x for x in mul(e, s)]

def chi12(n):
    n %= 12
    return {1: 1, 11: 1, 5: -1, 7: -1}.get(n, 0)
def chi20a(n):  # T(2,5): +1 at +-3, -1 at +-7
    n %= 20
    return {3: 1, 17: 1, 7: -1, 13: -1}.get(n, 0)
def chi20b(n):
    n %= 20
    return {1: 1, 19: 1, 9: -1, 11: -1}.get(n, 0)

# trefoil: Kontsevich–Zagier F = sum_k (q)_k
F1 = [Fr(0)]*M
for k in range(M):
    F1 = add(F1, poch[k])
T1 = theta_side(chi12, 12, 24, 1)
print('trefoil  F(e^-t) coeffs :', [str(x) for x in F1[:7]])
print('trefoil  theta side     :', [str(x) for x in T1[:7]])
print('trefoil match all', M, 'orders:', F1 == T1)

for a, chi, shift in [(1, chi20a, 9), (2, chi20b, 1)]:
    X = [Fr(0)]*M
    for k2 in range(M):
        for k1 in range(k2+1):
            e = k1*k1 + (k1 if a == 1 else 0)
            X = add(X, mul(mul(poch[k2], qpow(e)), qbinom(k2, k1)))
    T = theta_side(chi, 20, 40, shift)
    print(f'T(2,5) a={a} X coeffs   :', [str(x) for x in X[:6]])
    print(f'T(2,5) a={a} theta side :', [str(x) for x in T[:6]])
    print(f'T(2,5) a={a} match:', X == T, ' ratio X/T:', [str(X[i]/T[i]) if T[i] else 'inf' for i in range(5)])

# Glaisher-type integers: T_n = (-1)^n * L(-2n-1,chi12) * ... print first values of L(-2k-1, chi)
print('L(-2k-1,chi12), k=0..7 :', [str(Lneg(2*k+1, chi12, 12)) for k in range(8)])
print('L(-2k-1,chi20a), k=0..7:', [str(Lneg(2*k+1, chi20a, 20)) for k in range(8)])
