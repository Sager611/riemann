# explore14: F_q sanity test of the single-sign criterion. Coefficient sequences [u^n] of rational functions.
import mpmath as mp
mp.mp.dps = 60
q = mp.mpf(5); sq = mp.sqrt(q)
N = 90
def ser_poly(coeffs):           # polynomial -> series list
    s = [mp.mpf(0)]*(N+1)
    for i,c in enumerate(coeffs):
        if i <= N: s[i] = mp.mpf(c) if not isinstance(c, mp.mpf) else c
    return s
def mul(a,b): return [sum(a[i]*b[n-i] for i in range(n+1)) for n in range(N+1)]
def inv(a):
    b = [mp.mpf(0)]*(N+1); b[0] = 1/a[0]
    for n in range(1,N+1): b[n] = -sum(a[i]*b[n-i] for i in range(1,n+1))/a[0]
    return b
def subs_scale(a, lam): return [a[n]*lam**n for n in range(N+1)]           # f(lam u)
def subs_sq(a): # f(u^2)
    b = [mp.mpf(0)]*(N+1)
    for n in range(0, N//2+1): b[2*n] = a[n]
    return b
def ZC(P):  # Z_C(u) = P(u)/((1-u)(1-qu))
    return mul(ser_poly(P), inv(mul(ser_poly([1,-1]), ser_poly([1,-q]))))
def signs(seq, start=40):
    return ''.join('+' if x>0 else ('-' if x<0 else '0') for x in seq[start:])
def report(name, seq):
    print('%-44s n=40..%d: %s' % (name, N, signs(seq)))
def honest(P):   # G = 1/((1-sqrt(q)u)^2 Z_C(u))
    return mul(inv(mul(ser_poly([1,-sq]), ser_poly([1,-sq]))), inv(ZC(P)))
def literal(P):  # FRS template Z_C(u^2) Z_C(q u^2) / Z_C(u)
    z = ZC(P); z2 = subs_sq(z); zq2 = subs_sq(subs_scale(z, q))
    return mul(mul(z2, zq2), inv(z))
E = [1, 3, 5]                      # y^2=x^3+x+1 over F_5: #E=9, a=-3, P=1+3u+5u^2 (RH: |alpha|^2=5)
fake1 = [1, 5, 5]                  # FE-symmetric, |a|>2sqrt5: real roots off the circle (F_q "DH")
al = 3*mp.e**(1j*1); be = q/al     # genus-2 fake: roots 3e^{±i}, (5/3)e^{±i}: off-circle, complex
def poly_from_roots(rs):
    c = [mp.mpc(1)]
    for r in rs:
        c = [ (c[i] if i < len(c) else 0) - (r*c[i-1] if i>=1 else 0) for i in range(len(c)+1)]
    return [mp.re(x) for x in c]
fake2 = poly_from_roots([al, mp.conj(al), be, mp.conj(be)])
good2 = poly_from_roots([sq*mp.e**(1j*1), sq*mp.e**(-1j*1), sq*mp.e**(2j), sq*mp.e**(-2j)])
print('fake2 coeffs', [mp.nstr(x,6) for x in fake2]); print('good2 coeffs', [mp.nstr(x,6) for x in good2])
report('honest G, E/F_5 (RH true)', honest(E))
report('honest G, genus-2 on-circle (RH true)', honest(good2))
report('honest G, fake1 real off-circle (RH false)', honest(fake1))
report('honest G, fake2 complex off-circle (RH false)', honest(fake2))
report('literal FRS, E/F_5 (RH true)', literal(E))
report('literal FRS, genus-2 on-circle (RH true)', literal(good2))
# Fesenko-type elliptic-surface analog: zeta of a rational elliptic surface Z = 1/((1-u)(1-qu)^rho (1+qu)^m (1-q^2 u))
def surf(rho, m):
    d = mul(ser_poly([1,-1]), ser_poly([1,-q*q]))
    for _ in range(rho): d = mul(d, ser_poly([1,-q]))
    for _ in range(m): d = mul(d, ser_poly([1,q]))
    return inv(d)
def fesenko(rho, m):
    # square, then remove the polar part at u=1/q^2 (abscissa) and u=1 is irrelevant for large n; study coefficients/q^n
    z = surf(rho, m); z2 = mul(z, z)
    # polar part at u = 1/q^2 of order 2: subtract its expansion exactly via partial fractions numerically
    # residue data: z2 = R(u)/(1-q^2 u)^2 ; c2 = R(1/q^2), c1 from derivative
    def R(u): 
        return 1/(((1-u)*(1-q*u)**rho*(1+q*u)**m)**2)
    u0 = 1/(q*q); c2 = R(u0); c1 = mp.diff(R, u0)
    # (1-q^2u)^{-2} = sum (n+1) q^{2n} u^n ; R(u) ≈ c2 + c1 (u-u0) -> c2/(1-q^2u)^2 + c1(u-u0)/(1-q^2u)^2 = c2/(..)^2 - (c1/q^2)/(1-q^2u)
    pol = [ (c2*(n+1) - c1/(q*q))*q**(2*n) for n in range(N+1)]
    return [z2[n]-pol[n] for n in range(N+1)]
report('Fesenko analog, split surface rho=10, m=0', fesenko(10,0))
report('Fesenko analog, twisted surface rho=2, m=8', fesenko(2,8))
report('Fesenko analog, rho=5, m=5', fesenko(5,5))
