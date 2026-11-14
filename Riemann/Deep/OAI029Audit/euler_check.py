# Brute-force checks of 029 Prop 7.3 (05-poisson.tex l.444-654) and Lemma 3.4 finite-field identities.
# Model: degree-one exterior prime, completion Z_p, p = Q = 1 mod 12; chi = sextic character of F_p^x;
# psi°(y/p^n) = exp(2 pi i y/p^n) (any additive character of conductor Z_p; formulas are sign-agnostic).
import cmath, itertools, sys
PLANT = len(sys.argv) > 1 and sys.argv[1] == 'plant'
def prim_root(p):
    for g in range(2, p):
        if all(pow(g, (p-1)//q, p) != 1 for q in {q for q in range(2, p) if (p-1) % q == 0 and all(q % r for r in range(2, q))}):
            return g
def setup(p):
    g = prim_root(p); lg = {}
    x = 1
    for k in range(p-1): lg[x] = k; x = x*g % p
    z6 = cmath.exp(2j*cmath.pi/6)
    def chi(x, j=1, mask=True):      # chi_p(x)^j with unit mask (0 on nonunits even if 6 | j)
        x %= p
        if x == 0: return 0
        return z6**((lg[x]*j) % 6)
    return chi
def e(x): return cmath.exp(2j*cmath.pi*x)
def C_brute(p, chi, t, k, J):
    tot = 0
    dr = range(p**k) if k > 0 else [0]
    for d in dr:
        if k > 0 and (p**J - p**t*d) % p**k: continue
        cd = chi(d, k) if k > 0 else 1
        if cd == 0: continue
        vr = range(p**t) if t > 0 else [0]
        s = 0
        for v in vr:
            cv = chi(v, t) if t > 0 else 1
            if cv == 0: continue
            s += cv*e(((p**J - p**t*d)*v % p**(k+t))/p**(k+t))
        tot += cd*s
    return tot
def T(p, chi, v, y):
    return sum((chi(x, v) if True else 0)*e((y*x % p)/p) for x in range(1, p))
def C_formula(p, chi, t, k, J):
    Q = p
    tau = lambda v: T(p, chi, v, 1)
    Tv = lambda v, y: T(p, chi, v, y)
    if t == 0: return 1 if k == 0 else (1 if J == 0 else 0)
    if k == 0: return Q**(t-1)*Tv(t, p**(J-t+1)) if J >= t-1 else 0
    if k == 1:
        if J < t: return 0
        tau1 = tau(1) if not PLANT else tau(5)     # planted error: wrong Gauss sum
        return tau1*Q**(t-1)*Tv(t-1, p**(J-t))
    return 0
def main():
    worst = 0; n = 0
    for p in [13, 37]:
        chi = setup(p)
        for t in range(0, 4 if p == 13 else 3):
            for k in range(0, 3):
                if p**(t+k) > 13**4: continue
                for J in range(0, 8):
                    b = C_brute(p, chi, t, k, J); f = C_formula(p, chi, t, k, J)
                    worst = max(worst, abs(b-f)/max(1, abs(b))); n += 1
    print('C_{t,k}(J): %d cases, max rel diff %.2e' % (n, worst))
    # T_v values and tau_0 = -1, |tau_v| = sqrt(Q), tau1*tau5 = Q (chi(-1)=1)
    for p in [13, 37, 61, 73, 97]:
        chi = setup(p)
        G = [T(p, chi, v, 1) for v in range(6)]
        J = lambda a, b: sum(chi(x, a)*chi((1-x) % p, b) for x in range(2, p))
        r1 = abs(G[1]*G[5] - p); r2 = abs(G[0] + 1)
        # Lemma 3.4: (Q^-1/2 G2)^3 = Q^-1/2 J(rho,rho); G1 G2 = conj(chi(4)) J(rho,rho) G3
        l1 = abs((G[2]/p**0.5)**3 - J(2, 2)/p**0.5)
        l2 = abs(G[1]*G[2] - chi(4).conjugate()*J(2, 2)*G[3]) / p**1.5
        sgn = G[3]/p**0.5
        print('p=%d chi(-1)=%s tau1tau5-Q=%.1e tau0+1=%.1e L3.4a=%.1e L3.4b=%.1e gamma3=%.6f%+.6fi' % (p, chi(-1), r1, r2, l1, l2, sgn.real, sgn.imag))

if __name__ == '__main__':
    main()
