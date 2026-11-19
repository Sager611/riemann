# Falsification test of Wang lem_norm (2 pi F_I(x) = int |A_I(x,t)|^2 dt) and BGST (2.7) with OFF-LINE synthetic zeros.
# The identity needs the multiset {rho : gamma in I} invariant under rho -> 1 - conj(rho) (functional equation + conjugation).
import mpmath as mp, random
mp.mp.dps = 30
random.seed(1)
def integral27(a):  # int_R du / ((1+u^2)(1+(u+a)^2)) numerically vs 2 pi/(4+a^2)
    f = lambda u: 1/((1+u**2)*(1+(u+a)**2))
    return mp.quad(f, [-mp.inf, -20, 0, 20, mp.inf])
worst = 0
for k in range(6):
    a = mp.mpc(random.uniform(-5,5), random.uniform(-0.99,0.99))
    v = integral27(a); e = 2*mp.pi/(4+a**2)
    worst = max(worst, abs(v-e)/abs(e))
print('(2.7) 6 random a, |Im a|<1: max rel err', mp.nstr(worst,3))
a = mp.mpc(0.3, 1.2)  # outside the range: formula should FAIL
print('(2.7) at Im a = 1.2 (outside range): rel err', mp.nstr(abs(integral27(a)-2*mp.pi/(4+a**2))/abs(2*mp.pi/(4+a**2)),3))
def make_set(sym=True):
    Z = []
    for j in range(5):  # on-line zeros, some double
        g = mp.mpf(random.uniform(0, 12)); m = 2 if j == 0 else 1
        Z += [mp.mpc(0.5, g)]*m
    for j in range(3):  # off-line quadruple partners within I: 1/2+d+ig and 1/2-d+ig
        g = mp.mpf(random.uniform(0, 12)); d = mp.mpf(random.uniform(0.05, 0.45))
        Z += [mp.mpc(0.5+d, g)]
        Z += [mp.mpc(0.5-d, g)] if sym else []
    return Z
def F(Z, x):
    s = 0
    for r in Z:
        for rp in Z:
            z = r - rp
            s += mp.power(x, z) * 4/(4 - z**2)
    return s
def AI2int(Z, x):
    def A(t):
        s = 0
        for r in Z:
            d = r.real - 0.5; g = r.imag
            s += 2*mp.power(x, d + 1j*(g - t))/(1 + ((t-g) + 1j*d)**2)
        return abs(s)**2
    return mp.quad(A, [-mp.inf, -10, 0, 4, 8, 12, 22, mp.inf])
for sym in (True, False):
    Z = make_set(sym)
    for x in (mp.mpf(1), mp.mpf(3.7), mp.mpf(25)):
        lhs = 2*mp.pi*F(Z, x); rhs = AI2int(Z, x)
        print(f'sym={sym} x={mp.nstr(x,3)}: 2piF = {mp.nstr(lhs,12)}  int|A_I|^2 = {mp.nstr(rhs,12)}  rel diff {mp.nstr(abs(lhs-rhs)/abs(rhs),3)}')
