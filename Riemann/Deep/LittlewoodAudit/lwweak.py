# Sensitivity of Theorem 1 to the Dusart input: rerun lwcheck's independent covers with the Mertens error terms weakened to
# the classical Rosser-Schoenfeld shapes, eps_D(L) -> 1/(2L^2) (5x at L=13) and 0.2/L1 + 0.2/L1^2 -> 1/(2 L1) (2.2x at L1=8).
import sys
sys.argv = ['lwcheck.py', 'B']
exec(open('lwcheck.py').read().split('# ------------------------------------------------------------------ main')[0])
def F(ell, P):
    L, h, L1, h1 = P
    k0, k1 = kap(L, h, L1, h1)
    if not (k1.b < 1): return None
    kp = k0/(1 - k1); x = iv.exp(L); x1 = iv.exp(L1)
    Eh = ((1+h)*iv.log(1+h) - h)/h; eD = 1/(2*L**2)
    T0 = iv.exp(-L)/(L*(1 - iv.exp(-2*L))); T1 = iv.exp(-L1)/(1 - iv.exp(-2*L1))
    W = (L1*(1 + h1/2) + EUP + 1/(2*L1) + L1*taubar(L1) + (1+h1)*L1/(2*(x1 - 1)) + (ell - iv.log(PI))/2 - G/2 + T1)
    return C0 - iv.log(L) - Eh - (eD + 1/(x - 1) + taubar(L)/2 + T0) - kp*W
def Ff(ell, P):
    L, h, L1, h1 = P
    if L < 10 or L1 < 8 or not (0 < h <= 2) or not (0 < h1 <= 2): return -1e9
    k0, k1 = kapf(L, h, L1, h1)
    if k1 >= 1: return -1e9
    kp = k0/(1 - k1); Eh = ((1+h)*math.log1p(h) - h)/h; ex, ex1 = math.exp(-L), math.exp(-L1)
    W = (L1*(1 + h1/2) + EUPf + 0.5/L1 + L1*taubarf(L1) + (1+h1)*L1*ex1/(2*(1 - ex1)) + (ell - math.log(math.pi))/2 - Gf/2 + ex1/(1 - ex1*ex1))
    return C0f - math.log(L) - Eh - (0.5/L**2 + ex/(1 - ex) + taubarf(L)/2 + ex/L) - kp*W
ell3 = iv.exp(I(3)); r, P = optimise(float(ell3.a)); v = F(ell3, Pq(r))
print(f'weak: small q F(e^3) - log(0.11/3) = {float((v - iv.log(I(11, 300))).a):+.6f}')
L8 = I(8); om = EUP + 1/(2*L8) + L8*taubar(L8) + 3*L8/(2*(iv.exp(L8) - 1)) - iv.log(PI)/2 - G/2 + iv.exp(-L8)/(1 - iv.exp(-2*L8))
print(f'weak: omega(8)|h=2 <= {float(om.b):.4f} < 0: {om.b < 0}')
for lab, c, lo_, hi_ in (('A 0.11', (11, 100), 3000, 40000), ('B2 0.116', (116, 1000), 4000, 100000)):
    ok, worst, n = cover(c, lo_, hi_, 1000, wmax=256); print(f'weak: cover {lab} [{lo_/1000},{hi_/1000}]: {"PASS" if ok else "FAIL"} {n} cells worst {worst[0] if ok else worst}')
for st in range(3130, 3300, 5):
    ok, worst, n = cover((115, 1000), st, 40000, 1000, wmax=256)
    if ok: print(f'weak: 0.115 holds from loglog q >= {st/1000} (stated 3.13 needs Dusart)'); break
