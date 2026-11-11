# Spira's "RH for sections" (2405.12557 Conj, l. 164): Z_{N(t)}(t) = sum_{k<=N} cos(theta(t) - t ln k)/sqrt(k), N(t) = [t/2]
# Test near close zero pairs of Z: does the section have the two real zeros?
import sys
from mpmath import mp, mpf, siegelz, siegeltheta, cos, log, sqrt, floor, zetazero, nstr
mp.dps = 30
def ZN(t, N):
    th = siegeltheta(t)
    return sum(cos(th - t*log(k))/sqrt(k) for k in range(1, N+1))
def scan(a, b, steps):
    a = mpf(a); b = mpf(b); N = int(floor(a/2)); assert int(floor(b/2)) == N
    prevZ = prevS = None; chZ = chS = 0; minabsZ = (1e9, None); maxS = []
    for i in range(steps+1):
        t = a + (b-a)*i/steps
        z = siegelz(t); s = ZN(t, N)
        if prevZ is not None:
            chZ += (z*prevZ < 0); chS += (s*prevS < 0)
        prevZ, prevS = z, s
        if i % (steps//10) == 0: print('  t=%s Z=%s Z_N=%s' % (nstr(t,10), nstr(z,6), nstr(s,6)))
    print(f'[{a},{b}] N={N}: sign changes Z={chZ}, Z_N={chS}')
for n1 in [int(x) for x in sys.argv[1:]]:
    g1 = zetazero(n1).imag; g2 = zetazero(n1+1).imag
    print('pair', n1, nstr(g1,12), nstr(g2,12))
    a = g1 - (g2-g1); b = g2 + (g2-g1)
    scan(a, b, 200)
