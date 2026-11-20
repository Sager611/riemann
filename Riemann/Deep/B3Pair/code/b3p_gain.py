# b3pair: best gain G(k,l) = (1-l-3k)/(1+l-k) over all exponent pairs dominating TTY's beta table (Lemma beta-duality),
# vs A-images (G = (1-s)/(3+s), s = 2k'+l' = 3*Bhat(1/3)). Float LP on a grid, then exact check of the candidates.
from fractions import Fraction as F
import importlib.util, sys
spec = importlib.util.spec_from_file_location('chk', __file__.replace('b3p_gain.py', 'b3p_ttycheck.py'))
src = open(__file__.replace('b3p_gain.py', 'b3p_ttycheck.py')).read().split('def check')[0]
sys.argv = ['x', '10']; exec(src)          # defines T (table rows), cl (classical lines), A, B
rows = T + cl
def Bf(a):  # best tabulated upper bound for beta(alpha), alpha in [0,1]; reflection beta(1-a) = 1/2 - a + beta(a)
    if a > 0.5: return Bf(1 - a) + a - 0.5
    v = [float(r[0]) + float(r[1])*a for r in rows if float(r[2]) <= a < float(r[3]) or (a == 0.5 and r[4] == 'cl')]
    return min(v)
G = [i/20000 for i in range(20001)]
Bv = [Bf(a) for a in G]
best = None
for j in range(0, 20001):
    m = 0.3 + 0.7*j/20000            # slope l - k
    k = max(b - m*a for a, b in zip(G, Bv))
    l = k + m
    g = (1 - l - 3*k)/(1 + l - k)
    if best is None or g > best[0]: best = (g, k, l, m)
print('LP (float, grid 2e4 x 2e4): max G = %.7f at (k,l) = (%.7f, %.7f)' % best[:3])
kA, lA = A((F(89,1282), F(997,1282)))
print('A(89/1282,997/1282) =', (kA, lA), 'G =', (1-lA-3*kA)/(1+lA-kA), '=', float((1-lA-3*kA)/(1+lA-kA)))
kB, lB = A(A((F(13,84), F(55,84))))
print('A^2(13/84,55/84) =', (kB, lB), 'G =', (1-lB-3*kB)/(1+lB-kB), '=', float((1-lB-3*kB)/(1+lB-kB)))
# s* = 3 * concave envelope of B at 1/3 (float): max over a1<1/3<a2 of chord value
s3 = max(((a2-1/3)*Bv[i] + (1/3-a1)*Bv[j2])/(a2-a1) for i, a1 in enumerate(G[:6667:7]) for j2, a2 in [(jj, G[jj]) for jj in range(6667, 20001, 7)])
print('3*Bhat(1/3) (float, chords, step 7) = %.7f ; 1175/1282 = %.7f' % (3*s3, 1175/1282))
