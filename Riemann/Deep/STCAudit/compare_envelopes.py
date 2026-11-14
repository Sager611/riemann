# Where is Casana-Giner's |S| envelope (2.2) below Ismail's nine-row minimum (1.2)? And vs HSW/BW-type bounds.
import math
CG_S = [(0.20176, 0.844), (0.10823, 1.286), (0.09473, 1.946), (0.09349, 2.099), (0.08399, 4.318), (0.08382, 4.444)]
IS = [(0.140750,0.150951,1.377296),(0.134791,0.159606,1.505817),(0.129350,0.169319,1.643495),(0.122384,0.186030,1.859015),(0.118798,0.197172,1.996010),
      (0.115144,0.211179,2.163574),(0.111259,0.154570,2.730936),(0.108109,0.207558,2.795522),(0.105234,0.214408,3.175432)]
cg = lambda L: min(a*L + b for a, b in CG_S)
im = lambda L: min(a*L + b*math.log(L) + c for a, b, c in IS)
worst = (1e9, None); best = (-1e9, None)
L = 1.0
while L <= 200*math.log(10):
    d = im(L) - cg(L)
    worst = min(worst, (d, L)); best = max(best, (d, L)); L += 0.001
print('Ismail_min - CG_env over L in [1, 460.5]: min %.4f at L=%.3f (T=e^L), max %.4f at L=%.3f' % (worst[0], worst[1], best[0], best[1]))
for T in [math.e, 1e3, 3.061e10, 3e12, 1e20, 1e40, 1e100, 1e200]:
    L = math.log(T); print('T=%.3g  CG |S| <= %.4f   Ismail |S| <= %.4f' % (T, cg(L), im(L)))
