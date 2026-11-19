# Re-check BGST (2.7) with pole-aware breakpoints (near Im a = +-1 the pole -a +- i approaches the real axis).
import mpmath as mp, random
mp.mp.dps = 30; random.seed(1)
worst = 0
for k in range(12):
    a = mp.mpc(random.uniform(-5,5), random.uniform(-0.99,0.99))
    f = lambda u: 1/((1+u**2)*(1+(u+a)**2))
    c = -a.real
    pts = sorted(set([-mp.inf, c-1, c-0.1, c-0.01, c, c+0.01, c+0.1, c+1, -1, 0, 1, mp.inf]))
    v = mp.quad(f, pts, maxdegree=10); e = 2*mp.pi/(4+a**2)
    r = abs(v-e)/abs(e); worst = max(worst, r)
    print(f'a={mp.nstr(a,6)}: rel err {mp.nstr(r,3)}')
print('max rel err', mp.nstr(worst,3))
