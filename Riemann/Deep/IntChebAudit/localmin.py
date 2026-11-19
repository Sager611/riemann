import sys
sys.argv = ['indep.py', '../IntCheb/CERTIFICATE_LOWER_new.json', 'none']
exec(open('indep.py').read())
import bisect
mp.mp.dps = 30
x0 = mp.mpf('0.3160865709')
# nearest edges, nearest roots of active polynomials
k = bisect.bisect_left(E, float(x0))
print('edges around:', E[k-2:k+2], ' dist to nearest edge %.3e' % min(abs(e - float(x0)) for e in E))
near = []
for q, c in act:
    for r in roots_mp(q):
        near.append((abs(mp.mpc(r) - x0), mp.nstr(mp.mpc(r), 12), c, len(q) - 1))
near.sort(key=lambda t: t[0])
for t in near[:4]: print('root dist %s root %s c=%.3e deg %d' % (mp.nstr(t[0], 4), t[1], t[2], t[3]))
for d in ['-1e-5', '-3e-6', '-1e-6', '-3e-7', '-1e-7', '0', '1e-7', '3e-7', '1e-6', '3e-6', '1e-5']:
    print('A(x0%s) = %s' % (d if d.startswith('-') else '+' + d, mp.nstr(A_mp(x0 + mp.mpf(d)), 13)))
gr = (mp.sqrt(5) - 1) / 2
a, b = x0 - mp.mpf('1e-7'), x0 + mp.mpf('1e-7')
c1 = b - gr * (b - a); c2 = a + gr * (b - a); f1 = A_mp(c1); f2 = A_mp(c2)
for it in range(60):
    if f1 < f2: b, c2, f2 = c2, c1, f1; c1 = b - gr * (b - a); f1 = A_mp(c1)
    else: a, c1, f1 = c1, c2, f2; c2 = a + gr * (b - a); f2 = A_mp(c2)
xm = (a + b) / 2; fm = A_mp(xm)
print('LOCAL MIN A = %s at x = %s ; V = %r ; A_min - V = %s ; mirror A(1-x) = %s' % (mp.nstr(fm, 15), mp.nstr(xm, 15), V, mp.nstr(fm - V, 6), mp.nstr(A_mp(1 - xm), 15)))
