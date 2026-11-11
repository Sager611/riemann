# Jutila fallback: thresholds lambda*(theta) for theta = 1/2 (density slope 1-eps, Jutila 1983) and theta = 1/8 (Selberg 1946):
# (a) Lamzouri's Selberg-minorant family, closed form (exact up to float); (b) float SDP of ../SmallGap/sdp.py (uncertified).
import sys, math, time
sys.path.insert(0, '../SmallGap')
def J(l, th):
    return l + l*th**2 - 2/3*l*l*th**3 - th*math.cos(2*math.pi*l*th)/(2*math.pi**2) + math.sin(2*math.pi*l*th)/(4*math.pi**3*l)
def thr(th):
    lo, hi = 0.3, 1.0
    for _ in range(80):
        mid = (lo + hi)/2
        lo, hi = (lo, mid) if J(mid, th) > 1 else (mid, hi)
    return hi
for th in (1/8, 1/2, 4/7 - 1e-5):
    print(f'Selberg minorant: theta={th:.6f} lambda*={thr(th):.6f}', flush=True)
if len(sys.argv) > 1:
    from sdp import solve
    d = int(sys.argv[1]); th = 0.5
    best = (9, None)
    for s in (1.0, 1.1, 1.2, 1.3):
        t = time.time(); lo, hi = 0.79, 0.83
        while hi - lo > 2e-5:
            mid = (lo + hi)/2
            v = solve(d, s, mid, theta=th)[0]
            lo, hi = (lo, mid) if (v == v and v > 1) else (mid, hi)
        print(f'SDP d={d} s={s} theta={th}: lambda* <= {hi:.5f}  [{time.time()-t:.0f}s]', flush=True)
        best = min(best, (hi, s))
    print(f'SDP best d={d}: {best}', flush=True)
