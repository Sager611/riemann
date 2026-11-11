# explore30 qmc1: Weil's zero sum as a quadrature with nodes at the zero ordinates (node N1).
# W(f) = sum_rho K(gamma_rho, delta_rho),  K(t, d) = int g(y) cosh(d y) e^{ity} dy,  g = f * f^*,  rho = 1/2 + d + i gamma.
# One zero of an off-line pair (rho, 1 - conj rho) contributes h_d(t) = Re[F(t - id) conj F(t + id)] = K(t, d).
# (1) exactness int h_d = 2 pi ||f||^2 for every d; (2) Poisson: regular off-line lattices with step < pi/a are invisible;
# (3) adversary above T: pairs only in the negative phases of a two-bump f, gaps elsewhere, |S| bounded -> W < 0;
# (4) threshold scan in the half-width a.  Memory: everything chunked; self-abort above 250 MB RSS.
# usage: uv run --no-project --with numpy python3 qmc1_adversary.py [small|full]
import numpy as np, math, sys, time, resource
t0 = time.time()
CAP = 250e6
def memcheck(tag=''):
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss          # bytes on macOS
    if rss > CAP:
        print(f'ABORT memory {rss/1e6:.0f} MB at {tag}', flush=True); sys.exit(3)
    return rss

M = 8                                   # bump phi(s) = (1 - s^2)^M on [-1, 1]
sq, wq = np.polynomial.legendre.leggauss(160)
WP = wq * (1 - sq**2)**M
NORM_PHI2 = float(np.sum(wq * (1 - sq**2)**(2 * M)))
CH = 2048

def Phi(om):                            # int_{-1}^{1} phi(s) e^{i om s} ds, complex om; chunked (~5 MB per chunk)
    om = np.asarray(om, dtype=complex).ravel()
    out = np.empty(om.shape, dtype=complex)
    for i in range(0, len(om), CH):
        out[i:i + CH] = np.exp(1j * np.multiply.outer(om[i:i + CH], sq)) @ WP
    return out

class TwoBump:                          # f = phi((x - a')/w) + e^{i th} phi((x + a')/w), a' = a - w, supp f = [-a, a]
    def __init__(s, a, w, th=0.0):
        s.a, s.w, s.ap, s.th = a, w, a - w, th
        s.norm2 = 2 * w * NORM_PHI2
    def F(s, z):
        z = np.asarray(z, dtype=complex).ravel()
        return s.w * Phi(z * s.w) * (np.exp(1j * z * s.ap) + np.exp(1j * s.th) * np.exp(-1j * z * s.ap))
    def h(s, u, d):                     # one zero of an off-line pair at local height u
        u = np.asarray(u, dtype=float).ravel(); out = np.empty(u.shape)
        for i in range(0, len(u), 20000):
            uu = u[i:i + 20000]
            out[i:i + 20000] = np.real(s.F(uu - 1j * d) * np.conj(s.F(uu + 1j * d)))
        return out
    def h_sum(s, lo, hi, step, d):      # streaming trapezoid sum of h over [lo, hi]
        tot = 0.0; x = lo
        while x < hi:
            uu = np.arange(x, min(x + 20000 * step, hi), step)
            tot += float(np.sum(s.h(uu, d))) * step; x = uu[-1] + step
        return tot

def check_exactness(a=1.925, w=0.2):
    f = TwoBump(a, w); ref = 2 * math.pi * f.norm2
    return [(d, f.h_sum(-400.0, 400.0, 0.005, d), ref) for d in (0.0, 0.1, 0.3734)]

def check_lattice(D, a=1.925, w=0.2, d=0.3734, u0=0.123):
    f = TwoBump(a, w)
    j = np.arange(-int(300 * D), int(300 * D) + 1)
    W = 2 * float(np.sum(f.h(u0 + j * (2.0 / D), d)))
    return W / (D * 2 * math.pi * f.norm2), 2.0 / D, math.pi / a

def adversary(a, w, D, d, q, Lfac=40.0, th=0.0, grid=200, balance=False):
    """pairs only inside the negative intervals of h_d (central fraction q of each), counts by error diffusion so
    that N(u) tracks D*u; on-line lattice of density D outside the block.  Returns W / (D 2pi||f||^2) etc."""
    f = TwoBump(a, w, th)
    L = Lfac / w; P = math.pi / f.ap
    u = np.arange(-L, L, P / grid)
    neg = f.h(u, d) < 0
    memcheck('adversary grid')
    e = np.flatnonzero(np.diff(neg.astype(np.int8)))
    st = [k + 1 for k in e if not neg[k]]; en = [k for k in e if neg[k]]
    if neg[0]: st = [0] + st
    if neg[-1]: en = en + [len(u) - 1]
    ints = [(u[s_], u[e_]) for s_, e_ in zip(st, en) if e_ > s_]
    if len(ints) < 3:
        return None
    left = ints[0][0] - (ints[1][0] - ints[0][1]) / 2
    right = ints[-1][1] + (ints[-1][0] - ints[-2][1]) / 2
    bounds = [left] + [(ints[k][1] + ints[k + 1][0]) / 2 for k in range(len(ints) - 1)] + [right]
    placed = 0; offs = []; s_k = 0.0; I_k = 0.0
    for k, (x, y) in enumerate(ints):
        Pk = bounds[k + 1] - bounds[k]
        if not balance:
            npair = max(0, int(round((D * (bounds[k + 1] - left) - placed) / 2)))
            c = (x + y) / 2; r = (y - x) / 2 * q
            pos = [c] if npair == 1 else (list(c + r * np.linspace(-1, 1, npair)) if npair > 1 else [])
        else:                           # second-order control: keep S(u) = N - D(u - left) and int S near 0
            r = (y - x) / 2 * q; m = (y - x) / 2 - r
            n0 = int(round((D * Pk - s_k) / 2)); best = None
            for npair in (n0 - 1, n0, n0 + 1):
                if npair < 0: continue
                s_next = s_k + 2 * npair - D * Pk
                if abs(s_next) > 2.0: continue
                for sh in np.linspace(-m, m, 21):
                    c = (x + y) / 2 + sh
                    pos = [c] if npair == 1 else (list(c + r * np.linspace(-1, 1, npair)) if npair > 1 else [])
                    I_next = I_k + s_k * Pk + 2 * sum(bounds[k + 1] - pj for pj in pos) - D * Pk**2 / 2
                    key = abs(I_next) + 0.5 * abs(s_next)
                    if best is None or key < best[0]:
                        best = (key, npair, pos, s_next, I_next)
            _, npair, pos, s_k, I_k = best
        placed += 2 * npair
        offs.extend(pos)
    offs = np.array(offs)
    on_r = right + (np.arange(int(2 * L * D) + 1) + 0.5) / D
    on_l = left - (np.arange(int(2 * L * D) + 1) + 0.5) / D
    on_u = np.concatenate([on_l, on_r])
    base = D * 2 * math.pi * f.norm2
    w_on = float(np.sum(np.abs(f.F(on_u))**2)) / base
    w_off = 2 * float(np.sum(f.h(offs, d))) / base              # both members of each pair
    ev = np.repeat(np.sort(offs), 2)
    gu = np.linspace(left, right, 100001)
    S = np.searchsorted(ev, gu, side='right') - D * (gu - left)
    allz = np.sort(np.concatenate([on_l, offs, on_r]))
    Sc = S - S.mean(); S1 = np.cumsum(Sc) * (gu[1] - gu[0])
    memcheck('adversary end')
    return dict(W=w_on + w_off, w_on=w_on, w_off=w_off, swing=float(S.max() - S.min()), S1=float(np.max(np.abs(S1))),
                Sabs=float(max(S.max(), -S.min())), gap=float(np.max(np.diff(allz))), npairs=len(offs),
                frac_neg=float(np.mean(neg)), L=L)

def Sbounds(t):
    tr = 0.112 * math.log(t) + 0.278 * math.log(math.log(t)) + 2.510       # Trudgian 2014
    hsw = 0.1038 * math.log(t) + 0.2573 * math.log(math.log(t)) + 9.3675    # Hasanalizade-Shen-Wong 2022
    return min(tr, hsw), math.log(t / (2 * math.pi)) / (2 * math.pi)

if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'small'
    if mode == 'small':
        r = adversary(1.925, 0.3, 4.28, 0.3734, 0.7, Lfac=15.0, grid=60)
        print('small adversary a=1.925 w=0.3:', {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
        print(f'maxrss {memcheck()/1e6:.0f} MB, {time.time()-t0:.1f}s'); sys.exit(0)
    print('# (1) exactness: int h_d du = 2 pi ||f||^2 (a = 1.925, w = 0.2)')
    for d, I, ref in check_exactness():
        print(f'  d = {d:.4f}: int h = {I:.10f}, 2 pi ||f||^2 = {ref:.10f}, rel. err {abs(I-ref)/ref:.1e}')
    print('# (2) regular lattice of off-line pairs (d = 0.3734): W / (D 2pi||f||^2), exact iff pair step 2/D < pi/a')
    for D in (4.2799, 1.3, 1.0, 0.55):
        x, step, nyq = check_lattice(D)
        print(f'  D = {D}: step {step:.4f} (pi/a = {nyq:.4f}): {x:.12f}')
    print(f'  maxrss {memcheck("lattice")/1e6:.0f} MB', flush=True)
    for label, t in (('T = 3e12', 3e12 + 1e4), ('t = 1e30', 1e30), ('t = 1e100', 1e100)):
        Sb, D = Sbounds(t)
        Sallow = Sb - 1.0                                            # 1.0 reserved for the background S at the block ends
        print(f'# (3) {label}: D = {D:.4f}, min(Trudgian, HSW) |S| <= {Sb:.3f}; adversary allowed |S| <= {Sallow:.3f}, '
              f'allowed gap ~ 2|S|/D = {2*Sallow/D:.3f} (limit 4 pi 0.1038 = {4*math.pi*0.1038:.3f})', flush=True)
        for a in (0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 1.925, 2.131, 2.303, 2.6):
            best = None
            for wf in (0.06, 0.12, 0.2):
                for q in (0.3, 0.7, 1.0):
                    r = adversary(a, wf * a, D, 0.3734, q)
                    if r is None:
                        continue
                    ok = (r['W'] < 0) and (r['Sabs'] <= Sallow)
                    key = (not ok, r['Sabs'] if r['W'] < 0 else 1e9, r['W'])
                    if best is None or key < best[0]:
                        best = (key, wf, q, r)
            _, wf, q, r = best
            ok = r['W'] < 0 and r['Sabs'] <= Sallow
            print(f'  a = {a:5.3f} (c = {math.exp(2*a):8.2f}): best w = {wf}a, q = {q}: W/(D2pi||f||^2) = {r["W"]:+.4f} '
                  f'(on {r["w_on"]:.2e}, off {r["w_off"]:+.4f}), max|S| {r["Sabs"]:.2f}, max gap {r["gap"]:.3f}, '
                  f'pairs {r["npairs"]} -> {"ADMISSIBLE, W<0" if ok else "not admissible"}', flush=True)
        print(f'  maxrss {memcheck(label)/1e6:.0f} MB, {time.time()-t0:.0f}s', flush=True)
