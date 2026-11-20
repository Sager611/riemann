# ew_ritz.py (explore67): rigorous Rayleigh-Ritz UPPER bounds (Arb) for low eigenvalues of Weil's form Q and of its Markov part
# B = Q - pole on trial vectors in E_N (exact dyadic coefficients; Galerkin entries are closed-form Arb enclosures, weil_matrix.py).
# Prints upper bounds: QoddU >= E1odd(Q); Bev2U >= e2ev(B) (2-dim Ritz on B's two lowest Galerkin vectors); BoddU >= e1odd(B);
# Qev2U >= E2ev(Q) (2-dim Ritz). Usage: ew_ritz.py N bits c1 c2 ...
import sys, os, json, time
from flint import arb, acb, arb_mat, acb_mat, ctx, fmpq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from weil_matrix import WeilMat
def parse_c(s):
    if '/' in s: p, q = s.split('/'); return fmpq(int(p), int(q))
    return fmpq(int(s))
def vecs(T, k):   # k lowest approximate eigenvectors of mid(T) by shifted inverse iteration + Gram-Schmidt, exact dyadics
    n = T.nrows(); ev = sorted([e.real.mid() for e in acb_mat(T).eig(algorithm='approx')])
    out = []
    for i in range(k):
        gap = min([abs(ev[i] - ev[j]) for j in range(n) if j != i])
        sh = ev[i] - gap/1000
        Ts = arb_mat(T)
        for a in range(n): Ts[a, a] = Ts[a, a] - sh
        v = arb_mat(n, 1)
        for a in range(n): v[a, 0] = arb(1 + ((7*a + 3*i) % 5))/10
        for it in range(8):
            y = Ts.solve(v, algorithm='approx')
            nr = sum((y[a, 0]*y[a, 0] for a in range(n)), arb(0)).sqrt()
            for a in range(n): v[a, 0] = (y[a, 0]/nr).mid()
        w = [v[a, 0] for a in range(n)]
        for u in out:
            d = sum((w[a]*u[a] for a in range(n)), arb(0)); w = [(w[a] - d*u[a]).mid() for a in range(n)]
        nr = sum((x*x for x in w), arb(0)).sqrt(); w = [(x/nr).mid() for x in w]
        out.append(w)
    return out, ev
def form(T, u, v):
    n = len(u); s = arb(0)
    for a in range(n):
        t = arb(0)
        for b in range(n): t += T[a, b]*v[b]
        s += u[a]*t
    return s
def ritz2_upper(T, u, v):   # larger root of det(A - mu M) = 0, A = [u,v]^T T [u,v], M = Gram; rigorous enclosure
    A11, A12, A22 = form(T, u, u), form(T, u, v), form(T, v, v)
    dot = lambda x, y: sum((x[a]*y[a] for a in range(len(x))), arb(0))
    M11, M12, M22 = dot(u, u), dot(u, v), dot(v, v)
    a = M11*M22 - M12*M12; b = -(A11*M22 + A22*M11 - 2*A12*M12); c = A11*A22 - A12*A12
    return (-b + (b*b - 4*a*c).sqrt())/(2*a)
N = int(sys.argv[1]); ctx.prec = int(sys.argv[2])
for cs in sys.argv[3:]:
    t0 = time.time(); cq = parse_c(cs); W = WeilMat(cq); res = {'c': cs, 'N': N, 'prec': ctx.prec}
    for sec, ps in (('even', 2), ('odd', -2)):
        T, idx = W.block(N, sec); n = len(idx); pv = [W.pole_vec(j, sec) for j in idx]
        B = arb_mat(T)
        for a in range(n):
            for b in range(n): B[a, b] = T[a, b] - ps*pv[a]*pv[b]
        if sec == 'even':
            (v1, v2), evq = vecs(T, 2); res['Qev2U'] = arb(ritz2_upper(T, v1, v2).upper()).str(8, radius=False)
            res['Qev1U'] = arb((form(T, v1, v1)/sum((x*x for x in v1), arb(0))).upper()).str(8, radius=False)
            (w1, w2), evb = vecs(B, 2); res['Bev2U'] = arb(ritz2_upper(B, w1, w2).upper()).str(8, radius=False)
            res['Bev_galerkin'] = [x.str(6, radius=False) for x in evb[:3]]
        else:
            (v1,), evq = vecs(T, 1); res['QoddU'] = arb((form(T, v1, v1)/sum((x*x for x in v1), arb(0))).upper()).str(8, radius=False)
            (w1,), evb = vecs(B, 1); res['BoddU'] = arb((form(B, w1, w1)/sum((x*x for x in w1), arb(0))).upper()).str(8, radius=False)
            res['Bodd_galerkin'] = [x.str(6, radius=False) for x in evb[:2]]
    res['sec'] = round(time.time() - t0, 1)
    print(json.dumps(res), flush=True)
