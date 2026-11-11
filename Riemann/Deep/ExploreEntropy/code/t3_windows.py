# t3_windows.py (explore32): Beckner-type window bounds vs exact Galerkin minima vs the prime part, sectors even/odd.
# Parts in CCM's basis (closed forms of weil_matrix.py, copied from Riemann/Deep/WeilSimple): Arch, pole, primes.
import sys, json, numpy as np
from flint import arb, ctx, fmpq
from weil_matrix import WeilMat
ctx.prec = 128
def parts(cq, N, sector):
    s = WeilMat(cq)
    idx = list(range(0, N+1)) if sector == 'even' else list(range(1, N+1))
    n = len(idx); A = np.zeros((n, n)); P = np.zeros((n, n)); R = np.zeros((n, n))
    pv = [float(s.pole_vec(j, sector).mid()) for j in idx]; ps = 2 if sector == 'even' else -2
    Sv = {j: s.S(j) for j in idx}
    for a_, i in enumerate(idx):
        for b_ in range(a_, n):
            j = idx[b_]
            if i == j: ar = s.arch_diag(i, sector)
            else:
                wi, wj = s.w(i), s.w(j); sg = 1 if (i + j) % 2 == 0 else -1
                if sector == 'even': ar = -2*s.Nj(i)*s.Nj(j)*sg*(wj*Sv[j] - wi*Sv[i])/(wi*wi - wj*wj)
                else: ar = -2*(2/s.L)*sg*(wi*Sv[j] - wj*Sv[i])/(wi*wi - wj*wj)
            pr = arb(0)
            for _, ln, wt in s.atoms: pr += wt*s.C(i, j, ln, sector)
            A[a_, b_] = A[b_, a_] = float(ar.mid()); R[a_, b_] = R[b_, a_] = float((-2*pr).mid())
            P[a_, b_] = P[b_, a_] = ps*pv[a_]*pv[b_]
    return A, P, R
def K(u): return sum(np.exp(-(2*m + 0.5)*u)/(2*m + 0.5) for m in range(200000) if (2*m+0.5)*u < 60)
D = -5.372183419225666
lmin = lambda M: float(np.linalg.eigvalsh(M)[0])
out = []
cs = [fmpq(int(x*1000), 1000) for x in [1.2, 1.3, 1.38, 1.45, 1.6, 1.8, 2.0, 2.001, 2.1, 2.2, 2.5, 3.0, 4.0, 5.0]]
N = int(sys.argv[1]) if len(sys.argv) > 1 else 60
for cq in cs:
    c = float(cq); a = np.log(c)/2
    row = {'c': c, 'a': a, 'N': N, 'beckner_global': D, 'killing_even': D + 2*K(a), 'killing_odd': D + 2*K(a) - 2*(np.sinh(a) - a)}
    for sec in ['even', 'odd']:
        A, P, R = parts(cq, N, sec)
        row[sec] = {'arch': lmin(A), 'arch_pole': lmin(A + P), 'primes_min': lmin(R) if np.any(R) else 0.0, 'Q': lmin(A + P + R)}
    out.append(row)
    e, o = row['even'], row['odd']
    print('c=%.3f a=%.4f | kill e %.4f o %.4f | arch e %.4f o %.4f | arch+pole e %.5f o %.5f | primes e %.4f o %.4f | Q e %.3e o %.3e'
          % (c, a, row['killing_even'], row['killing_odd'], e['arch'], o['arch'], e['arch_pole'], o['arch_pole'], e['primes_min'], o['primes_min'], e['Q'], o['Q']), flush=True)
json.dump(out, open('../logs/t3_windows_N%d.json' % N, 'w'), indent=1)
