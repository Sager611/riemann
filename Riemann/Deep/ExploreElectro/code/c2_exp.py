# c2_exp.py: conditioning experiment for Hankel/Cauer pivots with synthetic moments m_k = sum_j x_j^(k+1), x_j = 1/gamma_j^2.
import sys, time, json
import mpmath as mp
N0 = int(sys.argv[1]) if len(sys.argv) > 1 else 120
mp.mp.dps = 30
t0 = time.time()
G = [mp.zetazero(n).imag for n in range(1, N0 + 1)]
print('zeros', N0, round(time.time() - t0, 1), 's', flush=True)
def pivots(x, nmax, dps):
    with mp.workdps(dps):
        xs = [mp.mpf(1) / mp.mpf(g) ** 2 for g in x]
        m = [mp.fsum(xx ** (k + 1) for xx in xs) for k in range(2 * nmax + 2)]
        out = {}
        for shift in (0, 1):
            A = [[m[i + j + shift] for j in range(nmax)] for i in range(nmax)]
            piv = []
            # LU without pivoting (Gaussian elimination), record pivots
            for k in range(nmax):
                p = A[k][k]
                piv.append(p)
                if p == 0: break
                for i in range(k + 1, nmax):
                    f = A[i][k] / p
                    if f != 0:
                        Ai = A[i]; Ak = A[k]
                        for j in range(k + 1, nmax):
                            Ai[j] -= f * Ak[j]
            out[shift] = piv
        return out
nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 40
res = {}
for dps in (60, 120, 240, 480):
    t0 = time.time(); res[dps] = pivots(G, nmax, dps); print('dps', dps, round(time.time() - t0, 1), 's', flush=True)
ref = res[480]
for dps in (60, 120, 240):
    good = 0
    for n in range(nmax):
        ok = all(abs(res[dps][s][n] / ref[s][n] - 1) < mp.mpf(10) ** -6 for s in (0, 1))
        if ok: good = n + 1
        else: break
    print('dps', dps, 'pivots accurate up to n =', good)
print('log10 pivots (shift 0) n=1..', nmax, [round(float(mp.log10(abs(p))), 1) for p in ref[0][:nmax:5]])
print('signs', all(p > 0 for p in ref[0]) and all(p > 0 for p in ref[1]))
