# hdx_path.py (explore64): "integer-cycle" ceiling. v_x(b) = chi_b(x) are unit vectors of L^2(nu) with <v_x, v_y> = M(x, y); if
# N' in (c', c) is an integer and N' = q_1 ... q_k (prime powers, each an atom), the path x -> x + ln q_1 -> ... -> x + ln N' ends at a point
# no link of diameter <= ln c' shares with x (M = 0), so by the triangle inequality  sum_i sqrt(1 - M(x_{i-1}, x_i)) >= 1  for every
# x in I = [-a, a - ln N'].  With u = alpha_i on I + (ln q_1 + ... + ln q_i), the Kronecker test gives
#   <u, P_M u>/|u|^2 >= [1 / sum_i 1/(2 w_{q_i} alpha_{i-1} alpha_i)] / sum_i alpha_i^2     (min of sum A_i p_i^2 s.t. sum p_i >= 1)
# for EVERY admissible localization M.  Rigorous (arb) at rational alpha.
import sys, math, itertools, numpy as np
from flint import arb, fmpq, ctx
from scipy.optimize import minimize
ctx.prec = 128

def ppow(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n; e = 0
            while m % p == 0: m //= p; e += 1
            return (p, e) if m == 1 else None
    return None

def factorizations(N, kmax=3):   # multisets of prime powers (>1) with product N, size <= kmax
    out = set()
    def rec(n, start, cur):
        if n == 1:
            if cur: out.add(tuple(cur))
            return
        if len(cur) >= kmax: return
        for d in range(start, n + 1):
            if n % d == 0 and ppow(d): rec(n//d, d, cur + [d])
    rec(N, 2, []); return sorted(out)

def bound_float(ws, al):
    k = len(ws); A = [2*ws[i]*al[i]*al[i + 1] for i in range(k)]
    if min(A) <= 0: return 0.0
    return (1/sum(1/x for x in A))/sum(x*x for x in al)

def best(ws):
    k = len(ws); bestv, besta = 0.0, None
    for perm in set(itertools.permutations(ws)):
        f = lambda z: -bound_float(perm, np.exp(z))
        r = minimize(f, np.zeros(k + 1), method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-12, 'maxiter': 4000})
        if -r.fun > bestv: bestv, besta = -r.fun, (perm, np.exp(r.x))
    return bestv, besta

def certify(perm, al, qs):
    den = 10**6; a = [fmpq(int(round(x/max(al)*den)), den) for x in al]
    ws = [arb(ppow(q)[0]).log()/arb(q).sqrt() for q in qs]
    A = [2*ws[i]*arb(a[i])*arb(a[i + 1]) for i in range(len(ws))]
    return (1/sum(1/x for x in A))/sum(arb(x)*arb(x) for x in a)

if __name__ == '__main__':
    pairs = [(fmpq(39), fmpq(227, 6)), (fmpq(40), fmpq(227, 6)), (fmpq(83, 2), fmpq(227, 6)), (fmpq(97, 2), fmpq(47)), (fmpq(50), fmpq(47)),
             (fmpq(72), fmpq(71)), (fmpq(10), fmpq(7943282, 10**6)), (fmpq(47), fmpq(3198, 100)), (fmpq(47), fmpq(44))]
    print('# c, c_link, best integer-cycle bound (N\' = factorization): float, arb lower bound; (prime-power N\' gives the single-atom bound w_N\'/... handled in hdx_ceiling)')
    for cq, cl in pairs:
        rows = []
        for Np in range(int(math.floor(float(cl))) + 1, int(math.ceil(float(cq)))):
            if not (fmpq(Np) > cl and fmpq(Np) < cq): continue
            for fac in factorizations(Np):
                ws = [math.log(ppow(q)[0])/math.sqrt(q) for q in fac]
                v, (perm, al) = best(ws)
                qs = [fac[[i for i in range(len(fac))][0]]]  # placeholder
                # map permuted weights back to prime powers (weights distinct per q up to ties; recompute by matching)
                pool = list(fac); qs = []
                for w in perm:
                    j = min(range(len(pool)), key=lambda t: abs(math.log(ppow(pool[t])[0])/math.sqrt(pool[t]) - w)); qs.append(pool.pop(j))
                cv = certify(perm, al, qs)
                rows.append((v, Np, qs, float(cv.lower())))
        rows.sort(reverse=True)
        top = ' | '.join(f"{Np}={'*'.join(map(str, qs))}: {v:.4f} >= {cv:.4f}" for v, Np, qs, cv in rows[:3])
        print(f'{float(cq):7.3f} {float(cl):8.4f}  {top if rows else "no integer in (c_link, c)"}', flush=True)
