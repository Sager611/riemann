# t5_pinning.py (explore39): interlacing of the Weil form Q with its Markov part M = Q - pole (wmat, CCM basis, mpmath).
# Even: Q = M + 2 C C^T  => mu_k <= lam_k <= mu_{k+1};  odd: Q = M - 2 s s^T => lam_k <= mu_k <= lam_{k+1}.
# S(c) = 2 sum_{n<c} Lambda(n)/sqrt(n) - (psi(1/4) - ln pi); E = M + S is the sub-Markov (killed jump) form; lambda_k(E) = mu_k + S.
import sys, math, json, time, mpmath as mp
sys.path.insert(0, '.scratch/explore39'); import wmat
from flint import ctx, fmpq
cs = [float(x) for x in sys.argv[1:]]
for c in cs:
    bits = int(4*math.pi*c*1.4427) + 128; ctx.prec = bits; mp.mp.prec = bits
    N = int(math.ceil(0.6*c*math.log(c))) + 20
    cq = fmpq(int(round(c*64)), 64)
    D = mp.digamma(mp.mpf(1)/4) - mp.log(mp.pi)
    S = -D + 2*sum(mp.log(wmat.ppow(n)[0])/mp.sqrt(n) for n in range(2, int(math.ceil(c))) if n < c and wmat.ppow(n))
    out = {'c': c, 'N': N, 'bits': bits, 'S': mp.nstr(S, 12)}
    t0 = time.time()
    for sec in ['even', 'odd']:
        F = wmat.Form(cq, 'zeta')
        T, idx = F.block(N, sec); lam = sorted(mp.mp.eigsy(wmat.to_mp(T, mp.mp), eigvals_only=True))
        F.pole = False
        T0, _ = F.block(N, sec); mu = sorted(mp.mp.eigsy(wmat.to_mp(T0, mp.mp), eigvals_only=True))
        out[sec] = {'lam': [mp.nstr(v, 6) for v in lam[:4]], 'mu': [mp.nstr(v, 6) for v in mu[:4]],
                    'E0=S+mu0': mp.nstr(S + mu[0], 10), 'E1-S=mu1': mp.nstr(mu[1], 6)}
        if sec == 'even':
            out[sec]['interlace_ok'] = bool(all(mu[k] <= lam[k] + mp.mpf(10)**(-bits//4) and lam[k] <= mu[k+1] + mp.mpf(10)**(-bits//4) for k in range(3)))
        else:
            out[sec]['interlace_ok'] = bool(all(lam[k] <= mu[k] + mp.mpf(10)**(-bits//4) and mu[k] <= lam[k+1] + mp.mpf(10)**(-bits//4) for k in range(3)))
    out['sec'] = round(time.time() - t0, 1)
    print(json.dumps(out), flush=True)
