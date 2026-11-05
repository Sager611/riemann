# t1_sieve_twists.py — explore6 (Möbius disjointness -> RH tree), Test A + digit histograms for T1.
# Segmented sieve for mu(n) and lambda(n), n < N = 2^27 (1.34e8). Partial sums of mu(n) f(n) for zero-entropy f,
# dyadic-window maxima, samples, and digit-sum histograms c_k(L) = sum_{n<2^L, s_2(n)=k} mu(n).
# Run: OMP_NUM_THREADS=1 perl -e 'alarm 1800; exec @ARGV' uv run --no-project --with numpy python3 t1_sieve_twists.py
import numpy as np, time, json, sys, math
LOGN = int(sys.argv[1]) if len(sys.argv) > 1 else 27
N = 1 << LOGN
B = 1 << 22 if LOGN >= 22 else N
t0 = time.time()
R = int(math.isqrt(N)) + 1
sv = np.ones(R + 1, bool); sv[:2] = False
for i in range(2, int(R**0.5) + 1):
    if sv[i]: sv[i*i::i] = False
primes = np.nonzero(sv)[0].astype(np.int64)
logp = np.log(primes.astype(np.float64)).astype(np.float32)
TWO64 = 2.0**64
def A64(a):  # alpha -> uint64 fixed point (alpha' = A/2^64)
    return np.uint64(int(round((a % 1.0) * TWO64)) % (1 << 64))
phi = (math.sqrt(5) - 1) / 2
A_phi, A_r2, A_r3 = A64(phi), A64(math.sqrt(2)), A64(math.sqrt(3))
r2 = math.sqrt(2)
def ephase(u64):  # e(u/2^64)
    ang = (u64 >> np.uint64(11)).astype(np.float64) * (2.0 * math.pi / 2.0**53)
    return np.exp(1j * ang)
names = ['mertens', 'liouville', 'rot_golden', 'rot_sqrt2', 'quad_n2sqrt2', 'bracket_heis', 'thue_morse',
         'rudin_shapiro', 'digit_0.1', 'sturmian_golden', 'random_pm1']
carry = {k: 0j for k in names}
dmax = {k: np.zeros(LOGN + 1) for k in names}       # max_{2^k<=y<2^{k+1}} |S(y)|
samp_y = sorted(set([int(round(2 ** (j / 8))) for j in range(0, 8 * LOGN + 1)] + [10**k for k in range(1, 9) if 10**k < N]))
samp_y = [y for y in samp_y if 1 <= y < N]
samples = {k: {} for k in names}
hist = {}            # L -> list c_k
hist_acc = np.zeros(LOGN + 2)
for lo in range(0, N, B):
    hi = lo + B
    n = np.arange(lo, hi, dtype=np.int64)
    mu = np.ones(B, np.int8); lam = np.ones(B, np.int8)
    accm = np.zeros(B, np.float32); accl = np.zeros(B, np.float32)
    for p, lp in zip(primes.tolist(), logp.tolist()):
        if p * p >= hi and p > hi: break
        s = (-lo) % p
        mu[s::p] *= -1; accm[s::p] += lp; lam[s::p] *= -1; accl[s::p] += lp
        pk = p * p
        if pk < hi:
            s2 = (-lo) % pk; mu[s2::pk] = 0
        while pk < hi:
            s2 = (-lo) % pk; lam[s2::pk] *= -1; accl[s2::pk] += lp
            pk *= p
    lnn = np.log(np.maximum(n, 1).astype(np.float64))
    big_m = accm < (lnn - 1.0); mu[big_m] *= -1
    big_l = accl < (lnn - 1.0); lam[big_l] *= -1
    if lo == 0:
        mu[0] = 0; lam[0] = 0
    muf = mu.astype(np.float64)
    nu = n.astype(np.uint64)
    pc = np.bitwise_count(nu).astype(np.int64)
    # twists are built one at a time (memory: one complex block at a time)
    def make_val(k):
        if k == 'mertens': return muf
        if k == 'liouville': return lam.astype(np.float64)
        if k == 'rot_golden': return muf * ephase(nu * A_phi)
        if k == 'rot_sqrt2': return muf * ephase(nu * A_r2)
        if k == 'quad_n2sqrt2': return muf * ephase((nu * nu) * A_r2)
        if k == 'bracket_heis':
            fl = np.floor(n.astype(np.float64) * r2).astype(np.uint64)
            return muf * ephase((fl * nu) * A_r3)
        if k == 'thue_morse': return muf * (1.0 - 2.0 * (pc & 1))
        if k == 'rudin_shapiro':
            rs = (np.bitwise_count(nu & (nu >> np.uint64(1))) & 1).astype(np.float64)
            return muf * (1.0 - 2.0 * rs)
        if k == 'digit_0.1': return muf * np.exp(2j * math.pi * 0.1 * pc)
        if k == 'sturmian_golden':
            st = np.floor((n + 1).astype(np.float64) * phi) - np.floor(n.astype(np.float64) * phi) - phi
            return muf * st
        if k == 'random_pm1':
            rng = np.random.default_rng(1000 + lo // B)
            return muf * (2.0 * rng.integers(0, 2, B) - 1.0)
        raise KeyError(k)
    # digit histograms at 2^L
    if lo == 0:
        for L in range(0, min(LOGN, 22) + 1):
            m = 1 << L
            if m <= B:
                hist[L] = np.bincount(pc[:m], weights=muf[:m], minlength=LOGN + 2)[:LOGN + 2].tolist()
        hist_acc = np.bincount(pc, weights=muf, minlength=LOGN + 2)[:LOGN + 2].astype(np.float64)
    else:
        hist_acc += np.bincount(pc, weights=muf, minlength=LOGN + 2)[:LOGN + 2]
        if (hi & (hi - 1)) == 0:
            hist[int(math.log2(hi))] = hist_acc.tolist()
    del accm, accl, lnn, big_m, big_l
    for k in names:
        cs = np.cumsum(make_val(k)) + carry[k]      # cs[i] = S(n_i) = sum_{m <= lo+i}
        carry[k] = cs[-1]
        a = np.abs(cs)
        for K in range(0, LOGN + 1):
            a0, a1 = max(1 << K, lo), min(1 << (K + 1), hi)
            if a0 < a1:
                dmax[k][K] = max(dmax[k][K], a[a0 - lo:a1 - lo].max())
        for y in samp_y:
            if lo <= y < hi:
                samples[k][y] = complex(cs[y - lo])
        del cs, a
    sys.stderr.write(f"block {lo//B} done {time.time()-t0:.1f}s\n")
# sanity: known M(10^k)
known = {10: -1, 100: 1, 1000: 2, 10000: -23, 100000: -48, 1000000: 212, 10000000: 1037, 100000000: 1928}
chk = {y: (samples['mertens'][y].real, known[y]) for y in known if y in samples['mertens']}
knownL = {10: 0, 100: -2, 1000: -14, 10000: -94, 100000: -288, 1000000: -530, 10000000: -842, 100000000: -3884}
chkL = {y: (samples['liouville'][y].real, knownL[y]) for y in knownL if y in samples['liouville']}
out = {'N': N, 'names': names, 'dmax': {k: dmax[k].tolist() for k in names},
       'samples': {k: {str(y): [v.real, v.imag] for y, v in samples[k].items()} for k in names},
       'hist': {str(L): hist[L] for L in sorted(hist)}, 'check_M': chk, 'check_L': chkL,
       'seconds': time.time() - t0}
json.dump(out, open('t1_sieve_out.json', 'w'))
print('check M(10^k):', chk)
print('check L(10^k):', chkL)
print('seconds', round(time.time() - t0, 1))
