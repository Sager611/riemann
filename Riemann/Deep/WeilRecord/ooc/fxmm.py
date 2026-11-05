"""fxmm.py -- exact fixed-point matrix arithmetic, out of core (weilooc, 2026-11-05).

Number format.  A real x is held as an integer X (the producer's rounding of x 2^F; errors are tracked by the caller)
in offset binary, Y = X + 2^(T-1) in [0, 2^T), as ND = T/16 little-endian uint16 digits.  Range |x| < 2^(T-F-1).

Exact products.  C = round(A B / 2^F) (round half up, |error| <= 2^-F-1 per entry, nothing else) is computed modulo
K primes p < 2^21 with float64 BLAS (every partial sum is an integer < k (p-1)^2 < 2^53 for inner dimension k <= 2048,
so dgemm is exact in any summation order; KMAX = 512 keeps sums < 2^51) and an explicit CRT:
    V = X + Off,  w_i = r_i (P/p_i)^-1 mod p_i,  V = sum_i w_i (P/p_i) - k P,  k = floor(sum_i w_i/p_i),
which is exact because the caller guarantees 0 <= V < P/4 (checked through magnitude bounds), so
sum_i w_i/p_i = k + V/P has fractional part >= 2^-30 > float error 2^-40.  The big constants P/p_i are split
into base-2^16 digits, sum_i w_i mu_is is another exact dgemm (|terms| < 2^37, 45 terms), then carries.
Every routine that returns rounded data returns digits only; the caller adds 2^-F-1 per rounding to its error budget.
"""
import os, math
import numpy as np

def _isprime(n):
    if n < 2: return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0: return n == p
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x*x % n
            if x == n - 1: break
        else: return False
    return True

KMAX = 512           # max inner dimension per dgemm call: every partial sum < KMAX (2^21)^2 = 2^51 (exact), and
                     # the one-step reduction (mod) is exact for |R| < 2^51 (see mod); all runs used inner chunks <= 512

class Fx:
    def __init__(self, F=448, T=480):
        assert F % 16 == 0 and T % 16 == 0 and T > F
        self.F, self.T, self.ND = F, T, T//16
        need = F + T + 6                       # P > 2^(F+T+5) >= 4 (2^(F+T+1) + offsets)
        ps, lg, n = [], 0.0, (1 << 21) - 1
        while lg < need:
            if _isprime(n): ps.append(n); lg += math.log2(n)
            n -= 2
        self.ps = ps; self.K = len(ps)
        P = 1
        for p in ps: P *= p
        self.P = P
        assert P > 1 << (F + T + 5)
        self.p = np.array(ps, dtype=np.float64)
        self.pinv = 1.0/self.p
        self.PW = np.array([[pow(2, 16*s, p) for p in ps] for s in range(self.ND)], dtype=np.float64)   # (ND, K)
        self.OFFIN = np.array([pow(2, T - 1, p) for p in ps], dtype=np.float64)
        # output offset: rounding 2^(F-1), offset binary 2^(F+T-1), CRT guard 2^(F+T)
        self.OFF = (1 << (F - 1)) + (1 << (F + T - 1)) + (1 << (F + T))
        self.OFFOUT = self.res_int(self.OFF)
        self.DP = (P.bit_length() + 15)//16 + 1
        self.CINV = np.array([pow(P//p % p, -1, p) for p in ps], dtype=np.float64)
        self.MU = np.array([[((P//p) >> (16*s)) & 0xFFFF for s in range(self.DP)] for p in ps], dtype=np.float64)   # (K, DP)
        self.PI = np.array([(P >> (16*s)) & 0xFFFF for s in range(self.DP)], dtype=np.float64)
        self.vmax_scale = 1 << (T - F - 1)     # |value| < 2^(T-F-1)

    # ---------- residues ----------
    def res_int(self, x):
        return np.array([x % p for p in self.ps], dtype=np.float64)

    def mod(self, R, g=slice(None)):
        """R (Kg, ...) float64 integers with -2^21 < R < 2^51 + 2^22 -> a representative of R mod p in [0, p], in place
        (primes self.ps[g]).  One step suffices: x = R/p < 2^30.1, the computed R * fl(1/p) differs from x by at most
        x 2^-52 (1 + 2^-50) < 2^-21.9 < 1/p, and frac(x) is 0 or >= 1/p, so q = floor(R * fl(1/p)) equals floor(x) unless
        R = 0 mod p, where the result is 0 or p.  Residues in [0, p] (not [0, p)) are fine everywhere: dgemm bounds use p, and the explicit
        CRT is invariant under w_i -> w_i + p_i (sum w_i P/p_i grows by P and floor(sum w_i/p_i) by 1)."""
        pg = self.p[g]; pig = self.pinv[g]
        sh = (len(pg),) + (1,)*(R.ndim - 1)
        p = pg.reshape(sh); pi = pig.reshape(sh)
        t = R*pi
        np.floor(t, out=t)
        t *= p
        R -= t
        return R

    def digits_f64(self, D):
        """digits (..., ND) -> float64 matrix (n, ND) (exact)."""
        return D.reshape(-1, self.ND).astype(np.float64)

    def res_from_f64(self, Dm, g, shp):
        """residues for primes g of the values whose digits are the rows of Dm (n, ND) -> (Kg,) + shp."""
        R = self.PW[:, g].T @ Dm.T                 # (Kg, n): < 30 2^16 2^21 < 2^42, exact
        R -= self.OFFIN[g][:, None]
        self.mod(R, g)
        return R.reshape((R.shape[0],) + tuple(shp))

    def to_res(self, D):
        """D (..., ND) uint16 digits -> residues (K, ...) of X = Y - 2^(T-1)."""
        shp = D.shape[:-1]
        Dm = D.reshape(-1, self.ND).astype(np.float64)
        R = self.PW.T @ Dm.T                 # (K, n): < 30 2^16 2^21 < 2^42, exact
        del Dm
        R -= self.OFFIN[:, None]
        self.mod(R)
        return R.reshape((self.K,) + shp)

    def matmul_acc(self, acc, A, B, coef=1, g=None):
        """acc[k] = (acc[k] + coef A[k] @ B[k]) mod p_k for the primes k in g (A, B hold those primes' residues);
        A (Kg, r, n), B (Kg, n, c), n <= KMAX, |coef| < 2^20."""
        assert A.shape[2] <= KMAX and B.shape[1] == A.shape[2]
        ks = range(self.K)[g] if g is not None else range(self.K)
        for idx, k in enumerate(ks):
            p = self.ps[k]; pv = self.p[k]; pi = self.pinv[k]
            t = A[idx] @ B[idx]                # exact integers <= KMAX p^2 < 2^52 (residues in [0, p])
            if coef != 1:
                t -= np.floor(t*pi)*pv         # [0, p]
                t *= float(coef % p)           # < 2^42
            a = acc[k]
            a += t                             # < 2^52 + 2^21
            a -= np.floor(a*pi)*pv             # [0, p]

    # ---------- CRT ----------
    def crt_digits(self, R, off):
        """R (K, n): residues of X.  Returns digits (n, DP) of V = X + off, which must satisfy 0 <= V < P/4."""
        n = R.shape[1]
        W = R*self.CINV[:, None]               # R in [0, p]: < 2^42
        W += (self.res_int(off)*self.CINV % self.p)[:, None]
        self.mod(W)                            # w_i in [0, p_i] (representatives; see mod)
        kk = np.floor(self.pinv @ W)           # sum_i w_i/p_i = k + V/P, V/P in [2^-30, 1/4): exact floor
        Y = self.MU.T @ W                      # (DP, n) exact (< 46 2^37)
        del W
        for s in range(self.DP):
            if self.PI[s] != 0: Y[s] -= self.PI[s]*kk
        carry = np.zeros(n)
        q = np.empty(n)
        for s in range(self.DP):
            row = Y[s]
            row += carry
            np.multiply(row, 1.0/65536, out=q); np.floor(q, out=q)
            row -= q*65536
            carry, q = q, carry
        if np.any(carry != 0): raise ArithmeticError("CRT: V outside [0, P)")
        return Y                                # (DP, n) float64 digits in [0, 65536)

    def crt_out(self, R, shape, chunk=32768):
        """R (K, ...) residues of X at scale 2^-2F -> digits of round(X / 2^F) (scale 2^-F), offset binary."""
        Rf = R.reshape(self.K, -1)
        n = Rf.shape[1]
        F16, T16 = self.F//16, self.T//16
        out = np.empty((n, self.ND), dtype=np.uint16)
        for a in range(0, n, chunk):
            b = min(n, a + chunk)
            Y = self.crt_digits(Rf[:, a:b], self.OFF)
            if np.any(Y[F16 + T16] != 1) or np.any(Y[F16 + T16 + 1:] != 0):
                raise OverflowError("fixed-point overflow: result outside the representable range")
            out[a:b] = Y[F16:F16 + T16].T
            del Y
        return out.reshape(tuple(shape) + (self.ND,))

    def crt_float(self, R, shift_bits):
        """R (K, n): residues of a small integer X (|X| < 2^(F+T)).  Returns float64 X / 2^shift_bits and the
        relative error bound of the conversion (applies to |X| as an upper bound: |X|/2^s <= |f| (1 + rel))."""
        off = 1 << (self.F + self.T + 1)
        Rf = R.reshape(self.K, -1)
        if Rf.shape[1] > 32768:
            parts = [self.crt_float(Rf[:, a:a + 32768], shift_bits)[0] for a in range(0, Rf.shape[1], 32768)]
            return np.concatenate(parts), 2.0**-45
        Y = self.crt_digits(Rf, off).T
        top = (self.F + self.T + 1)//16; bit = (self.F + self.T + 1) % 16
        assert bit == 1
        # V = X + 2^(F+T+1): digit 'top' holds 2 (X >= 0) or 1 (X < 0; borrow) ... handle via sign then magnitude
        hi = Y[:, top]
        neg = hi < 2
        mag = Y[:, :top].copy()
        mag[:, :] = np.where(neg[:, None], 65535 - mag, mag)
        # negative: |X| = 2^(16 top) - low = (~low) + 1 over 'top' digits; positive: X = low (+ (hi-2) 2^(16 top))
        carry = neg.astype(np.float64)
        for s in range(top):
            t = mag[:, s] + carry
            c = (t >= 65536).astype(np.float64)
            mag[:, s] = t - 65536*c; carry = c
        extra = np.where(neg, 0.0, hi - 2)
        f = extra.copy()
        for s in range(top - 1, -1, -1):
            f = f*65536.0 + mag[:, s]
        f = np.where(neg, -f, f)
        f = np.ldexp(f, -shift_bits)
        return f, 2.0**-45

    # ---------- digits <-> float, ints ----------
    def to_float(self, D):
        """digits (..., ND) -> float64 values x = X/2^F (relative error <= 2^-48; for diagnostics, solves and norms)."""
        shp = D.shape[:-1]
        Dm = D.reshape(-1, self.ND).astype(np.int64)
        neg = Dm[:, -1] < 0x8000
        mag = Dm.copy()
        mag[:, -1] -= 0x8000                          # X = Y - 2^(T-1): positive -> clear top bit
        # negative: |X| = 2^(T-1) - Y = complement over T-1 bits + 1
        cm = 0xFFFF - Dm
        cm[:, -1] = 0x7FFF - Dm[:, -1]
        carry = np.ones(len(Dm), dtype=np.int64)
        for s in range(self.ND):
            t = cm[:, s] + carry
            carry = t >> 16
            cm[:, s] = t & 0xFFFF
        mag = np.where(neg[:, None], cm, mag)
        f = np.zeros(len(Dm))
        for s in range(self.ND - 1, -1, -1):
            f = f*65536.0 + mag[:, s]
        f = np.where(neg, -f, f)
        return np.ldexp(f, -self.F).reshape(shp)

    def from_float(self, x):
        """float64 array -> digits of X = rint(x 2^F) (|x| < 2^(T-F-1)); exact apart from that one rounding."""
        x = np.asarray(x, dtype=np.float64)
        shp = x.shape
        xf = x.reshape(-1)
        if not np.all(np.abs(xf) < self.vmax_scale): raise OverflowError("from_float: value out of range")
        V = np.rint(np.ldexp(np.abs(xf), self.F))        # exact scaling by 2^F, then round to an integer
        out = np.zeros((len(xf), self.ND), dtype=np.int64)
        for s in range(self.ND - 1, -1, -1):              # V has <= 53 significant bits: every step is exact
            base = math.ldexp(1.0, 16*s)
            dgt = np.floor(V/base)
            V -= dgt*base
            out[:, s] = dgt.astype(np.int64)
        neg = xf < 0
        cm = 0xFFFF - out                                 # two's complement of |X| over T bits
        carry = np.ones(len(xf), dtype=np.int64)
        for s in range(self.ND):
            t = cm[:, s] + carry
            carry = t >> 16
            cm[:, s] = t & 0xFFFF
        res = np.where(neg[:, None], cm, out)
        res[:, -1] ^= 0x8000                              # +2^(T-1) (offset binary) == flip the top bit mod 2^T
        return res.astype(np.uint16).reshape(shp + (self.ND,))

    def int_digits(self, X):
        Y = X + (1 << (self.T - 1))
        if not (0 <= Y < 1 << self.T): raise OverflowError("int_digits: out of range")
        return np.array([(Y >> (16*s)) & 0xFFFF for s in range(self.ND)], dtype=np.uint16)

    def digits_int(self, d):
        Y = 0
        for s in range(self.ND - 1, -1, -1): Y = (Y << 16) | int(d[s])
        return Y - (1 << (self.T - 1))

    def _top64(self, D):
        """signed value of the top 64 bits: X = v 2^(16 (ND-4)) + (low part in [0, 2^(16 (ND-4))))."""
        top = D[..., -4:].astype(np.int64)
        return ((top[..., 3] - 32768)*65536 + top[..., 2])*65536*65536 + top[..., 1]*65536 + top[..., 0]

    def maxabs_bound(self, D):
        """rigorous upper bound of max |x| over a digit block (x = X/2^F): from the top 4 digits."""
        if D.size == 0: return 0.0
        v = self._top64(D)
        # x in [v, v + 1) 2^(16 (ND-4) - F)
        m = max(abs(int(v.min())), abs(int(v.max())) + 1)
        return math.ldexp(float(m) + 1.0, 16*(self.ND - 4) - self.F)*(1 + 2**-50)

    def sumsq_bound(self, D):
        """rigorous upper bound of sum x^2 over a digit block (x = X/2^F), via top 4 digits per entry."""
        if D.size == 0: return 0.0
        v = self._top64(D)
        a = np.maximum(np.abs(v.astype(np.float64)), np.abs(v.astype(np.float64) + 1)) + 2.0   # |X| / 2^(16(ND-4)) <= a
        s = float(np.sum(a*a))*(1 + 1e-9) + 1.0          # float summation: relative error < n 2^-53 < 1e-9
        return math.ldexp(s, 2*(16*(self.ND - 4) - self.F))


class FxMat:
    """Disk-backed fixed-point matrix: memmap uint16 (rows, cols, ND) + JSON sidecar (error bounds, notes)."""
    def __init__(self, fx, path, rows, cols, mode='r'):
        self.fx, self.path, self.rows, self.cols = fx, path, rows, cols
        if mode == 'w':
            self.a = np.lib.format.open_memmap(path, mode='w+', dtype=np.uint16, shape=(rows, cols, fx.ND))
        else:
            self.a = np.load(path, mmap_mode='r' if mode == 'r' else 'r+')
            assert self.a.shape == (rows, cols, fx.ND), (path, self.a.shape)

    @classmethod
    def open(cls, fx, path, mode='r'):
        a = np.load(path, mmap_mode='r')
        return cls(fx, path, a.shape[0], a.shape[1], mode)

    def flush(self):
        if hasattr(self.a, 'flush'): self.a.flush()
