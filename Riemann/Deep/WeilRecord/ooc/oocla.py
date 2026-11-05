"""oocla.py -- out-of-core block linear algebra on fixed-point matrices (weilooc, 2026-11-05).

Disk matrices are .npy files (uint16, shape (rows, cols, ND)) or raw .bin files with a known shape, accessed by explicit
seek/read/write of row segments (no memmap: the resident memory of a process is its working blocks only).
A View is a rectangular window of a disk matrix (or of an in-memory array), possibly transposed; Stack stacks views.
combine_block() evaluates, for one output block, an exact integer linear combination of
  ('prod', A, B, c)  : c * A B                    (A, B views, c small integer)
  ('mat',  M, c)     : (c / 2^F) * M              (c integer: a dyadic scalar times 2^F)
  ('diag', c)        : (c / 2^(2F)) * I           (c integer: a dyadic scalar times 2^(2F))
in residue arithmetic (primes in groups, to bound memory) and rounds ONCE to the output (|error| <= 2^-F-1 per entry).
A rigorous bound of the exact result (from the top digits of every input block) is checked against the representable
range before the CRT, so the CRT precondition 0 <= V < P/4 holds (fxmm.py).  Outputs are written block by block with a
progress file: an interrupted run resumes at the first unfinished block.
"""
import os, json, time, math
import numpy as np
from fxmm import Fx, KMAX

KC = min(KMAX, int(os.environ.get('OOC_KC', '384')))     # inner chunk
NG = int(os.environ.get('OOC_NGROUPS', '3'))             # prime groups (memory: residues of one group at a time)

class FileMat:
    """uint16 matrix (rows, cols, ND) in a .npy file (or raw file with shape given), row segments read on demand."""
    def __init__(self, path, mode='r', shape=None):
        self.path = path
        if shape is None:
            with open(path, 'rb') as f:
                version = np.lib.format.read_magic(f)
                rd = np.lib.format.read_array_header_1_0 if version == (1, 0) else np.lib.format.read_array_header_2_0
                shape, fortran, dtype = rd(f)
                assert not fortran and dtype == np.dtype('<u2'), (path, dtype, fortran)
                self.off = f.tell()
        else:
            self.off = 0
        self.shape = tuple(shape); self.rows, self.cols, self.nd = self.shape
        self.f = open(path, 'rb' if mode == 'r' else 'r+b')
    @classmethod
    def create(cls, path, rows, cols, nd, fill=None):
        with open(path, 'wb') as f:
            np.lib.format.write_array_header_1_0(f, {'descr': '<u2', 'fortran_order': False, 'shape': (rows, cols, nd)})
            off = f.tell()
            if fill is None:
                f.truncate(off + rows*cols*nd*2)
            else:
                row = np.tile(np.asarray(fill, dtype=np.uint16), (cols, 1))
                for _ in range(rows): f.write(row.tobytes())
        return cls(path, 'w')
    def read(self, r0, r1, c0, c1):
        out = np.empty((r1 - r0, c1 - c0, self.nd), dtype=np.uint16)
        if r1 <= r0 or c1 <= c0: return out
        rb = self.cols*self.nd*2
        if c0 == 0 and c1 == self.cols:
            self.f.seek(self.off + r0*rb)
            if self.f.readinto(memoryview(out).cast('B')) != out.nbytes: raise IOError(f"short read {self.path}")
        else:
            seg = (c1 - c0)*self.nd*2
            for i in range(r0, r1):
                self.f.seek(self.off + i*rb + c0*self.nd*2)
                if self.f.readinto(memoryview(out[i - r0]).cast('B')) != seg: raise IOError(f"short read {self.path}")
        return out
    def write(self, r0, c0, blk):
        blk = np.ascontiguousarray(blk, dtype=np.uint16)
        rb = self.cols*self.nd*2
        for i in range(blk.shape[0]):
            self.f.seek(self.off + (r0 + i)*rb + c0*self.nd*2); self.f.write(blk[i].tobytes())
    def flush(self):
        self.f.flush(); os.fsync(self.f.fileno())
    def close(self):
        self.f.close()

class View:
    """window [r0, r0 + rows) x [c0, c0 + cols) of a FileMat or an in-memory (rows, cols, ND) array; .T transposes."""
    def __init__(self, src, r0=0, c0=0, rows=None, cols=None, trans=False):
        self.src, self.r0, self.c0 = src, r0, c0
        R, C = src.shape[0], src.shape[1]
        self.nr = R - r0 if rows is None else rows; self.nc = C - c0 if cols is None else cols
        self.trans = trans
        self.shape = (self.nc, self.nr) if trans else (self.nr, self.nc)
    def _raw(self, a0, a1, b0, b1):
        if isinstance(self.src, FileMat): return self.src.read(self.r0 + a0, self.r0 + a1, self.c0 + b0, self.c0 + b1)
        return np.ascontiguousarray(self.src[self.r0 + a0:self.r0 + a1, self.c0 + b0:self.c0 + b1])
    def block(self, r0, r1, c0, c1):
        if self.trans: return np.ascontiguousarray(np.transpose(self._raw(c0, c1, r0, r1), (1, 0, 2)))
        return self._raw(r0, r1, c0, c1)
    def sub(self, r0, r1, c0, c1):
        """sub-view in logical coordinates."""
        if self.trans: return View(self.src, self.r0 + c0, self.c0 + r0, c1 - c0, r1 - r0, True)
        return View(self.src, self.r0 + r0, self.c0 + c0, r1 - r0, c1 - c0, False)
    @property
    def T(self):
        return View(self.src, self.r0, self.c0, self.nr, self.nc, not self.trans)

class Stack:
    """Views stacked along axis 0 (rows) or axis 1 (columns) of the logical matrix."""
    def __init__(self, srcs, axis=0):
        self.srcs, self.axis = srcs, axis
        if axis == 0:
            self.offs = np.cumsum([0] + [s.shape[0] for s in srcs]); self.shape = (int(self.offs[-1]), srcs[0].shape[1])
        else:
            self.offs = np.cumsum([0] + [s.shape[1] for s in srcs]); self.shape = (srcs[0].shape[0], int(self.offs[-1]))
    def block(self, r0, r1, c0, c1):
        parts = []
        lo, hi = (r0, r1) if self.axis == 0 else (c0, c1)
        for s, a, b in zip(self.srcs, self.offs[:-1], self.offs[1:]):
            x0, x1 = max(lo, int(a)), min(hi, int(b))
            if x0 >= x1: continue
            if self.axis == 0: parts.append(s.block(x0 - int(a), x1 - int(a), c0, c1))
            else: parts.append(s.block(r0, r1, x0 - int(a), x1 - int(a)))
        return np.concatenate(parts, axis=self.axis) if len(parts) > 1 else parts[0]
    @property
    def T(self):
        return Stack([s.T for s in self.srcs], 1 - self.axis)

def inner_splits(srcA):
    """inner-dimension chunk boundaries: <= KC, aligned to the pieces of a column-stacked A."""
    n = srcA.shape[1]
    cuts = [0]
    if isinstance(srcA, Stack) and srcA.axis == 1:
        for a, b in zip(srcA.offs[:-1], srcA.offs[1:]):
            x = int(a)
            while x < b:
                y = min(int(b), x + KC); cuts.append(y); x = y
        return cuts
    x = 0
    while x < n:
        x = min(n, x + KC); cuts.append(x)
    return cuts

def groups(K):
    g = max(1, NG); step = -(-K//g)
    return [slice(i, min(K, i + step)) for i in range(0, K, step)]

def combine_block(fx, r0, r1, c0, c1, terms, stats=None, raw=False):
    K = fx.K; F = fx.F
    acc = np.zeros((K, r1 - r0, c1 - c0))
    G = groups(K)
    bound = 0.0              # rigorous bound of |exact result| (value scale)
    for term in terms:
        kind = term[0]
        if kind == 'prod':
            _, A, B, coef = term
            assert A.shape[1] == B.shape[0], (A.shape, B.shape)
            cuts = inner_splits(A)
            for k0, k1 in zip(cuts[:-1], cuts[1:]):
                DA = A.block(r0, r1, k0, k1); DB = B.block(k0, k1, c0, c1)
                bound += abs(coef)*(k1 - k0)*fx.maxabs_bound(DA)*fx.maxabs_bound(DB)
                FA = fx.digits_f64(DA); FB = fx.digits_f64(DB); del DA, DB
                for g in G:
                    RA = fx.res_from_f64(FA, g, (r1 - r0, k1 - k0))
                    RB = fx.res_from_f64(FB, g, (k1 - k0, c1 - c0))
                    fx.matmul_acc(acc, RA, RB, coef, g)
                    del RA, RB
                del FA, FB
        elif kind == 'mat':
            _, M, c = term
            if c == 0: continue
            DM = M.block(r0, r1, c0, c1)
            bound += abs(c)/2.0**F*fx.maxabs_bound(DM)*(1 + 1e-12)
            FM = fx.digits_f64(DM); del DM
            for g in G:
                R = fx.res_from_f64(FM, g, (r1 - r0, c1 - c0))
                R *= fx.res_int(c)[g][:, None, None]
                fx.mod(R, g)
                acc[g] += R; del R
                fx.mod(acc[g], g)
            del FM
        elif kind == 'diag':
            _, c = term
            lo, hi = max(r0, c0), min(r1, c1)
            if lo < hi and c != 0:
                bound += abs(c)/2.0**(2*F)*(1 + 1e-12)
                rc = fx.res_int(c)
                for gg in range(lo, hi):
                    acc[:, gg - r0, gg - c0] += rc
                fx.mod(acc)
        else:
            raise ValueError(kind)
    if not (bound < fx.vmax_scale*(1 - 1e-9) - 1):
        raise OverflowError(f"combine: result bound {bound:.3e} exceeds the fixed-point range {fx.vmax_scale:.3e}")
    if stats is not None: stats['bound'] = max(stats.get('bound', 0.0), bound)
    if raw: return acc                 # residues of the exact result at scale 2^-2F
    return fx.crt_out(acc, (r1 - r0, c1 - c0))

def blocks(n, b):
    return [(i, min(n, i + b)) for i in range(0, n, b)]

class Progress:
    def __init__(self, path):
        self.path = path
        self.done = set()
        if os.path.exists(path):
            with open(path) as fh: self.done = set(tuple(x) for x in json.load(fh)['done'])
    def mark(self, key):
        self.done.add(tuple(key))
        tmp = self.path + '.tmp'
        with open(tmp, 'w') as fh: json.dump({'done': sorted(self.done)}, fh)
        os.replace(tmp, self.path)

class OutOfTime(Exception):
    pass

RSS_LIMIT_MB = float(os.environ.get('OOC_RSS_LIMIT_MB', '600'))
def rss_now_mb():
    """current resident set size of this process (MB), from ps."""
    import subprocess
    out = subprocess.run(['ps', '-o', 'rss=', '-p', str(os.getpid())], capture_output=True, text=True).stdout.strip()
    return int(out)/1024.0 if out else 0.0
def check_budget(deadline, tag):
    """raise OutOfTime (the runner restarts the process, which resumes) when the time budget is used or when the
    resident memory exceeds RSS_LIMIT_MB (freed numpy/BLAS buffers are not always returned to the OS)."""
    if deadline is not None and time.time() > deadline: raise OutOfTime(tag + ' (time)')
    r = rss_now_mb()
    if r > RSS_LIMIT_MB: raise OutOfTime(f"{tag} (memory {r:.0f} MB)")

def run_blocks(fx, out_path, rows, cols, b, terms_fn, sym=False, deadline=None, log=print, tag=''):
    """Compute every output block (I, J) (only I <= J when sym; the mirror is written as the transpose).
    terms_fn(r0, r1, c0, c1) -> terms.  Resumable via out_path + '.prog'.  Returns stats."""
    part = out_path + '.part.npy'
    if os.path.exists(out_path): return {'resumed_complete': True}
    out = FileMat(part, 'w') if os.path.exists(part) else FileMat.create(part, rows, cols, fx.ND)
    assert out.shape == (rows, cols, fx.ND)
    prog = Progress(out_path + '.prog')
    stats = {}
    t0 = time.time()
    RB = blocks(rows, b); CB = blocks(cols, b)
    todo = [(i, j) for i in range(len(RB)) for j in range(len(CB)) if (i <= j or not sym)]
    for (i, j) in todo:
        if (i, j) in prog.done: continue
        try: check_budget(deadline, tag)
        except OutOfTime: out.close(); raise
        (r0, r1), (c0, c1) = RB[i], CB[j]
        tb = time.time()
        D = combine_block(fx, r0, r1, c0, c1, terms_fn(r0, r1, c0, c1), stats)
        out.write(r0, c0, D)
        if sym and i != j: out.write(c0, r0, np.transpose(D, (1, 0, 2)))
        out.flush()
        prog.mark((i, j))
        log(f"    {tag} block ({i},{j}) of {len(RB)}x{len(CB)}{' sym' if sym else ''}: {time.time()-tb:.1f}s (stage {time.time()-t0:.0f}s)")
    out.close()
    os.replace(part, out_path)
    if os.path.exists(out_path + '.prog'): os.remove(out_path + '.prog')
    return stats

def fm(path): return FileMat(path)
def view(path): return View(FileMat(path))

# ---------------- norms (rigorous upper bounds from digits) ----------------
def _rows(fmat, b):
    for r0, r1 in blocks(fmat.rows, b):
        yield fmat.read(r0, r1, 0, fmat.cols)

def frob_bound(fx, fmat, b=128):
    s = 0.0
    for D in _rows(fmat, b): s += fx.sumsq_bound(D)
    return math.sqrt(s)*(1 + 1e-9)

def maxabs(fx, fmat, b=128):
    return max(fx.maxabs_bound(D) for D in _rows(fmat, b))

def inf_bound(fx, fmat, b=128):
    """max row sum of |x| (upper bound): per entry |x| <= (|v| + 1) 2^(16 (ND-4) - F) + ..."""
    best = 0.0
    sc = math.ldexp(1.0, 16*(fx.ND - 4) - fx.F)
    for D in _rows(fmat, b):
        v = fx._top64(D).astype(np.float64)
        a = (np.maximum(np.abs(v), np.abs(v + 1)) + 2.0)*sc
        best = max(best, float(a.sum(axis=1).max()))
    return best*(1 + 1e-9)

def to_float_mat(fx, fmat, b=64):
    out = np.empty((fmat.rows, fmat.cols))
    for r0, r1 in blocks(fmat.rows, b):
        out[r0:r1] = fx.to_float(fmat.read(r0, r1, 0, fmat.cols))
    return out

def write_float(fx, path, M, b=64):
    rows, cols = M.shape
    out = FileMat.create(path + '.tmp.npy', rows, cols, fx.ND)
    for r0, r1 in blocks(rows, b): out.write(r0, 0, fx.from_float(M[r0:r1]))
    out.flush(); out.close()
    os.replace(path + '.tmp.npy', path)
