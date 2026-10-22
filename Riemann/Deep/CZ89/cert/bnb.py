
# Rigorous branch-and-bound certificate for Lemma 5.2 (seven-point inequality) of
# Cobeli-Zaharescu arXiv 2610.23467:  S(g) >= 11/250 for all g in [0,inf)^6.
# If sum g >= 8.8 the linear term alone gives S >= 0.044, so we cover R = {g >= 0, sum g <= 8.8}.
# Inputs: tables*.npz from table.py (Arb enclosures of q=K_*^2 and q' at grid points x_m = m*h_t;
# per-cell [x_m, x_{m+1}] lower bound of q and upper bound of |q''|).
# Box coordinates are integers in units h_b = 2 h_t (dyadic), so box centres lie on the table grid.
# Tests (both rigorous lower bounds for S on a box):
#   L1 = sum_ij w_ij min_{cells covering J_ij} qmin + p*sum_k l_k
#   L2 = S(c) - sum_k |D_k| r_k - sum_ij w_ij (M2_ij/2) r_ij^2,  M2_ij = max |q''| on the pair range
# Floating point: every sum has < 40 terms; we add gamma = 64u times the sum of absolute values (u = 2^-53).
import os, sys, time, json, hashlib, numpy as np
from fractions import Fraction
TABLE = os.environ.get("CZ_TABLE", "tables.npz")
T = np.load(TABLE)
qlo, q1lo, q1hi, cq2, cqmin = T["qlo"], T["q1lo"], T["q1hi"], T["cq2"], T["cqmin"]
HT = float(T["HT"]); HB = 2*HT
NC = len(cqmin)
SUMB = int(Fraction(44, 5) / Fraction(HB)) + 1       # sum g <= 8.8  <=>  sum(g/h_b) <= 8.8/h_b < SUMB
P = 1/200.0
TARGET = float(os.environ.get("CZ_TARGET", "0.044"))
GLOBALM2 = float(os.environ.get("CZ_GLOBALM2", "0"))  # >0: replace local |q''| bound by this constant (paper uses 40)
GAM = 64 * 2.0**-53
rows = {0:[0.306,0.645,1.000,1.000,1.000,2.000],1:[0.347,0.047,0.000,0.000,1.000],2:[0.347,0.616,0.000,1.000],3:[0.347,0.047,1.000],4:[0.347,0.645],5:[0.306]}
pairs = []; w = []
for i,v in rows.items():
    for t,val in enumerate(v):
        pairs.append((i, i+1+t)); w.append(val)
w = np.array(w); NPAIR = len(pairs)
assert NPAIR == 21 and abs(w.sum() - 12) < 1e-12
for dlt in range(1, 7):
    assert abs(sum(wv for (i,j), wv in zip(pairs, w) if j - i == dlt) - 2) < 1e-12
A = np.zeros((NPAIR,6))
for p_,(i,j) in enumerate(pairs): A[p_, i:j] = 1.0      # pair (i,j) spans gaps g_{i+1..j}
AT = A.T.copy()
def build(arr, fn):
    lv = [arr]; k = 1
    while (1 << k) <= len(arr):
        prev = lv[-1]; h = 1 << (k-1)
        lv.append(fn(prev[:len(prev)-h], prev[h:])); k += 1
    out = np.full((len(lv), len(arr)), np.nan)
    for k_, a_ in enumerate(lv): out[k_, :len(a_)] = a_
    return out
STMIN = build(cqmin, np.minimum); STMAX = build(cq2, np.maximum)
LOG2 = np.zeros(NC+2, dtype=np.int64)
for n in range(2, NC+2): LOG2[n] = LOG2[n//2] + 1
def rmin(lo, hi):
    k = LOG2[hi - lo + 1]
    return np.minimum(STMIN[k, lo], STMIN[k, hi - np.left_shift(1, k) + 1])
def rmax(lo, hi):
    k = LOG2[hi - lo + 1]
    return np.maximum(STMAX[k, lo], STMAX[k, hi - np.left_shift(1, k) + 1])
def tighten(L, U):
    tot = L.sum(axis=1, keepdims=True)
    return np.minimum(U, SUMB - (tot - L))
def evaluate(L, U):
    Lf = L.astype(np.float64); Uf = U.astype(np.float64)
    SL = Lf @ AT; SU = Uf @ AT                        # exact integers
    totL = Lf.sum(axis=1, keepdims=True)
    Jhi = np.minimum(SU, SUMB - (totL - SL))          # pair range clipped by the sum constraint
    lo = (2*SL).astype(np.int64); hi = np.maximum((2*Jhi).astype(np.int64) - 1, lo)
    qm = rmin(lo, hi)
    t1 = qm @ w; t2 = P*HB*totL[:,0]
    L1 = t1 + t2
    err1 = GAM*(np.abs(t1) + t2) + 1e-300
    ok = (2*SU).max(axis=1) <= NC
    C = Lf + Uf                                       # centres in h_t units
    s = np.minimum((C @ AT).astype(np.int64), NC)
    lo2 = np.minimum((2*SL).astype(np.int64), NC-1)
    hi2 = np.minimum(np.maximum((2*SU).astype(np.int64)-1, lo2), NC-1)
    M2 = rmax(lo2, hi2)
    if GLOBALM2 > 0: M2 = np.full_like(M2, GLOBALM2)
    Sc = qlo[s] @ w + P*HT*C.sum(axis=1)
    Dlo = P + (q1lo[s]*w) @ A; Dhi = P + (q1hi[s]*w) @ A
    eD = GAM*(P + (np.maximum(np.abs(q1lo[s]), np.abs(q1hi[s]))*w) @ A)
    Dabs = np.maximum(np.abs(Dlo), np.abs(Dhi)) + eD
    r = (Uf - Lf)*HT
    rij = r @ AT
    lin = (Dabs*r).sum(axis=1)
    quad = (rij*rij*M2) @ w * 0.5
    L2 = Sc - lin - quad
    err2 = GAM*(np.abs(Sc) + lin + quad) + 1e-300
    L2 = np.where(ok, L2, -np.inf)
    return L1 >= TARGET + err1, L2 >= TARGET + err2, L1, L2
def run(L0, U0, batch=8192, maxbox=10**12, tmax=1e9, verbose=True, label=""):
    t0 = time.time()
    stackL = [np.array(L0, dtype=np.int64).reshape(-1,6)]; stackU = [np.array(U0, dtype=np.int64).reshape(-1,6)]
    st = dict(processed=0, acc_L1=0, acc_L2=0, empty=0, split=0, unresolved=0)
    unres = []; minacc = np.inf; lastp = 0; status = "COMPLETE"
    while stackL:
        L = stackL.pop(); U = stackU.pop()
        if len(L) > batch:
            stackL.append(L[batch:]); stackU.append(U[batch:]); L = L[:batch]; U = U[:batch]
        st["processed"] += len(L)
        U = tighten(L, U)
        ne = (L.sum(axis=1) < SUMB) & np.all(U >= L, axis=1)
        st["empty"] += int((~ne).sum()); L = L[ne]; U = U[ne]
        if len(L) == 0: continue
        a1, a2, L1, L2 = evaluate(L, U)
        acc = a1 | a2
        st["acc_L1"] += int(a1.sum()); st["acc_L2"] += int((a2 & ~a1).sum())
        if acc.any():
            minacc = min(minacc, float(np.maximum(np.where(a1, L1, -np.inf), np.where(a2, L2, -np.inf))[acc].min()))
        F = ~acc
        if F.any():
            L = L[F]; U = U[F]
            wd = U - L; k = wd.argmax(axis=1); ar = np.arange(len(k)); mw = wd[ar, k]
            bad = mw <= 1                              # cannot bisect further on this grid
            if bad.any():
                st["unresolved"] += int(bad.sum())
                for i in np.where(bad)[0][:3]: unres.append(((L[i]*HB).tolist(), (U[i]*HB).tolist()))
                L = L[~bad]; U = U[~bad]; k = k[~bad]; ar = np.arange(len(k))
            st["split"] += len(L)
            mid = L[ar, k] + (U - L)[ar, k]//2
            La = L.copy(); Ua = U.copy(); Ua[ar, k] = mid
            Lb = L.copy(); Ub = U.copy(); Lb[ar, k] = mid
            stackL.append(np.concatenate([La, Lb])); stackU.append(np.concatenate([Ua, Ub]))
        if verbose and st["processed"] - lastp >= 5_000_000:
            lastp = st["processed"]
            print(label, st, "stack", sum(len(x) for x in stackL), "t %.1f" % (time.time()-t0), flush=True)
        if st["processed"] > maxbox or time.time()-t0 > tmax:
            status = "STOPPED_LIMIT"; break
        if st["unresolved"] > 50:
            status = "STOPPED_UNRESOLVED"; break
    left = sum(len(x) for x in stackL)
    if status == "COMPLETE" and (st["unresolved"] > 0 or left > 0): status = "FAILED"
    res = dict(label=label, status=status, table=TABLE, table_sha256=hashlib.sha256(open(TABLE,'rb').read()).hexdigest()[:16],
               h_t=HT, h_b=HB, SUMB=SUMB, target=TARGET, globalM2=GLOBALM2, stack_left=left,
               min_accepted_lower_bound=minacc, seconds=round(time.time()-t0,1), **st, unresolved_examples=unres[:3])
    print(json.dumps(res), flush=True)
    return res
if __name__ == "__main__":
    mode = sys.argv[1]; tmax = float(sys.argv[2]) if len(sys.argv) > 2 else 1e9
    maxbox = int(float(sys.argv[3])) if len(sys.argv) > 3 else 10**12
    if mode == "full":
        run(np.zeros(6, dtype=np.int64), np.full(6, SUMB, dtype=np.int64), tmax=tmax, maxbox=maxbox, label="full R={g>=0,sum<=8.8}")
    elif mode == "omega1":
        run(np.full(6, int(np.floor(0.743/HB)), dtype=np.int64), np.full(6, int(np.ceil(5.085/HB)), dtype=np.int64), tmax=tmax, maxbox=maxbox, label="Omega'")
    elif mode == "local":
        g = np.array([1.0837079364,1.8939190725,1.0413580806,1.0248802119,1.0289099949,1.0255587059])
        rad = float(sys.argv[4])
        run(np.floor((g-rad)/HB).astype(np.int64), np.ceil((g+rad)/HB).astype(np.int64), tmax=tmax, maxbox=maxbox, label="cube radius %g at g*" % rad)

