
import math, json
exec(open('/Users/adriansager/playground/riemann/.scratch/oddrlean/genlib.py').read())
def msel(q, thr=F(1, 10**13)):
    a = q-1
    for m in range(1, 40):
        if ltail(a, m) < thr: return m
    return 40
def series_lo(q, m, D=14):
    a = q - 1; s = lsum(a, m); V = floorD(s, D)
    assert V <= s
    return V, f"hL.ge_of (q := {lq(q)}) (a := {lq(a)}) (by norm_num) (by norm_num) {m} (by norm_num [lsum, Finset.sum_range_succ])"
def series_hi(q, m, D=14):
    a = q - 1; s = lsum(a, m) + ltail(a, m); U = ceilD(s, D)
    assert U >= s
    return U, f"hL.le_of (q := {lq(q)}) (a := {lq(a)}) (by norm_num) (by norm_num) {m} (by norm_num [lsum, ltail, Finset.sum_range_succ])"
HDR = lambda imports: "\n".join("import " + i for i in imports) + "\n\nnamespace OddRobin\n\nopen ArithmeticFunction\nopen scoped ArithmeticFunction.sigma\n\nvariable {L : ℝ → ℝ}\n"
files = {}
def emit(fname, imports, body):
    files[fname] = HDR(imports) + "\n" + "\n".join(body) + "\n\nend OddRobin\n"
# ---------- NumA: L 2, prime logs, eps
B = []; W = B.append
q = F(2); m = msel(q, F(1,10**14))
V2, pl = series_lo(q, m, 13); U2, pu = series_hi(q, m, 13)
L2LO, L2HI = V2, U2
W(f"theorem L2_ge (hL : LogLike L) : {lq(V2)} ≤ L 2 := {pl}")
W(f"theorem L2_le (hL : LogLike L) : L 2 ≤ {lq(U2)} := {pu}")
LOGP = {}
for p in [3,5,7,11,13,61,67,71,73]:
    k = p.bit_length()-1; r = F(p, 2**k)
    m = msel(r)
    Lr, pfl = series_lo(r, m); Ur, pfu = series_hi(r, m)
    lo = floorD(k*L2LO + Lr, 12); hi = ceilD(k*L2HI + Ur, 12)
    LOGP[p] = (lo, hi)
    W(f"theorem L{p}_eq (hL : LogLike L) : L ({p} : ℝ) = {k} * L 2 + L {lq(r)} := by")
    W(f"  have h := hL.two_pow_mul {k} (r := {lq(r)}) (by norm_num)")
    W(f"  rw [show (2 : ℝ) ^ {k} * {lq(r)} = {p} by norm_num] at h")
    W(f"  exact_mod_cast h")
    W(f"theorem L{p}_ge (hL : LogLike L) : {lq(lo)} ≤ L ({p} : ℝ) := by")
    W(f"  have h : {lq(Lr)} ≤ L {lq(r)} := {pfl}")
    W(f"  rw [L{p}_eq hL]; linarith [L2_ge hL]")
    W(f"theorem L{p}_le (hL : LogLike L) : L ({p} : ℝ) ≤ {lq(hi)} := by")
    W(f"  have h : L {lq(r)} ≤ {lq(Ur)} := {pfu}")
    W(f"  rw [L{p}_eq hL]; linarith [L2_le hL]")
q = F(68,67); m = msel(q, F(1,10**16))
L68, pf68l = series_lo(q, m, 16); U68, pf68u = series_hi(q, m, 16)
W(f"theorem L6867_ge (hL : LogLike L) : {lq(L68)} ≤ L (68 / 67 : ℝ) := {pf68l}")
W(f"theorem L6867_le (hL : LogLike L) : L (68 / 67 : ℝ) ≤ {lq(U68)} := {pf68u}")
lo67, hi67 = LOGP[67]
EPSLO = floorD(L68/hi67, 14); EPSHI = ceilD(U68/lo67, 14)
assert EPSLO*hi67 <= L68 and U68 <= EPSHI*lo67
W(f"theorem L67_pos (hL : LogLike L) : 0 < L (67 : ℝ) := by linarith [L67_ge hL]")
W(f"theorem eps_ge (hL : LogLike L) : {lq(EPSLO)} ≤ eps L := by")
W(f"  unfold eps; rw [le_div_iff₀ (L67_pos hL)]; linarith [L6867_ge hL, L67_le hL]")
W(f"theorem eps_le (hL : LogLike L) : eps L ≤ {lq(EPSHI)} := by")
W(f"  unfold eps; rw [div_le_iff₀ (L67_pos hL)]; linarith [L6867_le hL, L67_ge hL]")
W(f"theorem eps_pos (hL : LogLike L) : 0 < eps L := by linarith [eps_ge hL]")
emit("NumA", ["Riemann.Deep.OddRobin.Lean.Cost"], B)
# ---------- NumB / NumC: sigma values and inc bounds
def sig(p,a): return (p**(a+1)-1)//(p-1)
INC = {}
def inc_block(W, p, a):
    s = sig(p,a); q = 1 + F(1, p*s); mq = msel(q, F(1,10**15))
    Lq, pfl = series_lo(q, mq, 16); Uq, pfu = series_hi(q, mq, 16)
    lo_p, hi_p = LOGP[p]
    lo = floorD(Lq - EPSHI*hi_p, 14); hi = ceilD(Uq - EPSLO*lo_p, 14)
    INC[(p,a)] = (lo, hi)
    W(f"theorem sig_{p}_{a} : σ 1 ({p} ^ {a}) = {s} := by")
    W(f"  rw [sigma_one_apply_prime_pow (by norm_num)]; norm_num [Finset.sum_range_succ]")
    W(f"theorem inc_{p}_{a}_eq : inc L {p} {a} = L {lq(q)} - eps L * L ({p} : ℝ) := by")
    W(f"  rw [inc_val {p} {a} {s} sig_{p}_{a}]; norm_num")
    W(f"theorem inc_{p}_{a}_ge (hL : LogLike L) : {lq(lo)} ≤ inc L {p} {a} := by")
    W(f"  have h1 : {lq(Lq)} ≤ L {lq(q)} := {pfl}")
    W(f"  have h2 : eps L * L ({p} : ℝ) ≤ {lq(EPSHI)} * {lq(hi_p)} :=")
    W(f"    mul_le_mul (eps_le hL) (L{p}_le hL) (by linarith [L{p}_ge hL]) (by norm_num)")
    W(f"  rw [inc_{p}_{a}_eq]; linarith")
    W(f"theorem inc_{p}_{a}_le (hL : LogLike L) : inc L {p} {a} ≤ {lq(hi)} := by")
    W(f"  have h1 : L {lq(q)} ≤ {lq(Uq)} := {pfu}")
    W(f"  have h2 : {lq(EPSLO)} * {lq(lo_p)} ≤ eps L * L ({p} : ℝ) :=")
    W(f"    mul_le_mul (eps_ge hL) (L{p}_ge hL) (by norm_num) (eps_pos hL).le")
    W(f"  rw [inc_{p}_{a}_eq]; linarith")
B = []
for (p,a) in [(3,3),(3,4),(3,5),(5,1),(5,2),(5,3),(7,1),(7,2)]: inc_block(B.append, p, a)
emit("NumB", ["Riemann.Deep.OddRobin.Lean.NumA"], B)
B = []
for (p,a) in [(11,0),(11,1),(11,2),(67,1),(71,0),(71,1)]: inc_block(B.append, p, a)
W = B.append
W("theorem sig_67_0 : σ 1 (67 ^ 0) = 1 := by")
W("  rw [sigma_one_apply_prime_pow (by norm_num)]; norm_num")
W("theorem inc_67_0 (hL : LogLike L) : inc L 67 0 = 0 := by")
W("  rw [inc_val 67 0 1 sig_67_0]")
W("  have : (1 : ℝ) + 1 / (((67 : ℕ) : ℝ) * ((1 : ℕ) : ℝ)) = 68 / 67 := by norm_num")
W("  rw [this]; unfold eps; push_cast; rw [div_mul_cancel₀ _ (L67_pos hL).ne']; ring")
UNI = {}
for (name, q) in [("u61", F(62,61)), ("u182", F(183,182)), ("u73", F(74,73))]:
    mq = msel(q, F(1,10**15))
    Lq, pfl = series_lo(q, mq, 16); Uq, pfu = series_hi(q, mq, 16)
    UNI[name] = (Lq, Uq)
    W(f"theorem L_{name}_ge (hL : LogLike L) : {lq(Lq)} ≤ L {lq(q)} := {pfl}")
    W(f"theorem L_{name}_le (hL : LogLike L) : L {lq(q)} ≤ {lq(Uq)} := {pfu}")
emit("NumC", ["Riemann.Deep.OddRobin.Lean.NumA"], B)
for k,v in files.items():
    open(f'/Users/adriansager/playground/riemann/Riemann/Deep/OddRobin/Lean/{k}.lean','w',encoding='utf-8').write(v)
for k,v in INC.items(): print(k, float(v[0]), float(v[1]))
print('eps', float(EPSLO), float(EPSHI), 'L2', float(L2LO), float(L2HI))
print('u61 inc lower', float(UNI['u61'][0] - EPSHI*LOGP[61][1]))
print('u182 inc upper', float(UNI['u182'][1] - EPSLO*LOGP[13][0]))
print('u73 inc upper', float(UNI['u73'][1] - EPSLO*LOGP[73][0]))
json.dump({'EPSLO':[EPSLO.numerator,EPSLO.denominator],'EPSHI':[EPSHI.numerator,EPSHI.denominator],
  'L2LO':[L2LO.numerator,L2LO.denominator],'L2HI':[L2HI.numerator,L2HI.denominator],
  'LOGP':{str(p):[[v[0].numerator,v[0].denominator],[v[1].numerator,v[1].denominator]] for p,v in LOGP.items()},
  'INC':{f"{p}_{a}":[[v[0].numerator,v[0].denominator],[v[1].numerator,v[1].denominator]] for (p,a),v in INC.items()},
  'UNI':{k:[[v[0].numerator,v[0].denominator],[v[1].numerator,v[1].denominator]] for k,v in UNI.items()}},
  open('/Users/adriansager/playground/riemann/.scratch/oddrlean/num3.json','w'))

