
import math, json
exec(open('/Users/adriansager/playground/riemann/.scratch/oddrlean/genlib.py').read())
d = json.load(open('/Users/adriansager/playground/riemann/.scratch/oddrlean/num3.json'))
EPSLO = F(*d['EPSLO']); EPSHI = F(*d['EPSHI']); L2LO = F(*d['L2LO']); L2HI = F(*d['L2HI'])
def msel(q, thr=F(1, 10**13)):
    a = q-1
    for m in range(1, 40):
        if ltail(a, m) < thr: return m
    return 40
def series_lo(q, m, D=14):
    a = q - 1; s = lsum(a, m); V = floorD(s, D); assert V <= s
    return V, f"hL.ge_of (q := {lq(q)}) (a := {lq(a)}) (by norm_num) (by norm_num) {m} (by norm_num [lsum, Finset.sum_range_succ])"
def series_hi(q, m, D=14):
    a = q - 1; s = lsum(a, m) + ltail(a, m); U = ceilD(s, D); assert U >= s
    return U, f"hL.le_of (q := {lq(q)}) (a := {lq(a)}) (by norm_num) (by norm_num) {m} (by norm_num [lsum, ltail, Finset.sum_range_succ])"
def short_lo(x, digs):
    e = math.floor(math.log10(float(x))); D = digs - 1 - e; assert 0 <= D <= 40
    return F(math.floor(x * F(10)**D), 10**D)
def short_hi(x, digs):
    e = math.floor(math.log10(float(x))); D = digs - 1 - e; assert 0 <= D <= 40
    return F(math.ceil(x * F(10)**D), 10**D)
HDR = lambda imports, ns_extra="": "\n".join("import " + i for i in imports) + "\n\nnamespace OddRobin\n\nopen ArithmeticFunction\nopen scoped ArithmeticFunction.sigma\n\nvariable {L : ℝ → ℝ}\n" + ns_extra
files = {}
def emit(fname, imports, body):
    files[fname] = HDR(imports) + "\n" + "\n".join(body) + "\n\nend OddRobin\n"
P61 = [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61]
def aN(p): return {3:4,5:3,7:2}.get(p,1)
def sig(p,a): return (p**(a+1)-1)//(p-1)
def val(f):
    n=1; s=1
    for p,a in f.items(): n*=p**a; s*=sig(p,a)
    return n,s
fN = {p:aN(p) for p in P61}
f3 = dict(fN); f3[3]=5
f11 = dict(fN); f11[11]=2
f67 = dict(fN); f67[5]=2; f67[67]=1
NN, SN = val(fN); m3, S3 = val(f3); m11, S11 = val(f11); m67, S67 = val(f67)
c0 = 67*NN
assert NN == 277093800961598968947975375 and c0 == 18565284664427130919514350125
# ---------------- SigVals: sigma of N and the three candidates (pure arithmetic)
B = []; W = B.append
W("lemma sig_mul {a b sa sb : ℕ} (h : Nat.Coprime a b) (ha : σ 1 a = sa) (hb : σ 1 b = sb) :")
W("    σ 1 (a * b) = sa * sb := by")
W("  rw [isMultiplicative_sigma.map_mul_of_coprime h, ha, hb]")
needed = set()
for f in [fN, f3, f11, f67]:
    for p,a in f.items(): needed.add((p,a))
for (p,a) in sorted(needed):
    W(f"lemma sv_{p}_{a} : σ 1 {p**a} = {sig(p,a)} := by")
    W(f"  have h := sigma_one_apply_prime_pow (p := {p}) (i := {a}) (by norm_num)")
    W(f"  norm_num [Finset.sum_range_succ] at h")
    W(f"  exact h")
def chain(name, f):
    items = sorted(f.items())
    # nested right-assoc product: x1 * (x2 * (... * xk))
    def prod_expr(its): 
        if len(its) == 1: return f"{its[0][0]**its[0][1]}"
        return f"{its[0][0]**its[0][1]} * ({prod_expr(its[1:])})" if len(its) > 2 else f"{its[0][0]**its[0][1]} * {its[1][0]**its[1][1]}"
    def sig_expr(its):
        if len(its) == 1: return f"{sig(*its[0])}"
        return f"{sig(*its[0])} * ({sig_expr(its[1:])})" if len(its) > 2 else f"{sig(*its[0])} * {sig(*its[1])}"
    k = len(items)
    prev = f"sv_{items[-1][0]}_{items[-1][1]}"
    for i in range(k-2, -1, -1):
        its = items[i:]
        W(f"lemma {name}_c{i} : σ 1 ({prod_expr(its)}) = {sig_expr(its)} :=")
        W(f"  sig_mul (by norm_num) sv_{items[i][0]}_{items[i][1]} {prev}")
        prev = f"{name}_c{i}"
    return prod_expr(items), prev
for (name, X, S, f) in [("NN", NN, SN, fN), ("m3", m3, S3, f3), ("m11", m11, S11, f11), ("m67", m67, S67, f67)]:
    pe, last = chain("ch" + name, f)
    if name != "NN":
        W(f"/-- {'3N' if name=='m3' else ('11N' if name=='m11' else '67N/5')}. -/")
        W(f"def {name} : ℕ := {X}")
    W(f"lemma sig_{name} : σ 1 {name} = {S} := by")
    W(f"  rw [show {name} = {pe} by norm_num [{name}], {last}]")
emit("SigVals", ["Riemann.Deep.OddRobin.Lean.Defs", "Mathlib.Tactic.NormNum.Prime", "Mathlib.Tactic.NormNum.GCD"], B)
# ---------------- log chains
RES = {}
def qlog(W, name, q, D=13):
    q = F(q); k = 0
    while q / 2**(k+1) >= 1: k += 1
    r = q / 2**k; m = msel(r)
    Lr, pfl = series_lo(r, m); Ur, pfu = series_hi(r, m)
    lo = floorD(k*L2LO + Lr, D); hi = ceilD(k*L2HI + Ur, D)
    W(f"theorem q_{name}_eq (hL : LogLike L) : L {lq(q)} = {k} * L 2 + L {lq(r)} := by")
    W(f"  have h := hL.two_pow_mul {k} (r := {lq(r)}) (by norm_num)")
    W(f"  rw [show (2 : ℝ) ^ {k} * {lq(r)} = {lq(q)} by norm_num] at h")
    W(f"  exact_mod_cast h")
    W(f"theorem q_{name}_ge (hL : LogLike L) : {lq(lo)} ≤ L {lq(q)} := by")
    W(f"  have h : {lq(Lr)} ≤ L {lq(r)} := {pfl}")
    W(f"  rw [q_{name}_eq hL]; linarith [L2_ge hL]")
    W(f"theorem q_{name}_le (hL : LogLike L) : L {lq(q)} ≤ {lq(hi)} := by")
    W(f"  have h : L {lq(r)} ≤ {lq(Ur)} := {pfu}")
    W(f"  rw [q_{name}_eq hL]; linarith [L2_le hL]")
    return lo, hi
def bignum(W, name, X, S, cast):
    k = X.bit_length()-1; r = F(X, 2**k)
    rlo = short_lo(r, 14); rhi = short_hi(r, 14)
    assert 2**k * rlo <= X <= 2**k * rhi
    m = msel(rhi)
    Lr, pfl = series_lo(rlo, m); Ur, pfu = series_hi(rhi, m)
    A_lo = floorD(k*L2LO + Lr, 10); A_hi = ceilD(k*L2HI + Ur, 10)
    W(f"theorem L{name}_ge (hL : LogLike L) : {lq(A_lo)} ≤ L ({cast} : ℝ) := by")
    W(f"  have h1 : (2 : ℝ) ^ {k} * {lq(rlo)} ≤ ({cast} : ℝ) := by norm_num [{name}]")
    W(f"  have h2 := hL.mono _ _ (by norm_num) h1")
    W(f"  rw [hL.two_pow_mul {k} (by norm_num)] at h2")
    W(f"  have h3 : {lq(Lr)} ≤ L {lq(rlo)} := {pfl}")
    W(f"  push_cast at h2; linarith [L2_ge hL]")
    W(f"theorem L{name}_le (hL : LogLike L) : L ({cast} : ℝ) ≤ {lq(A_hi)} := by")
    W(f"  have h1 : ({cast} : ℝ) ≤ (2 : ℝ) ^ {k} * {lq(rhi)} := by norm_num [{name}]")
    W(f"  have h2 := hL.mono _ _ (by norm_num [{name}]) h1")
    W(f"  rw [hL.two_pow_mul {k} (by norm_num)] at h2")
    W(f"  have h3 : L {lq(rhi)} ≤ {lq(Ur)} := {pfu}")
    W(f"  push_cast at h2; linarith [L2_le hL]")
    lo1 = qlog(W, f"{name}_A_lo", A_lo)[0]; hi1 = qlog(W, f"{name}_A_hi", A_hi)[1]
    B_lo = floorD(lo1, 10); B_hi = ceilD(hi1, 10)
    W(f"theorem LL{name}_ge (hL : LogLike L) : {lq(B_lo)} ≤ L (L ({cast} : ℝ)) := by")
    W(f"  have h := hL.mono _ _ (by norm_num) (L{name}_ge hL)")
    W(f"  linarith [q_{name}_A_lo_ge hL]")
    W(f"theorem LL{name}_le (hL : LogLike L) : L (L ({cast} : ℝ)) ≤ {lq(B_hi)} := by")
    W(f"  have h := hL.mono _ _ (by linarith [L{name}_ge hL]) (L{name}_le hL)")
    W(f"  linarith [q_{name}_A_hi_le hL]")
    lo2 = qlog(W, f"{name}_B_lo", B_lo)[0]; hi2 = qlog(W, f"{name}_B_hi", B_hi)[1]
    C_lo = floorD(lo2, 10); C_hi = ceilD(hi2, 10)
    W(f"theorem LLL{name}_ge (hL : LogLike L) : {lq(C_lo)} ≤ L (L (L ({cast} : ℝ))) := by")
    W(f"  have h := hL.mono _ _ (by norm_num) (LL{name}_ge hL)")
    W(f"  linarith [q_{name}_B_lo_ge hL]")
    W(f"theorem LLL{name}_le (hL : LogLike L) : L (L (L ({cast} : ℝ))) ≤ {lq(C_hi)} := by")
    W(f"  have h := hL.mono _ _ (by linarith [LL{name}_ge hL]) (LL{name}_le hL)")
    W(f"  linarith [q_{name}_B_hi_le hL]")
    Q = F(S, X); Qlo = short_lo(Q, 13); Qhi = short_hi(Q, 13)
    lq1 = qlog(W, f"{name}_Q_lo", Qlo)[0]; hq1 = qlog(W, f"{name}_Q_hi", Qhi)[1]
    W(f"theorem LQ{name}_le (hL : LogLike L) : L ((σ 1 {name} : ℝ) / ({cast} : ℝ)) ≤ {lq(hq1)} := by")
    W(f"  rw [sig_{name}]")
    W(f"  have h1 : (({S} : ℕ) : ℝ) / ({cast} : ℝ) ≤ {lq(Qhi)} := by norm_num [{name}]")
    W(f"  have h2 := hL.mono _ _ (by norm_num [{name}]) h1")
    W(f"  linarith [q_{name}_Q_hi_le hL]")
    W(f"theorem LQ{name}_ge (hL : LogLike L) : {lq(lq1)} ≤ L ((σ 1 {name} : ℝ) / ({cast} : ℝ)) := by")
    W(f"  rw [sig_{name}]")
    W(f"  have h1 : {lq(Qlo)} ≤ (({S} : ℕ) : ℝ) / ({cast} : ℝ) := by norm_num [{name}]")
    W(f"  have h2 := hL.mono _ _ (by norm_num) h1")
    W(f"  linarith [q_{name}_Q_lo_ge hL]")
    RES[name] = dict(A=(A_lo,A_hi), B=(B_lo,B_hi), C=(C_lo,C_hi), Q=(lq1,hq1))
imps = ["Riemann.Deep.OddRobin.Lean.NumA", "Riemann.Deep.OddRobin.Lean.SigVals"]
for (fn, name, X, S) in [("BigN", "NN", NN, SN), ("Big3", "m3", m3, S3), ("Big11", "m11", m11, S11), ("Big67", "m67", m67, S67)]:
    B = []
    bignum(B.append, name, X, S, name)
    emit(fn, imps, B)
# ---------------- Final numeric facts
GLO = F(5772145, 10**7); GHI = F(5776, 10**4)
B = []; W = B.append
dN = RES['NN']
d0hi = dN['Q'][1] - GLO + L2HI - dN['C'][0]
print('delta0 upper', float(d0hi)); assert d0hi < F(13,10000)
W("theorem delta0_lt (hL : LogLike L) {g : ℝ} (hg : (5772145 / 10000000 : ℝ) ≤ g) :")
W("    L ((σ 1 NN : ℝ) / (NN : ℝ)) - g + L 2 - L (L (L (NN : ℝ))) < 13 / 10000 := by")
W("  linarith [LQNN_le hL, L2_le hL, LLLNN_ge hL]")
for name in ['m3','m11','m67']:
    dd = RES[name]
    marg = dd['Q'][1] - (GLO - L2HI + dd['C'][0])
    print(name, 'margin', float(marg)); assert marg < 0
    W(f"theorem OR_{name} (hL : LogLike L) {{E g : ℝ}} (hE : 0 < E) (hLE : L E = g)")
    W(f"    (hg : (5772145 / 10000000 : ℝ) ≤ g) : ORf L E {name} := by")
    W(f"  rw [ORf_iff_log hL hE (by norm_num [{name}]) (by linarith [LL{name}_ge hL]), hLE]")
    W(f"  linarith [LQ{name}_le hL, L2_le hL, LLL{name}_ge hL]")
margN = dN['Q'][0] - (GHI - L2LO + dN['C'][1])
print('NN violation margin', float(margN)); assert margN > 0
W("theorem not_OR_NN (hL : LogLike L) {E g : ℝ} (hE : 0 < E) (hLE : L E = g)")
W("    (hg : g ≤ (5776 / 10000 : ℝ)) : ¬ ORf L E NN := by")
W("  rw [ORf_iff_log hL hE (by norm_num [NN]) (by linarith [LLNN_ge hL]), hLE]")
W("  exact not_lt.2 (by linarith [LQNN_ge hL, L2_ge hL, LLLNN_le hL])")
T = F(16289, 250)
assert 94*L2HI <= T and c0 <= 2**94
W("theorem Lc0_le (hL : LogLike L) : L (c0 : ℝ) ≤ (16289 / 250 : ℝ) := by")
W("  have h1 : (c0 : ℝ) ≤ (2 : ℝ) ^ 94 * 1 := by norm_num [c0]")
W("  have h2 := hL.mono _ _ (by norm_num [c0]) h1")
W("  rw [hL.two_pow_mul 94 (by norm_num), hL.map_one] at h2")
W("  push_cast at h2; linarith [L2_le hL]")
loT, hiT = qlog(W, "T", T)
T2 = ceilD(hiT, 6)
W(f"theorem LT_le (hL : LogLike L) : L (16289 / 250 : ℝ) ≤ {lq(T2)} := by linarith [q_T_le hL]")
assert EPSHI * T * T2 <= 1
print('eps T T2', float(EPSHI*T*T2))
W(f"theorem eps_TT2 (hL : LogLike L) : eps L * ((16289 / 250 : ℝ) * {lq(T2)}) ≤ 1 := by")
W(f"  have := eps_le hL")
W(f"  nlinarith")
W("theorem LNN_gt_one (hL : LogLike L) : 1 < L (NN : ℝ) := by linarith [LNN_ge hL]")
emit("Final", ["Riemann.Deep.OddRobin.Lean.BigN", "Riemann.Deep.OddRobin.Lean.Big3", "Riemann.Deep.OddRobin.Lean.Big11", "Riemann.Deep.OddRobin.Lean.Big67"], B)
json.dump({'T2': [T2.numerator, T2.denominator]}, open('/Users/adriansager/playground/riemann/.scratch/oddrlean/num4.json','w'))
for k,v in files.items():
    open(f'/Users/adriansager/playground/riemann/Riemann/Deep/OddRobin/Lean/{k}.lean','w',encoding='utf-8').write(v)
print('T2', T2)
for k,v in RES.items(): print(k, {kk: (float(a), float(b)) for kk,(a,b) in v.items()})

