
import json
exec(open('/Users/adriansager/playground/riemann/.scratch/oddrlean/genlib.py').read())
d = json.load(open('/Users/adriansager/playground/riemann/.scratch/oddrlean/num3.json'))
L2LO = F(*d['L2LO']); L2HI = F(*d['L2HI'])
def msel(q, thr=F(1, 10**13)):
    a = q-1
    for m in range(1, 40):
        if ltail(a, m) < thr: return m
    return 40
r = F(21,16); m = msel(r)
a = r - 1
Lr = floorD(lsum(a, m), 14); Ur = ceilD(lsum(a,m)+ltail(a,m), 14)
# L 2 bounds proven here directly (to keep this file independent of NumA)
q2 = F(2); m2 = msel(q2, F(1,10**14))
V2 = floorD(lsum(F(1), m2), 13); U2 = ceilD(lsum(F(1),m2)+ltail(F(1),m2), 13)
assert V2 == L2LO and U2 == L2HI
H = F(7381,2520); M = 10**7
dhi = H - (3*L2LO + Lr)  # d_10 upper
dlo = H - (3*L2HI + Ur)  # d_10 lower
ghi = dhi + F(1, 2*M); glo = dlo - F(1,2646) - F(1, 2*M+1)
GLO = F(5772145, 10**7); GHI = F(5776, 10**4)
print(float(glo), float(ghi)); assert glo >= GLO and ghi <= GHI
def ge(q,a,m): return f"hL.ge_of (q := {lq(q)}) (a := {lq(a)}) (by norm_num) (by norm_num) {m} (by norm_num [lsum, Finset.sum_range_succ])"
def le(q,a,m): return f"hL.le_of (q := {lq(q)}) (a := {lq(a)}) (by norm_num) (by norm_num) {m} (by norm_num [lsum, ltail, Finset.sum_range_succ])"
txt = f'''theorem gL2_ge (hL : LogLike L) : {lq(V2)} ≤ L 2 := {ge(F(2),F(1),m2)}

theorem gL2_le (hL : LogLike L) : L 2 ≤ {lq(U2)} := {le(F(2),F(1),m2)}

lemma L_21_16_ge (hL : LogLike L) : {lq(Lr)} ≤ L (21 / 16 : ℝ) := {ge(r,a,m)}

lemma L_21_16_le (hL : LogLike L) : L (21 / 16 : ℝ) ≤ {lq(Ur)} := {le(r,a,m)}

lemma dSeq_ten (hL : LogLike L) : dSeq L 10 = 7381 / 2520 - (3 * L 2 + L (21 / 16 : ℝ)) := by
  rw [← L_21_2_eq hL]
  unfold dSeq
  rw [harmonic_ten]
  norm_num

/-- 0.5772145 ≤ g ≤ 0.5776 from Mathlib's two bounds for gamma (at M = 10^7). -/
theorem gamma_bounds (hL : LogLike L) {{g : ℝ}}
    (hlo : ∀ n : ℕ, (harmonic n : ℝ) - L ((n : ℝ) + 1) < g)
    (hup : ∀ n : ℕ, 0 < n → g < (harmonic n : ℝ) - L n) :
    (5772145 / 10000000 : ℝ) ≤ g ∧ g ≤ (5776 / 10000 : ℝ) := by
  have h1 := le_g hL hlo (n := 10) (M := 10000000) (by norm_num)
  have h2 := g_le hL hup (n := 10) (M := 10000000) (by norm_num) (by norm_num)
  have hc : cSeq 10 = 1 / 2646 := by unfold cSeq; norm_num
  rw [dSeq_ten hL, hc] at h1
  rw [dSeq_ten hL] at h2
  have a1 := gL2_ge hL
  have a2 := gL2_le hL
  have a3 := L_21_16_ge hL
  have a4 := L_21_16_le hL
  push_cast at h1 h2
  constructor <;> linarith'''
t = open('/Users/adriansager/playground/riemann/.scratch/oddrlean/GammaL.tmpl',encoding='utf-8').read().replace('GAMMA_NUMERIC', txt)
open('/Users/adriansager/playground/riemann/Riemann/Deep/OddRobin/Lean/GammaL.lean','w',encoding='utf-8').write(t)
print('ok')

