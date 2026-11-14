
from fractions import Fraction as Fr
import math
D=10**20
def dn(x): return Fr(math.floor(x*D), D)
def up(x): return Fr(math.ceil(x*D), D)
def L(x): 
    x=Fr(x)
    return '(%d / %d : ℝ)'%(x.numerator,x.denominator) if x.denominator!=1 else '(%d : ℝ)'%x.numerator
slo=Fr(141421356237,10**11); shi=Fr(141421356238,10**11)
xlo=Fr(353553,10**6); a=xlo/8
sl=dn(a-a**3/6-a**5/100); sh=up(a-a**3/6+a**5/100)
cl=dn(1-a**2/2-a**4*Fr(5,96)); ch=up(1-a**2/2+a**4*Fr(5,96))
st=[(sl,sh,cl,ch)]
for i in range(3):
    assert sl>=0 and cl>=0
    nsl=dn(2*sl*cl); nsh=up(2*sh*ch); ncl=dn(1-2*sh**2); nch=up(1-2*sl**2)
    assert nsl<=2*sl*cl and 2*sh*ch<=nsh and ncl<=1-2*sh**2 and 1-2*sl**2<=nch
    sl,sh,cl,ch=nsl,nsh,ncl,nch; st.append((sl,sh,cl,ch))
U=up(ch/sl)
CLO=dn(Fr(7,4)-U/slo)
print('U',float(U),'CLO',float(CLO))
A=['353553 / 8000000','353553 / 4000000','353553 / 2000000','353553 / 1000000']
o=[]
o.append('')
o.append('/-- sin, cos at 353553/10⁶ (≤ 1/(2√2)): Taylor at 353553/(8·10⁶), then three doublings. -/')
o.append('theorem trig_xlo : %s ≤ F.sin (353553 / 1000000) ∧ F.cos (353553 / 1000000) ≤ %s := by'%(L(st[3][0]),L(st[3][3])))
o.append('  have b0 := base_bounds F (a := 353553 / 8000000) (by norm_num) (by norm_num)')
o.append('  have s0 : %s ≤ F.sin (353553 / 8000000) ∧ F.sin (353553 / 8000000) ≤ %s :='%(L(st[0][0]),L(st[0][1])))
o.append('    ⟨le_trans (by norm_num) b0.1.1, le_trans b0.1.2 (by norm_num)⟩')
o.append('  have c0 : %s ≤ F.cos (353553 / 8000000) ∧ F.cos (353553 / 8000000) ≤ %s :='%(L(st[0][2]),L(st[0][3])))
o.append('    ⟨le_trans (by norm_num) b0.2.1, le_trans b0.2.2 (by norm_num)⟩')
prev=('s0','c0')
for i in range(1,4):
    s=st[i]
    o.append('  have d%d := dbl F (sl\' := %s) (sh\' := %s) (cl\' := %s) (ch\' := %s) %s %s'%(i,L(s[0]),L(s[1]),L(s[2]),L(s[3]),prev[0],prev[1]))
    o.append('    (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)')
    o.append('  rw [show (2 : ℝ) * (%s) = %s by norm_num] at d%d'%(A[i-1],A[i],i))
    prev=('d%d.1'%i,'d%d.2'%i)
o.append('  exact ⟨d3.1.1, d3.2.2⟩')
o.append('')
o.append('/-- c(1/2) ≥ c_lo (E1.md §4: c(1/2) = −0.16596385…; c_lo = %.12f). -/'%float(CLO))
o.append('theorem cW_half_ge : %s ≤ cW F (1 / 2) := by'%L(CLO))
o.append('  obtain ⟨hs, hss, hs1⟩ := sqrt2_facts F')
o.append('  obtain ⟨hsl, hsh⟩ := sqrt2_bounds F')
o.append('  have ht := trig_xlo F')
o.append('  have hcot : cot F (353553 / 1000000) ≤ %s := by'%L(U))
o.append('    unfold cot')
o.append('    calc F.cos (353553 / 1000000) / F.sin (353553 / 1000000) ≤ %s / %s :='%(L(st[3][3]),L(st[3][0])))
o.append('          div_le_div₀ (by norm_num) ht.2 (by norm_num) ht.1')
o.append('      _ ≤ %s := by norm_num'%L(U))
o.append('  set x₀ := (1 / 2 : ℝ) / F.sqrt 2 with hx₀')
o.append('  have hxlo : (353553 / 1000000 : ℝ) ≤ x₀ := by rw [hx₀, le_div_iff₀ hs]; linarith')
o.append('  have hx1 : x₀ ≤ 1 := by rw [hx₀, div_le_iff₀ hs]; linarith')
o.append('  have hcotx : cot F x₀ ≤ %s := by'%L(U))
o.append('    rcases hxlo.lt_or_eq with hlt | heq')
o.append('    · have := cot_sub_gt F (by norm_num) hlt hx1; linarith')
o.append('    · rw [← heq]; exact hcot')
o.append('  unfold cW')
o.append('  rw [← hx₀]')
o.append('  have h1 : cot F x₀ / F.sqrt 2 ≤ %s / F.sqrt 2 := div_le_div_of_nonneg_right hcotx hs.le'%L(U))
o.append('  have h2 : %s / F.sqrt 2 ≤ %s / (141421356237 / 100000000000) :='%(L(U),L(U)))
o.append('    div_le_div_of_nonneg_left (by norm_num) (by norm_num) hsl')
o.append('  have h3 : %s ≤ 2 - 1 / 2 / 2 - %s / (141421356237 / 100000000000) := by norm_num'%(L(CLO),L(U)))
o.append('  linarith')
open('/Users/adriansager/playground/riemann/.scratch/e1lean/numB.lean','w').write('\n'.join(o)+'\n')
open('/Users/adriansager/playground/riemann/.scratch/e1lean/CLO.txt','w').write(L(CLO)+'\n'+L(U)+'\n')
print(len(o), max(len(x) for x in o))

