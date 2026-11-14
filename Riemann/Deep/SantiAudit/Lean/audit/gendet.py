
import json
d=json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/SantiAudit/params_exp010.json'))
ps=[p for p in d['parameter_sets'] if p['nu']=='999/10000'][0]
assert ps['r_P']=='7707/10000' and ps['K']==100 and len(ps['x_j'])==100
out=[]
out.append('')
out.append('/-- x_j, j = 1..100, of the frozen EXP-010 set ν = 999/10000 (params_exp010.json; SL\'s Zenodo 23135354), verbatim. -/')
out.append('noncomputable def xs₁ : List ℝ := [')
items=[]
for x in ps['x_j']:
    if '/' in x:
        n,dd=x.split('/'); items.append('(%s / %s : ℝ)'%(n,dd))
    else: items.append('(%s : ℝ)'%x)
for i,it in enumerate(items):
    out.append('  '+it+(',' if i<len(items)-1 else ']'))
out.append('')
out.append('/-- The detector Q of E1.md §4 (degree 201): Q₁ = Qfam xs₁. -/')
out.append('noncomputable def Q₁ : ℝ → ℝ := Qfam xs₁')
out.append('')
out.append('/-- The detector P of E1.md §4 (degree 13): P₁(x) = S(r_P x)/S(r_P), r_P = 7707/10000. -/')
out.append('noncomputable def P₁ : ℝ → ℝ := Pfam (7707 / 10000)')
out.append('')
out.append('/-- Detector admissibility (E1.md §3 last row; SL Thm 1.2, 1.3 hypotheses), exact. -/')
out.append('theorem detector_admissible : IsPoly P₁ ∧ IsPoly Q₁ ∧ P₁ 0 = 0 ∧ P₁ 1 = 1 ∧ Q₁ 0 = 1 ∧ Q₁ 1 = 0 ∧')
out.append('    ∀ y, Q₁ y + Q₁ (1 - y) = 1 :=')
out.append('  ⟨isPoly_Pfam _, isPoly_Qfam _, Pfam_zero _, Pfam_one (by norm_num), Qfam_zero _, Qfam_one _, Qfam_symm _⟩')
out.append('')
out.append('end E1Lean')
open('/Users/adriansager/playground/riemann/.scratch/e1lean/detdata.lean','w').write('\n'.join(out)+'\n')
print(len(out), max(len(l) for l in out))

