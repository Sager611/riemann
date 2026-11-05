
import re
rows={}
for l in open('t1_out_raw.txt'):
    m=re.match(r'(\S+) s=(\S+) A=(\S+)',l)
    if m: rows.setdefault(m.group(1),{})[round(float(m.group(2)),2)]=float(m.group(3))
S=[0.68,0.70,0.71,0.72,0.75,0.78,0.80]
names={'H0':'H0 = GM Prop 11.1 (check: 15/(3+5σ))','H1':'H1: HB third term removed from the energy bound (R309 idea)',
'H1+C73':'H1 + Conj 7.3 also in S2 (HB third term gone everywhere)','C0.1':'C0.1: GM bottleneck \\|W\\|N^{4-4σ} cut by N^{-0.1}',
'C0.2':'C0.2: cut by N^{-0.2}','H3':'H3: minimal possible energy \\|W\\|^2+\\|W\\|^4/T','H3+C73':'H3 + Conj 7.3','H2':'H2: union of short APs, E ≤ \\|W\\|^2 N^{2-2σ} (alone)'}
print('| energy / Heath-Brown input | '+' | '.join('σ=%.2f'%s for s in S)+' | sup A |')
print('|---|'+'---|'*(len(S)+1))
for H in ['H0','H1','H1+C73','C0.1','C0.2','H3','H3+C73','H2']:
    d=rows[H]; sup=max(d.values())
    print('| '+names[H]+' | '+' | '.join('%.3f'%d[s] for s in S)+' | %.4f |'%sup)
print('| GM formula 15/(3+5σ); Ingham 3/(2−σ) | '+' | '.join('%.3f / %.3f'%(15/(3+5*s),3/(2-s)) for s in S)+' | 30/13 = %.4f |'%(30/13))
