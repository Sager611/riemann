import sys
R='/Users/adriansager/playground/riemann/'
T=open(R+'TREE.md',encoding='utf-8').read(); NG=open(R+'NOGO.md',encoding='utf-8').read()
L='[CeilAudit9](Riemann/Deep/CeilAudit9/REPORT.md)'
rerun = sys.argv[1]
E=[]
E.append(('N',"(orchestrator reran cue4.py from a copy: identical) [ExploreMalliavin](Riemann/Deep/ExploreMalliavin/REPORT.md)",
 "(orchestrator reran cue4.py from a copy: identical); ceilaudit9: ROBUST in 𝒮 + normalized microscopic four-point means (uniform in the shifts) + kurtosis + one-point law of Z; its witness fails macroscopic K(α) (K_ζ min 0.29–0.34 at α ≈ 3, t = 10⁴–10⁸; Cox-placed energy windows repair it, argued) and the one-point law of Z′ (the bulk carries 42–54% of ⟨Z′²⟩; repaired only under a hump-level CLT) "+L+" [ExploreMalliavin](Riemann/Deep/ExploreMalliavin/REPORT.md)"))
E.append(('N',"N5's worst-case-error form is a re-encoding [ExploreQMC](Riemann/Deep/ExploreQMC/REPORT.md)",
 "N5's worst-case-error form is a re-encoding; ceilaudit9: stated too weakly: |g(y)| ≤ g(0)/2 for |y| ≥ a makes two bumps optimal for periodic blocks, and LP placement under explore30's own constraints gives admissible integer blocks with W/(D·2π‖f‖²) = −0.10 on c = 13.46 and −0.39 on c = 20.1 at 3·10¹² (floats; LP thresholds c ≈ 5.2–11.0 by end and rounding model; anchored LP at t → ∞ negative from c = 13.46"+rerun+"): zero verification certifies no window c ≥ 13.46 in this model "+L+" [ExploreQMC](Riemann/Deep/ExploreQMC/REPORT.md)"))
E.append(('N',"so zero-atomicity buys ≤ 17% [ExploreSymm](Riemann/Deep/ExploreSymm/REPORT.md)",
 "so zero-atomicity buys ≤ 17%; ceilaudit9: ROBUST for subharmonic inequalities + line mean squares of a fixed VΨ: at T = 10⁴ twists B = 1 + Σ_{k≤6}β_k k^{−s} remove ≤ 4.1% of the concavity loss (a plain Levinson mollifier of the same length: 8.1%), two-piece splits ≤ 1.6%; atoms add nothing beyond T3 "+L+" [ExploreSymm](Riemann/Deep/ExploreSymm/REPORT.md)"))
E.append(('T',"open: mollified moments, macroscopic K(α), arithmetic inputs, [ExploreMalliavin](Riemann/Deep/ExploreMalliavin/REPORT.md), NOGO).",
 "open: mollified moments, macroscopic K(α), arithmetic inputs, [ExploreMalliavin](Riemann/Deep/ExploreMalliavin/REPORT.md), NOGO; ceilaudit9: ROBUST in that class; its witness misses macroscopic K(α) (repaired by Cox-placed windows, argued) and the one-point law of Z′ (repaired only under a hump-level CLT), "+L+")."))
E.append(('T',"makes Weil's form negative on every window with c ≥ 30 (explore30, NOGO)",
 "makes Weil's form negative on every window with c ≥ 30 (explore30; ceilaudit9: in floats already from c = 13.46, "+L+"; NOGO)"))
for w,o,n in E:
    s = T if w=='T' else NG
    assert s.count(o)==1, (w, s.count(o), o[:60])
assert 'ceilaudit9 (Orca' not in NG
for w,o,n in E:
    if w=='T': T=T.replace(o,n)
    else: NG=NG.replace(o,n)
NG = NG.rstrip('\n')+'\n'+"- ceilaudit9 (Orca, falsification of three barriers, all tested: explore62 N8, explore30 N1, explore45 N8/N7): all ROBUST in their exact classes; explore30 was stated too weakly (zero verification certifies no Weil window c ≥ 13.46 in its model, floats, not only c ≥ 30); explore62's witness misses macroscopic K(α) and the one-point law of Z′ "+L+"\n"
open(R+'TREE.md','w',encoding='utf-8').write(T); open(R+'NOGO.md','w',encoding='utf-8').write(NG); print('ok')
