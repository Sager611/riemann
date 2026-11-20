NG = open('NOGO.md', encoding='utf-8').read(); T = open('TREE.md', encoding='utf-8').read()
R = open('Riemann/Deep/ExplorePositivity/REPORT.md', encoding='utf-8').read().split('\n')
L = [l for l in R if l.startswith('- explore63 (Orca, exploration: special-function positivity')]
assert len(L) == 1
line = L[0]
tail = ' [ExplorePositivity](Riemann/Deep/ExplorePositivity/REPORT.md)'
assert line.endswith(tail)
line = line[:-len(tail)] + ' (orchestrator reran mp_c2n.py from a copy: identical c₀, c₂, c₄; ratio 0.3858)' + tail
t_old = "Barriers: uniformity in d ≡ RH; asymptotic methods give only n ≫ d."
t_new = ("Barriers: uniformity in d ≡ RH; asymptotic methods give only n ≫ d; orthogonal-basis Pólya–Turán transfers (Hermite, Meixner–Pollaczek, continuous Hahn) certify nothing: "
 "Newton fails at every order N ≥ 2 for every positive Φ, their approximants have non-real zeros from N = 2, and Jensen's damping is Maclaurin-extremal among finite multipliers (explore63, 20 Nov 2026, [ExplorePositivity](Riemann/Deep/ExplorePositivity/REPORT.md), NOGO).")
assert 'explore63 (Orca' not in NG and T.count(t_old) == 1
if not NG.endswith('\n'): NG += '\n'
open('NOGO.md', 'w', encoding='utf-8').write(NG + line + '\n'); open('TREE.md', 'w', encoding='utf-8').write(T.replace(t_old, t_new, 1)); print('ok')
