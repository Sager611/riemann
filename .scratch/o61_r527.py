R='/Users/adriansager/playground/riemann/'
T=open(R+'TREE.md',encoding='utf-8').read()
o="Jensen's damping is Maclaurin-extremal among finite multipliers (explore63, 20 Nov 2026, [ExplorePositivity](Riemann/Deep/ExplorePositivity/REPORT.md), NOGO)."
n=o+(" Total positivity of Riemann's kernel: Φ has exact Pólya-frequency order 4 (every translation minor det[Φ(x_i − y_j)] of order ≤ 4 at ordered nodes is > 0, some order-5 minor is < 0; J. Rainstar, Zenodo 21487371, July 2026, Lean 4: PFOrderExactly globalRiemannKernel 4, kernel 2Φ(t/2), unrefereed; "
 "answers Michalowski's Problem 1, arXiv 2602.20313, which certified the order-5 failure; replayed by D. Finlay, Zenodo 24082460; R527: definitions faithful up to that PF-preserving normalisation, static axiom scan clean (build not rerun here), order-5 witness reproduced (−1.9·10⁻⁷ of the diagonal product), 780 adversarial order-4 minors positive); no RH implication (Φ cannot be PF∞).")
assert T.count(o)==1
open(R+'TREE.md','w',encoding='utf-8').write(T.replace(o,n,1)); print('ok')
