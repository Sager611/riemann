R='/Users/adriansager/playground/riemann/'
T=open(R+'TREE.md',encoding='utf-8').read()
o="ceilaudit3's float lead (zeros near β = B charged in bulk, not one zero) gave the record ([CeilAudit3](Riemann/Deep/CeilAudit3/REPORT.md))."
n=o+(" Sibling criterion, Widder rungs (W_k(x) = (2k−1)! Σ_{[ρ]} m_ρ (q_ρ/(x + q_ρ)²)^k, q_ρ = ρ(1 − ρ); RH ⟺ W_k(x) ≥ 0 for all k ≥ 1, x > 0, by Widder's Stieltjes-transform theorem): "
 "W_k(x) > 0 for all x > 0 and k ≤ 0.4425 H² from RH to H = 3·10¹² and the HSW count (constant 9.4925, BW's typo note), 0.7660 H² with Bellotti–Wong's count, 0.8603 H² with their centred band, each sharp for its inputs (point sets with negative rungs at 0.4426, 0.7661, 0.8604 H²) "
 "(J. Huckstead, Zenodo 24210319, 14 Nov 2026, unrefereed; R526: inputs read correctly, Proposition E's point set legitimate, certificate r4 (85 items) rerun byte-identical; no implication either way with Li's range found).")
assert T.count(o)==1
open(R+'TREE.md','w',encoding='utf-8').write(T.replace(o,n,1)); print('ok')
