/-!
# Round 445 — node N1 (Weil positivity), round 10, check of a new version (R. Mori, Zenodo 24147570 = v6 of concept 24055205, 12 Nov 2026, 'The smallest eigenvalue of windowed Weil forms decays like e^{−4πμ} under RH and hypotheses on close pairs of zeros'): the bounds we audited are unchanged; Theorem 4's extremal measure σ* and its potential check numerically (p(σ*) = 2π); Theorem 3 (RH alone: λ_min ≥ exp(−exp(μ^ϖ)), ϖ > 1/2) rests on Carneiro–Chandee–Milinovich's S(t) bound; the per-pair cost in Prop. 7.6 has an implicit constant, so our lost-sample law cannot be compared; Zhu's conjectured −ln λ* = 2π²N(T*)/ln N(T*)(1 + o(1)) ∼ 2π²μ (2608.24827v2) fails under RH and Mori's hypotheses, though both laws fit at L ≤ 2; Zhu's certified upper bounds agree with Chekmachev's lower bounds; hypocoercivity, optimal quantization and prophet inequalities transfer no tool

* Node: N1 (10th round; R434 before). Round type: [check] a new version of a paper on the node line (Zenodo sweep, 12 Nov
  05:30 UTC; our audit morirate used v4 = 24116425 and paper III 24047156).
* [check] v4 -> v6 (pdftotext of both, difflib on the statements): Thm 1 remainder 4πe^{2/s}μ + c₃s³μ^{ϖB} -> 4πμ + c₃sμ^{ϖB};
  Thm 2 gains a rate; Props 6.1 and 7.6 take a free gap η ≤ 1/log μ (v4: 2/s); Cor 7.7's close-pair budget is weakened; Thm 5
  remainder -> 4πμ + C_τ log μ. Unchanged: λ_min ≤ 1.362·10⁹μ²³e^{−4πμ} (μ ≥ 5) and 1.24·10²³μ⁵(log μ)⁵e^{−4πμ} (μ ≥ 50),
  as on our N1 line: HOLDS.
* [falsify] Thm 4 (extremal problem, l.250-262; Lemmas 3.4-3.5): σ*(y) = −log|1 + φ(y)²| (φ the inverse Joukowski map) has
  ∫σ* = −7·10⁻³⁴ and ∫|σ*| = 4.000000; its potential at x = 0, 0.5, 1, 1.5, 1.99, 2.5, 4 is −6.283185, −4.013176, −2.151530,
  −0.750132, −0.002095, 0, 0, matching −π(√(4 − x²) − |x| arccos(|x|/2)) on |x| ≤ 2 to 6 digits, so p(σ*) = 2π: HOLDS.
* [check] Thm 3 (RH alone, l.223-233; threshold not explicit): Lemma 7.1 (at most (h/2π)log(t + 2) + (1/2 + ε)log t/log log t
  zeros in [t, t + h], 0 < h ≤ 1) follows in two lines from N(t) = ϑ(t)/π + 1 + S(t) and Carneiro–Chandee–Milinovich's
  |S(t)| ≤ (1/4 + o(1))log t/log log t under RH (Math. Ann. 2013): HOLDS (Goldston–Gonek 2007 not reread: arXiv API throttled).
* [check] Prop 7.6 / Cor 7.7 (l.1376-1470): each close labelled pair at scale h ≤ 2π/(log μ)³ costs c₄[log(1/h) + log log μ
  + log(1/η)] in −log λ_min with c₄ = c₄(C₀) implicit; our certified lost-sample cost (4 log(2 + √3) = 5.268 at y = 1, larger
  for low zeros) cannot be compared: NOT TESTABLE.
* [check] Zhu's conjecture (2608.24827v2, Conj. conj:law, TeX l.1926-1934): −ln λ*(L) = 2π²N(T*)/ln N(T*)(1 + o(1)),
  T* = 2πe^{2L}, with a fitted constant (Remark rem:rival: R₁ = 1.015-1.027 at L = 1.4-2.0, drifting up), so −ln λ* ∼ 2π²μ.
  Mori v6 §1.6: λ_min ≤ λ*, so under RH and (SEP_B), or (EXT) and (CENT), −ln λ* ≤ (4π + o(1))μ and the conjecture fails.
  At L = 2 Zhu's certified −ln λ* ≥ 650.47 lies between 2π²N/ln N = 638 and 4πμ = 686: finite sizes do not separate the laws.
* [falsify] Zhu's unconditional certified upper bounds (Table tab:main) against Chekmachev's certified lower bounds (ChekAudit):
  ε₂₅ = 1.17·10⁻¹²³ (half-width 1.609) ≤ λ_min(1.6) ≤ 8.61·10⁻¹²¹ = U(1.6), and ε₄₅ ≤ U(1.8) = 7.86·10⁻¹⁸⁶: consistent
  (ChekAudit's own bracket at c = 25 is tighter, 8.4·10⁻¹²³).
* Exploration (tool transfers, tested): Dolbeault–Mouhot–Schmeiser hypocoercivity (no generator with a coercive decomposition
  for ζ; kinetic theory tried in R369, R421, R427), Zador/Lloyd quantization (atomic zero measures; R29), Bruss's odds theorem
  (no independent trials; R376): no transfer.
-/

namespace Fayed.R445
end Fayed.R445
