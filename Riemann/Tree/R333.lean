/-!
# Round 333 — node N3 (Li/Keiper), round 4: Fabbian's inputs rechecked (RH to T + HSW count only); our Li certificate moved off Trudgian 2014's flawed zero count (numbers unchanged) and extended to 106.25 T^2 with the Lean-checked half-plane (routine)

* Node: N3 Li / Keiper-Li (fourth round; R315 before; dirli works on the Dirichlet-L branch).
* [check] Fabbian (Zenodo 22839423) remains the latest certified range; its Theorem 1 uses only RH to height T
  (Platt-Trudgian) and the Hasanalizade-Shen-Wong zero count ('No hypothesis beyond these inputs', sec1_intro l.50-59);
  Brown's hypothesis B_r(T) is a separate weaker variant that holds if RH holds to T (l.69-76). The adversarial scout's
  'GAP' (a strip condition above T) misread that variant: HOLDS.
* Bug found and fixed: Fabbian's Remark rem:trudgian (sec2_prelim l.62-68) records that Trudgian 2014's explicit N(T) bound
  has a reported error (HSW 2022 fn. 2, following BMOR 2021). Our Li200k certificate used it. verify.py now uses HSW 2022
  Cor. 1.2 (|N - Ntilde| <= 0.1038 log t + 0.2573 log log t + 9.3675 + 7/8); rerun (80 s): every certified N is unchanged
  (A: 78.713642 T0^2; C, 7/8 half-plane: 952,080,909,647,480,761,027,133,440 = 105.774404 T0^2). RamRepair's zero-count use
  (log X + 40 for X >= 1e10) already holds with HSW. dirli told to use BMOR for Dirichlet L.
* Computation (routine, COROLLARIES): new variant F with the Lean-checked half-plane a = 0.37338196661320752: lambda_n > 0
  for all n <= 956,318,469,361,119,813,265,522,688 = 106.245188 T0^2 (the scout's 1/a scaling 106.2328 T0^2 was close).
* State of the art: no 2024-2026 unconditional Li-coefficient asymptotic or new range on arXiv (Johansson's certified values
  stand).
* Exploration: OPUC / Caratheodory positivity needs a positive moment representation of lambda_n that would already encode
  RH (no idea); matrix models / topological recursion (2609.04908 builds finite Hermitian matrices from known zeros) and
  Connes-Kreimer renormalization (2001.08722, 1804.05568; R298): no idea.
-/

namespace Fayed.R333
end Fayed.R333
