/-!
# Round 125 — quick audits of 8 October math.NT papers (scouts at concurrency 4; the one flaw re-checked in the TeX)

* 2610.27343 Salinas et al., unconditional Λ_χ ≤ 953/3200 = 0.2978125 for all 242 primitive real characters
  with 3 ≤ |d| ≤ 400 (AI-assisted). The decisive Arb zero-free-region checks are asserted but no code or
  certificates are published, so it is not reproducible. RECORD TARGET: certify or improve Λ_χ for small |d|
  (e.g. χ₋₄) with our de Bruijn–Newman pipeline. LATER (after dbn5).
* 2610.03787 Lavrov, Dirichlet spectrum (μ ≤ 18.4809…; AI-assisted): Lemma ("d′ ≤ α") assumes without
  justification that (z_n, z_{n+1}) = (3,4). Pairs such as (4,1) occur (the partial quotients lie in
  {1,2,3,4}). Their products are easily below M₃M₄ (second factor < 2), and (4,3) is the mirror case, since the
  forbidden words are closed under reversal. Exposition gap; the conclusion stands.
* No flaw found:
  - 2610.08569 Johnston–Naik, least prime with given primitive root (GRH; small cases checked);
  - 2610.14257 Garcia–Grenié–Lee–Molteni, GRH residue bounds (checked on h = 1 fields);
  - 2610.15989 Ye, short-strip zero density (exponent arithmetic checked);
  - 2610.17675 Xu–Sang, sextic large sieve lower bound;
  - 2610.14127 Radomskii, additive complements;
  - 2610.10582, rough numbers in prime gaps (tables recomputed).
-/

namespace Fayed.R125
end Fayed.R125
