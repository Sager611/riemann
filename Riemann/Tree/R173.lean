/-!
# Round 173 — Wormald's connected case survives; an OEIS growth conjecture looks false; a claimed 4D packing proof flagged

* Wormald's conjecture, connected case (open): the 16-vertex Sylvester graph, 96 bridged block graphs (n ≤ 28) and
  2000 random connected cubic graphs (n = 16, 20, 24) all split into two isomorphic spanning linear forests.
* OEIS A399309 and A399861 (a(n) = a(n−1) ± a(n−2), sign + iff A191232(n−1) = 1; conjectured
  lim |a(n)|^(1/n) = Viswanath's constant 1.1319882…): exactly, |a(n)|^(1/n) = 1.16968 at n = 99,999 and 1.15532 at
  n = 999,999, far above the constant for such n; growth fits are ill-conditioned. Next: the structure of A191232 and
  a Lyapunov-exponent estimate.
* 2609.25120 (Bhattacharjee et al., v2 24 Sep 2026) claims Δ₄ = π²/16, i.e. D₄ is an optimal packing in R⁴: m ≤ 22
  contacts by covering, m = 23 by a volume identity plus an SDP certificate, m = 24 by a classification of kissing
  configurations; both certificates are external data not in the source (TeX l. 553–557, 608–611, 2493–2507).
  A deep audit follows.
* Clean: 2609.15330 (Erdős Problem 767; n ≤ 6 agrees with the paper's formula); 2609.29250 (34 wild Erdős–Straus
  families, 1,001,000 instances, n ≤ 1.44·10⁸); 2609.29477 (an(C_{4k}) = 2k² − 1 for k ≤ 5); 2609.27788
  (R(K₄, Z₃) > 6 certificate).
-/

namespace Fayed.R173
end Fayed.R173
