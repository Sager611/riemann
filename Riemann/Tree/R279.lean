/-!
# Round 279 — node N7 (zero-free regions, zero density): two newer records found; Liu's constant reproduced

* Node: N7 zero-free regions and zero-density constants.
* State of the art (scouts, orchestrator read the abstracts): the classical-shape record is R_0 = 4.896 for t >= 3
  (Bellotti-Trudgian-Yang, arXiv 2603.21490: zeta(sigma + it) != 0 for sigma >= 1 - 1/(4.896 log t)); our [check] value
  5.558691 (Mossinghoff-Trudgian-Yang) was outdated; TREE.md N7 updated. Chang-Wang (arXiv 2611.00924, Nov 2026) claim
  the zero-free half-plane boundary B = B* - 1/160000 = 0.87446205769353 (Liu's B* = 0.874886647040), stating that the
  conclusion is conditional on cited general-field theta, detector, prior-bound and Mellin input theorems: to audit.
* Adversarial: Liu 2610.26835's constant B = 0.87488664703996 reproduced from its stated lemmas (quartic root; kappa =
  2B - 1 = 0.74977, 6 kappa - 1 = 3.4986 >= 10/3), so our R131 flag does not affect the arithmetic.
* Wild: a degree-16 non-negative cosine polynomial with a_1/a_0 = 1.7502 > 1.7412 (the MTY polynomial's ratio) was found,
  but a_1/a_0 is not the objective that sets R_0; no improved zero-free constant follows.
* Exploration idea (decoupling / large values): use the short-strip density bound of 2610.15989
  (N(sigma; T, T + T^a) << T^{15a(1-sigma)/(3+5 sigma)+eps}) in the dyadic bookkeeping of explicit density-to-zero-free
  arguments (N7). Homogeneous-dynamics and Boolean-analysis scouts and the N(T) computation timed out.
-/

namespace Fayed.R279
end Fayed.R279
