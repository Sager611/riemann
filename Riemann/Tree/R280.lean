/-!
# Round 280 — node N9 (moments, Lindelöf, subconvexity): 13/84 still the record; Lindelöf above 7/8 recorded as a corollary

* Node: N9 moments, Lindelöf, subconvexity.
* State of the art: mu(1/2) <= 13/84 (Bourgain) is still the best pointwise bound (2026 math.NT listings Jan-Nov scanned; 2602.05731
  says its own "Lindelöf proof" is not genuine). Moments: sixth moment << T^{5/4 - 1/60 + eps} and twelfth << T^{2 -
  520/445121 + eps} (de Faveri-Pandey, arXiv 2610.01035). TREE.md N9 updated.
* Adversarial: Chang-Wang 2611.00924 Prop. gain: the algebraic numbers reproduce (root r = 0.45057074672406, B* =
  0.87446830769353, endpoint B* - 1/160000 = 0.87446205769353, reserve -197/20000000): consistent; full audit by Orca cwaudit.
* Corollary (orchestrator; a scout's objection was wrong): Titchmarsh's section 14.2 argument needs only zero-freeness, so the
  2026 half-plane gives mu(sigma) = 0 for sigma > 7/8 and, by convexity with 13/84, mu(sigma) <= (13/84)(7/8 - sigma)/(3/8)
  on [1/2, 7/8] (mu(3/4) <= 13/252). Routine: logged in COROLLARIES.md.
* Exploration ideas: the uniform GL(2) bound C(pi)^{29/120 + eps} (2607.04476) in amplified fourth moments for Hecke families
  (N9/N11); Wooley's moment bounds over Z[t]/t^2 (2611.03487) in Vinogradov-Korobov mean values (N7/N9). Two scouts timed out.
-/

namespace Fayed.R280
end Fayed.R280
