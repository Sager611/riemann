/-!
# Round 158 — record re-checks; four record assessments; new-submissions finder

* Weil positivity (our half-width 1.55): nothing larger. A scout reported "Liu, support 2.125"; that is the window LENGTH
  2.125 = 2·17/16 (Liu's half-width 17/16, already in our list); github.com/jose-paul-b/weil-positivity-certificates
  (2026-10-10) confirms only Zhu 0.8 and Mori ½log 15.
* NT entries (least primes in classes, ring-class density, CM analogue): nothing newer (math.NT pastweek/new, Zenodo).
* Gale (our 0.8906210556 in R³): arXiv 2610.28828 still v1; its repository (latest commits 2026-10-22) still has
  0.889281162393 for R³, and also lists verified examples in R⁴ (61 points, ratio 0.913148536990) and R⁵ (32 points,
  0.892304779206) → our refinement method may improve those (Orca gale45).
* Assessments (chance that a few CPU-hours of search improve the record): 2610.19980 stabbing unit squares 14/9 (2/5 by plain
  local search, 4/5 with staircase moves; conjectured (35 − √5)/20); 2609.21968 Moser's worm, convex universal cover of unit
  arcs, 0.239 ≤ α ≤ 0.246332237 (3/5, upper bound; code public) → Orca moser; 2610.06891 7-universal tournament of order 14,
  13 open (2/5); 2610.12122 Ramsey R(6,8) ≥ 135, R(8,10) ≥ 345 (2/5; data public).
* Finder: 2610.29291 ("a partition is Latin iff wide"): 138 partitions of size ≤ 10, no counterexample.
-/

namespace Fayed.R158
end Fayed.R158
