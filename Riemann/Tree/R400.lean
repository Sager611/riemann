/-!
# Round 400 — node N3 (Li/Keiper), round 9: records stand (Zenodo, HAL and listings since 6 Nov: nothing new; no newer Dirichlet GRH heights than Platt 2016 found); the limiting term of our 106.25 T^2 range is a per-zero worst case, so explicit zero density cannot raise it; the newest Zenodo Li-criterion item claims no proof (known); sandpile groups and n-body central configurations give no route, the Kac-Moody scout returned nothing

* Node: N3 Li/Keiper (ninth round; R383 before). Three scouts timed out (arXiv API throttling, a 56k-line Zenodo PDF); I ran
  the [check] and the GRH-height query myself.
* [check] lst.py 'Li coefficient|Keiper|Li criterion' (math.NT, math.CV, Oct-Nov) and Zenodo '"Li coefficients"': no title hit,
  newest Zenodo items are Sep 2026 or older working notes; HAL '"Li coefficients"': two hits (2013, 2019): HOLDS.
* Wild: our 106.25 T^2 (Li200k REPORT, Theorem and proof 1-2; verify.py l. 10-14, 126-145) bounds each zero above T by
  q_n = 4 sinh^2(n a/(2 gamma^2)) + 4 n a^2/gamma^3; the limit is one worst-case zero with beta near B, which a zero-density bound
  cannot exclude, so Kadiri-Lumley-Ng-type density does not raise the constant: HOLDS (NOGO).
* Adversarial: Zenodo RH claims via Li's criterion, Keiper-Li or Bombieri-Lagarias since 1 Nov: the newest, 23891944, states 'No
  proof of RH is claimed' and is already in NOGO: KNOWN.
* Input check (Dirichlet ranges, COROLLARIES line 20): an arXiv API query for newer GRH verifications of all Dirichlet
  L-functions found nothing beyond Platt (Math. Comp. 2016): HOLDS (weak: one query).
* Harvest: timed out on the throttled arXiv API: none.
* Exploration: abelian sandpile groups (the tree number is one symmetric function of the Laplacian spectrum, Kirchhoff and
  Northshield's Ihara residue, so it cannot certify that every eigenvalue lies in the Ramanujan interval), n-body central
  configurations (BKK and tropical counts need finite polynomial systems; KAM and Arnold diffusion need near-integrable
  Hamiltonian dynamics, unlike the dissipative de Bruijn flow): no route; affine Kac-Moody denominator identities: the scout
  returned no answer (untested, may be retried).
-/

namespace Fayed.R400
end Fayed.R400
