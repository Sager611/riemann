/-!
# Round 483 — node N6 (Jensen-polynomial hyperbolicity), round 13, [check]/[falsify] of Spark's Jensen trace matrices (Zenodo 24242312): the sharp quartet criterion sin²(2 arg ρ) ≤ (4d − 6)/(d(d − 1)) holds (discriminant, Borcea–Brändén, 1,800 random tests) and gives J^{d,0} hyperbolic for d ≤ 6.398·10^25 at every t ≥ 0, 4× our sector-lemma range but at shift 0 only; with the Lean θ₁₄₁ strip d ≤ 6.4568·10^25 (routine); optimal transport, Riesz energy and Meyer sets transfer no tool

* Node: N6 (13th round; R477 before). Round type: [check]/[falsify] of a fresh unrefereed Zenodo claim beyond our node record (lesson 147):
  Spark, 'Jensen Trace Matrices in the de Bruijn–Newman Flow: Certified Bounds and Last-Contact Asymptotics' (Zenodo 24242312, 15 Nov 2026),
  Theorem ba-large-finite-degree: the fixed-node Hermite trace matrices Q_{n,t} of the degree-2n Jensen polynomials of Ξ_t are PSD with
  nullity one for 1 ≤ n ≤ 3199·10^22 and every t ≥ 0, i.e. J^{d,0}(Ξ_t) hyperbolic for d ≤ 6.398·10^25, from RH to 3·10^12 and the 7/8 strip;
  our N6 record is d ≤ 1.6144·10^25 for every shift (R463, R466).
* [falsify] Proposition quartet-preserver (TeX l.1800–1860): q_ρ(D)x^d = x^{d−4}(x⁴ − 2(d)₂Re(ρ^{−2})x² + (d)₄|ρ|^{−4}) is real-rooted iff
  sin²(2 arg ρ) ≤ (4d − 6)/(d(d − 1)) (the discriminant, re-derived); preservation on all real polynomials of degree ≤ d is Borcea–Brändén's
  finite-degree symbol test (q(D)(z + w)^d = (q(D)x^d)(z + w)); numerically the critical angle matches for d = 4–30 and 1,800 random
  real-rooted inputs (double roots, clusters) at 0.999 of it stay real-rooted: HOLDS (orchestrator; the scout returned nothing).
* [check] transfer (Proposition ba-wall-transfer, l.1925–2000): Q_{n,t} controls the unshifted J^{d,0} only (d = 2n); θ_j = arg z_j² = 2 arg z_j;
  the walls give B_j/A_j < 1/55000 and the reserved real factors keep strictness: HOLDS, shift 0 only (scout transfer).
* [check] the finite-degree proof (l.2373–2445): (1 + 1/15000000)^{10^9} > 2^{60} > 32U (log gap 25.08), the moment bound
  (log ratio −0.108), 2(3199·10^22)ϑ² = 719775095970003199/180000000000000000 < 4 − 2/99999, 4/ϑ² = 6.3999991·10^25: HOLDS (scout finite).
* [check] our R463 route (Kim–Lee Thm 3 / Chasse sector: dδ² ≤ 1, δ = 2a√(R + a²)/R) loses the asymptotic factor 4 against the quartet
  criterion; the quartet transfer does not directly cover every shift (scout ours463: IMPROVABLE at shift 0 only).
* Consequence (routine, COROLLARIES): with the Lean θ₁₄₁ strip, every zero above U has θ_j < 2a/U, a = θ₁₄₁ − 1/2, so J^{d,0}(Ξ) is
  hyperbolic for d ≤ (U/a)² − 1 ≈ 6.4568·10^25; Spark's 6.398·10^25 holds for every t ≥ 0; all shifts stay at 1.6144·10^25.
  Prior art of the quartet criterion not settled (scout prior6 timed out).
* Exploration: optimal transport (no canonical zero measure without explicit-formula density input; R459 had Kantorovich–Rubinstein),
  minimal Riesz energy (zero triples are not s = 1 energy minimisers: 3.62466 > 3.46410), Meyer sets (finite gaps cannot prove uniform
  discreteness of Z − Z; R2.6 had Fourier quasicrystals): no transfer.
-/

namespace Fayed.R483
end Fayed.R483
