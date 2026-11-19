/-!
# Round 485 — node N2 (de Bruijn–Newman), round 14, [check]/[falsify] of Moki & Julio's degree-2 bounds (Zenodo 24318337) and Missiaen v15.3.0: 0.449938 and 0.489938 are superseded by the Lean-checked 31/32 strip (Λ_f ≤ 225/512 for every weight-one newform, routine); their product lemma (Lieb–Sokal) holds and with Freire et al.'s per-character pairs gives Λ_{Q(√15)} ≤ 0.2278 (conditional); Λ_DH bounds unchanged; LGV, Budan–Fourier and Grace–Walsh–Szegő transfer no tool

* Node: N2 (14th round; R467 before). Round type: [check]/[falsify] of fresh unrefereed Zenodo claims next to the node line (lessons 147, 201, 283).
  Moki & Julio, 'Explicit de Bruijn–Newman bounds below 1/2 in degree 2, and a product lemma for Dedekind zeta functions' (Zenodo 24318337
  v1.0.1, 17 Nov 2026; v1.0 = 23794434, 30 Oct), Theorem 1: Λ ≤ 0.449938 for L(s, η(z)η(23z)) and 0.489938 for L(s, η(2z)η(22z)), from a
  certified winding region of H_{0.38} (y ∈ [0.374, 0.5], x ≤ 12,500; 13.4 h on 9 cores), a tail estimate and de Bruijn's theorem twice.
* [check] (orchestrator) L(s, η(z)η(23z)) = ζ_K/ζ for K = Q(x³ − x − 1) (d_K = −23) and L(s, η(2z)η(22z)) = ζ_K/ζ for K = Q(x³ − x² + x + 1)
  (d_K = −44): a_p = r_p − 1 at all 429 and 428 good primes p < 3000 (q-expansions by Euler's pentagonal theorem). So L(s, f) ≠ 0 on
  Re s > 31/32 (Cureton–Thorner, Lean-checked, every ζ_K) and, in their normalisation (H_0(z) = ξ((1 + iz)/2); their Step 1: Δ = 1 gives
  √(1 − 0.76) at t = 0.38), Λ ≤ (15/16)²/2 = 225/512 = 0.439453 for both: Theorem 1 is superseded and its Steps 2–4 are redundant
  (the 31/32 strip confines the zeros of H_{0.38} to |Im z| ≤ 0.34483 < 0.374). Every weight-one newform: Brauer induction (COROLLARIES).
* [check] Dobner's class (arXiv 2005.05142 TeX l.194–215, 249–254, 317–325): weight-one L(s, f) ∈ S♯ with Γ_C, root number absorbed by a
  phase; a zero-free half-plane Re s > θ gives Λ_F ≤ 2(θ − 1/2)²: HOLDS (scout artin).
* [falsify] Theorem P (Λ(Φ ∗ Ψ) ≤ max(Λ(Φ), Λ(Ψ), 0)) in its finite-degree form Λ(PQ) ≤ max(Λ(P), Λ(Q)), Λ(P) = min{t : e^{−tD²}P real-
  rooted}: exact Λ((z² + 1)²) = 1/3 and Λ((z² + 1)(z² + 4)) = 5/12 + √46/24 ≤ 2; 400 random pairs of degrees 2–5: max excess 0 (ratio min
  0.202, median 0.636): HOLDS (scout thmp). Proof (NOTE.md l.201–256): exp(−2t∂₁∂₂) with a = 2t ≥ 0 is the stability-preserving sign,
  class-D growth gives order < 2, the diagonal restriction keeps real-rootedness: HOLDS (scout thmp2).
* [check] prior art: no earlier explicit bound below 1/2 for a degree-2 or Dedekind Λ on arXiv (Dobner; Stopple 1301.3158, a lower bound;
  Andrade–Chang–Miller 1310.3477, existence) and no earlier statement of Theorem P (scout prior); our Λ_K ≤ 225/512 for every number
  field (R357, 5 Nov) postdates their v1.0 (30 Oct).
* [check] Freire et al., arXiv 2610.27343, Table certified_pairs_m400 (D = 1 reference pairs (t₀, y₀), λ_pair = t₀ + y₀²/2 rounded up;
  H_0(z) = ξ((1 + iz)/2), flow e^{tu²}): 131 of the 242 real primitive characters of conductor 3–400 have λ_pair < 0.2787756662, minimum
  0.2278 at d = 60 (t₀ = 43/200, y₀ = 4/25), maximum 0.29782 (their 953/3200); the N2 line had only the uniform 953/3200 (scout freire;
  orchestrator re-extraction identical). With Theorem P, Λ_{Q(√d)} ≤ max(Λ_ζ, λ_pair(d)): best Λ_{Q(√15)} ≤ 0.2278 (COROLLARIES,
  conditional on the unrefereed table).
* [check] Missiaen v15.3.0 (Zenodo 24242050, 15 Nov 2026): 28825920237/125000000000 = 0.230607361896 < Λ_DH < 2(σ* − 1/2)², the N2
  line's v15.2.0 values; the verification suite reruns: 'ALL PAPER XV CERTIFICATES VERIFIED 100% IN 91.14s' (memcap peak 29 MB,
  orchestrator); the new Pólya–Schur multiplier statements do not touch the node line.
* Exploration: Lindström–Gessel–Viennot (no planar-network model for Li coefficients; det(λ_{i+j+1})_{i,j<2} = −0.00373), Budan–Fourier
  (certifies finite truncations only), Grace–Walsh–Szegő (no stable multiaffine lift of Ξ; degree 2 is a routine check): no transfer.
* Scouts: 8, all answered.
-/

namespace Fayed.R485
end Fayed.R485
