/-!
# Round 443 — node N10 (explicit ψ, Mertens and Liouville bounds), round 14, check of a new Zenodo paper (Bhattacharjee–Bhattacharya, 'Liouville Energy and the Rightmost Zeta Zero', 24150355 v3, 12 Nov 2026): the energy exponent log Q_λ(N)/log N → 2Θ and the one-sequence criterion (Thm 3.1, Cor 4.4, via Pintz's single-pole kernel, Prop 4.1) and the Beurling-type prime-perturbation threshold β > Θ (Thm C) look new; the weighted L² Nyman–Beurling form (Thm 13.5) differs from our L^p discovery; |L(x)| ≤ x/2641 + 8.3√x is Cohen–Dress–El Marraki through λ = μ ∗ 1_□ and holds to 3·10^7; the Lean project proves only finite identities (energy_nonneg states less than its table claims); our Mertens record gives |L(x)| ≤ 1.74·10^5 x^0.9382 (routine); random projections, multiple testing and symplectic integrators transfer no tool

* Node: N10 (14th round; R437 before). Round type: [check] a new unrefereed paper with a Lean project (Zenodo sweep, 12 Nov; no arXiv
  listing: 11 Nov 2026 is an arXiv holiday, next listing Fri 13 Nov).
* [check] Thm 3.1, Cor 3.2, Thm 4.2, Cor 4.3–4.4 (paper.txt l. 418–730): Q_λ(N) = N⁻¹Σ_t(Σ_{t<n≤t+N} λ(n)1_{[1,N]}(n))², V(N) = Σ_{j≤N} L(j)²;
  log Q/log N → 2Θ, and Q(N_k) ≤ N_k^{2σ} along one sequence (σ ≥ 1/2) implies Θ ≤ σ; the all-large-N lower bound uses Pintz's
  single-pole kernel (Prop 4.1) and V(N) ≤ (2N+1)Q(N): NEW as stated (no prior matching theorem found), routine in method.
* [falsify] |L(x)| ≤ x/2641 + 8.3√x (Thm 15.1, App. I.1): derived from Cohen–Dress–El Marraki's |M(x)| ≤ x/4345 via λ = μ ∗ 1_□
  (ζ(2)/4345 < 1/2641) plus 8.29293√x for small terms; checked for every x ≤ 3·10^7: largest ratio |L(x)|/bound = 0.158 at
  x = 96,862 (L = −414): HOLDS, routine. Our Mertens record gives a power saving: |L(x)| ≤ 1.74·10^5 x^0.9382 and |m(x)| ≤
  1.68·10^6 x^−0.0618 for all x ≥ 1 (COROLLARIES, updated from Mertens78's constants).
* [check] Lean (LiouvilleEnergy/, scripts/Axioms.lean): sorry- and axiom-free, 10 + 11 finite theorems (Fejér/Toeplitz identities,
  ζ(2) constant, a log threshold); Thms 3.1, 4.2, 5.6 and the L bound are PDF-only; energy_nonneg proves 0 ≤ energy while its
  table also claims the equality case (minor FLAW of description).
* [check pre-emption] Thm 13.5 (l. 3196–3203, 3338–3350): d_{M,σ} → 0 ⟺ ζ ≠ 0 on Re s > σ in the weighted L² norm
  ∫|u|²x^{2σ−1}dx, called by the authors the Hilbert counterpart of Beurling's L^p scale; it does not imply our natural-number L^p
  criterion (discoveries/nyman-beurling-lp.md): DIFFERENT.
* [check] Thm C (l. 76) and Prop 5.7 (l. 1090): changing λ on a prime set with π_S(x) ≍ x^β breaks the exponent 2Θ when β > Θ and
  x^{1/2+o(1)} changes preserve it; liminf_{u↓0} u^{2m−1}J(β + u/2) ≥ |a_ρ|² per pole; no flaw; explore33's planted Beurling zero
  is a different model.
* Exploration (tool transfers, tested): Johnson–Lindenstrauss (preserves distances of sampled zero vectors, controls neither Re ρ
  nor d_N), Simes/Benjamini–Hochberg (zero deviations are no calibrated p-values), Störmer–Verlet backward error analysis
  (stability thresholds depend on the step, not on ζ): no transfer.
-/

namespace Fayed.R443
end Fayed.R443
