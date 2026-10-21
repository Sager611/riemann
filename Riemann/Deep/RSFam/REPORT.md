DISCOVERY: for every A > 3, all but O(X^{1−c}) fundamental discriminants D with |D| ≤ X, negative or positive, have in every ideal class (narrow if D > 0) primes p with a_E(p) > 0 and primes with a_E(p) < 0 in every interval (x, x + x|D|^{−κ}], x ≥ |D|^A (E/Q any fixed curve without CM; any newform f in place of E). So every primitive binary quadratic form of discriminant D represents primes p, p′ ≤ 2|D|^A with a_E(p) > 0 > a_E(p′). One-sided conditions a_E(p) > −α√p (α > 0) need only A > 9/4. No result of this kind with any power of |D| was found, per D or on average.

# RSFam report (2026-10-21)
Proofs: [NOTE.md](NOTE.md). Brief: `.scratch/orca/rsfam-brief.md`. Scratch: `.scratch/rsfam/`. Paper proofs, not refereed, no Lean, no 2026 half-planes.

## Results (complete proofs in the NOTE)
- **(a) Lemma A, (b) Lemma B (routine).** The project's weighted large sieve over primitive ideals carries the Rankin–Selberg coefficients of Sym^m f ⊗ θ_χ (Deligne bounds); the second moment Σ_χ |L(½+iv, Sym^m f ⊗ θ_χ)·M_V|² ≪ |D|^{(m+1)/2+ε}V^{1+ε} saves the family size h over convexity.
- **(c) Theorem C (new statement).** Σ over the class group characters χ of D of the zeros of L(s, Sym^m f ⊗ θ_χ) with β ≥ σ, |γ − t| ≤ 1: ≪ (|D|^{m+2}T^{2m+2})^{(1−σ)/(3−2σ)+ε}. For f itself (degree 4, conductor ≍ |D|², about |D|^{1/2} members): |D|^{3(1−σ)/(3−2σ)}, nontrivial for σ > 3/4. Best published (Pascadi–Thorner 2025, Thm 1.3): |D|^{(m+7)(1−σ)/(3−2σ)}, nontrivial only for σ > 13/14 (m = 1).
- **(c′) Proposition F (new statement).** Over all |D| ≍ X: (X^{m+10}T^{2m+10})^{(1−σ)/(3−2σ)+ε}, against X^{m+11} from Pascadi–Thorner; the gain is Lemma B in the Humphries–Thorner detector.
- **(d) Theorem D, Corollary D1, Theorem D⁺.** For a set S ⊂ [−2,2] of Hecke eigenvalues with a Chebyshev minorant P = c₀ + c₁x + c₂(x²−1) ≤ 1_S, c₀ > 0: for A > 3 (A > 9/4 if c₂ = 0), almost all D have ≥ c₀ν_D H/(32h) · (log weight) primes p ∤ ND with λ_f(p) ∈ S in each class and each (x, x+H], x ≥ |D|^A, H ≥ x|D|^{−κ}. Covers λ_f(p) > 0, < 0, ≠ 0, |λ_f(p)| < β (β > 1), |λ_f(p)| > β (β < 1), λ_f(p) > α (α < ½); and primes a² + nb², a² − nb² with a_E(p) > 0 for almost all n.

## How it works
- Detector: λ_f(p) > 0 ⇐ P(λ) = λ(λ+2)/8 > 0, i.e. the twists by 1, f and Sym²f (degree 2, 4, 6) for all class group characters at once.
- Non-genus characters: zeros with β < 5/6 by the trivial count per D; the rest by Proposition F and Markov over D. Threshold: 6/(A+1) − ½ < 1, i.e. A > 3.
- Genus characters give twists Sym^m f ⊗ χ_d: large d by the Brumley–Thorner–Zaman log-free density (Markov over D); small d by the same estimate with its Landau–Siegel repulsion, plus no Landau–Siegel zeros on GL₂/GL₃ (Hoffstein–Ramakrishnan, Banks). The m = 0 part is RCApps2 (A > 3/2).
- Constants c, κ effective but not numerical (BTZ and Siegel-zero constants are not explicit); thresholds exact.

## Thresholds (minorant degree M)
| M | condition | this note | from published densities |
|---|---|---|---|
| 0 | none | 3/2 (RCApps2) | 8/3 |
| 1 | a_E(p) > −α√p | 9/4 | 3 |
| 2 | sign, a_E(p) ≠ 0, … | 3 | 10/3 |
| ≥ 3 | small Sato–Tate intervals | (M+7)/3, conditional on no Siegel zeros of Sym^m f ⊗ χ_d | (M+8)/3 |

## Limits (NOGO with present tools)
- Sign needs Sym²f: no degree-one minorant has positive mean. At A = 3, σ = 5/6 the trivial per-D count, Theorem C, Proposition F and the Markov bar all equal X^{3/2}.
- Below 3: needs a family large sieve below X^{7/2} or Lindelöf on average for the degree-6 family (A₀ = 5/3, resp. 1). Michel–Venkatesh subconvexity, Hölder with fourth moments and per-D sieves inside the family detector give nothing; a Lindelöf second moment for f ⊗ θ_χ is an averaged mixing conjecture.

## Novelty check
35 arXiv API queries and 4 evidence-API claims (`.scratch/rsfam/axq*.out`, `gxl.out`); TeX of Humphries–Thorner, Brumley–Thorner–Zaman, Lemke Oliver–Thorner, Pascadi–Thorner read. Closest: Pascadi–Thorner (general densities, weaker here), Lemke Oliver–Thorner and Thorner 2021 (Sato–Tate least primes in the level aspect), Nordentoft, Van Order, Burungale–Tian (central values over class group characters), Vaishya et al. (signs of λ_f(n) at integers represented by one fixed form). Nothing on zero density for this family or on Hecke-eigenvalue conditions for primes in classes of growing class groups.

## Significance (for the reviewer)
A new theorem: the first polynomial bound in |D| for primes in ideal classes with a Sato–Tate condition, for almost all D. The method is RCApps2's (per-D large sieve + Humphries–Thorner detector + Markov over D) extended to Rankin–Selberg twists; new ingredients are the twisted second moment and the genus-twist analysis. Theorem C and Proposition F are new statements by a classical method (could go to COROLLARIES.md if judged routine).
