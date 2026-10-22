DISCOVERY: for every A > 3, all but O(X^{1−c}) fundamental discriminants D with |D| ≤ X, negative or positive, have in every ideal class (narrow if D > 0) primes p with a_E(p) > 0 and primes with a_E(p) < 0 in every interval (x, x + x|D|^{−κ}], x ≥ |D|^A (E/Q any fixed curve without CM; any newform f in place of E). So every primitive binary quadratic form of discriminant D represents primes p, p′ ≤ 2|D|^A with a_E(p) > 0 > a_E(p′). One-sided conditions a_E(p) > −α√p (α > 0) need only A > 9/4. No result of this kind with an explicit power of |D| was found, per D or on average (a per-D bound with an inexplicit power of order 10^10 looks derivable from published inputs but is stated nowhere: [NOVELTY.md](NOVELTY.md)).

# RSFam report (2026-10-21; repaired after review the same day)
Proofs: [NOTE.md](NOTE.md) (repairs listed in §10). Briefs: `.scratch/orca/rsfam-brief.md`, `.scratch/orca/rsfamfix-brief.md`. Scratch: `.scratch/rsfam/`, `.scratch/rsfamfix/`. Paper proofs, not refereed, no Lean, no 2026 half-planes.

## Results (complete proofs in the NOTE)
- **(a) Lemma A, (b) Lemma B (routine).** The weighted large sieve over primitive ideals carries the base-change coefficients λ_m(𝔞) of Sym^m f to K (written out prime by prime in (R1′); they are not λ_{Sym^m f}(N𝔞)), using only |λ_m(𝔞)| ≤ d_{m+1}(𝔞). Second moment Σ_χ |L(½+iv, Sym^m f ⊗ θ_χ)·M_V|² ≪ |D|^{(m+1)/2+ε}V^{1+ε}(|v|+3)^{m+1+ε}, uniformly in V and v; it saves the family size h over convexity.
- **(c) Theorem C (new statement).** Σ over the class group characters χ of D of the zeros of L(s, Sym^m f ⊗ θ_χ) with β ≥ σ in the unit window |γ − t| ≤ 1: ≪ (|D|^{m+2}T^{2m+2})^{(1−σ)/(3−2σ)+ε}; all |γ| ≤ T costs an extra factor T. For f (degree 4, conductor ≍ |D|², about |D|^{1/2} members): |D|^{3(1−σ)/(3−2σ)} per window, nontrivial for σ > 3/4. Best published (Pascadi–Thorner 2025, Thm 1.3): |D|^{(m+7)(1−σ)/(3−2σ)}, nontrivial only for σ > 13/14 (m = 1).
- **(c′) Proposition F (new statement).** Over all |D| ≍ X: (X^{m+10}T^{2m+10})^{(1−σ)/(3−2σ)+ε}, against X^{m+11} from Pascadi–Thorner; the gain is Lemma B in the Humphries–Thorner detector.
- **(d) Theorems D, D⁺, D_M, Corollary D1.** For a set S of Hecke eigenvalues with a Chebyshev minorant P = Σ_{m≤M} c_mU_m(x/2) ≤ 1_S, c₀ > 0: for A > (M+7)/3 (M ≥ 2; A > 9/4 for M = 1), almost all D have ≥ c₀ν_DH/(32h) (log weight) primes p ∤ ND with λ_f(p) ∈ S in each class and each (x, x+H], x ≥ |D|^A, H ≥ x|D|^{−κ}. Covers λ_f(p) > 0, < 0, ≠ 0, |λ_f(p)| < β (β > 1), |λ_f(p)| > β (β < 1), λ_f(p) > α (α < ½), short Sato–Tate intervals (M ≥ 3); primes a² + nb², a² − nb² with a_E(p) > 0 for almost all n.

## How it works
- Detector: λ_f(p) > 0 ⇐ P(λ) = λ(λ+2)/8 > 0, i.e. the twists by 1, f and Sym²f (degree 2, 4, 6) for all class group characters at once. Positivity comes only from the m = 0 lower bound hS^{(0)} ≥ ν_DH/8 together with m = 1, 2 errors ≤ η₀ν_DH, η₀ = c₀/(16(|c₁|+|c₂|)) (NOTE §5, table of constants).
- Non-genus characters: zeros with β < 5/6 by the trivial count per D; the rest, including any exceptional zero, by Proposition F and Markov over D: such zeros are removed by discarding Proposition F's exceptional set of D, never by a zero-free region (NOTE §5.6). Threshold: 6/(A+1) − ½ < 1, i.e. A > 3.
- Genus characters give Dirichlet twists Sym^m f ⊗ χ_d: large d by the Brumley–Thorner–Zaman log-free density (Markov over D); small d by the same estimate with its Landau–Siegel repulsion (needed because the m = 0 main term may shrink to ν_DH), plus no Landau–Siegel zeros (Hoffstein–Ramakrishnan, Banks; Thorner, IMRN 2025, for every m). The m = 0 part is RCApps2 (A > 3/2).
- Constants c, κ effective but not numerical (BTZ and Siegel-zero constants are not explicit); thresholds exact.

## Thresholds (minorant degree M)
| M | condition | this note | from published densities |
|---|---|---|---|
| 0 | none | 3/2 (RCApps2) | 8/3 |
| 1 | a_E(p) > −α√p | 9/4 | 3 |
| 2 | sign, a_E(p) ≠ 0, … | 3 | 10/3 |
| ≥ 3 | small Sato–Tate intervals | (M+7)/3, unconditional (Thorner, IMRN 2025, arXiv 2404.06482, eq. (5.1), removes the Siegel-zero condition of the first version) | (M+8)/3 |

## Repairs after review (NOTE §10)
1. Theorem C is a unit-window bound; (C′) and every use state the extra factor T (window count 3T₁ in §5.4(b)). 2. Lemma B: the old tail O(|D|^{−8}) was not uniform in V, v; the cutoff is now enlarged by (|D|V(|v|+3))^ε, statement unchanged, uses re-checked. 3. Base-change coefficients written out ((R1′)) with Rankin's identity for the brief's λ_f(N𝔞)χ(𝔞); the large sieve and Deligne bounds re-verified for them. 4. Non-genus exceptional zeros are removed through Proposition F's exceptional set of D (§5.6). 5. Positivity needs the m = 0 lower bound and small m = 1, 2 errors; constants checked (`params2.py`: 85 checks OK).

## Limits (NOGO with present tools)
- Sign needs Sym²f: no degree-one minorant has positive mean. At A = 3, σ = 5/6 the trivial per-D count, Theorem C, Proposition F and the Markov bar all equal X^{3/2}.
- Below 3: needs a family large sieve below X^{7/2} or Lindelöf on average for the degree-6 family (A₀ = 5/3, resp. 1). Subconvexity, Hölder with fourth moments and per-D sieves inside the family detector give nothing.

## Novelty check
First audit: 35 arXiv queries, 4 evidence-API claims, TeX of HT, BTZ, LOT, PT. Second audit ([NOVELTY.md](NOVELTY.md), first line NEW): 65 more arXiv queries, OpenAlex (searches and author lists of Thorner, Wong, Zaman, Lemke Oliver), Zenodo, Hexagon (323 works), openai/math (749 preprints), 5 evidence-API claims, TeX of Thorner 2025, Thorner–Zhao 2026, Akbary–Wong 2021. Nothing with an ideal-class or form condition plus a sign or Sato–Tate condition, uniform in D.

## Significance (for the reviewer)
A new theorem: the first explicit power of |D| (3; 9/4 one-sided) for primes in ideal classes with a Sato–Tate condition, for almost all D, in short intervals. A per-D bound with an astronomically large exponent seems assemblable from published pieces (Linnik's method with BTZ and recent exceptional-zero results) but is not stated. Method: RCApps2's per-D large sieve + Humphries–Thorner detector + Markov over D, extended to Rankin–Selberg twists. Theorem C and Proposition F are new statements by a classical method (could go to COROLLARIES.md if judged routine).
