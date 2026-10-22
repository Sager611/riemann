DISCOVERY: Sato–Tate in every ideal class, for almost all discriminants. Fix a newform f without CM (e.g. f = f_E, E/Q without CM) and an interval J ⊂ [−1,1] with nonempty interior. For every A > A(J) := 𝒜(M*(J)) — M*(J) the least degree of a polynomial minorant of 1_J with positive Sato–Tate mean, 𝒜(M) = 3, 10/3, 11/3, (M+3)/3 for M = 2, 3, 4–8, ≥ 8 (𝒜(1) = 9/4), so A(J) ≤ max{11/3, 25.2/L + 1.7} with L the length of the θ-interval — all but O(X^{1−c}) fundamental D with |D| ≤ X, negative or positive, have in every (narrow) ideal class of Q(√D) and in every (x, x + x|D|^{−κ}] with x ≥ |D|^A degree-one primes 𝔭 with cos θ_p ∈ J; and inside every class the angles θ_p of these primes are Sato–Tate distributed up to discrepancy 10/(3A−3) (A > 11/3, smooth weights). With published inputs only (without the Pascadi–Thorner 2025 preprint) 11/3 becomes 5.

# ChebST report (2026-10-21)
Proofs: [NOTE.md](NOTE.md). Novelty: [NOVELTY.md](NOVELTY.md) (first line NEW). Brief: `.scratch/orca/chebst-brief.md`. Paper proofs, not refereed, no Lean.

## Results (complete proofs in the NOTE)
- **Minorants (§1).** Chebyshev–Markov–Stieltjes polynomials for the Sato–Tate measure, written out: for every interval J and n ≥ 1, P⁻ ≤ 1_J ≤ P⁺ of degree 2n with means μ_ST(J) ∓ (λ_n(α) + λ_n(β)), λ_n(cos θ) = 4sin³θ/((2n+3)sin θ − sin(2n+3)θ) ≤ 2/(n+1). Positive mean once n ≥ ⌈37.7/L⌉, so M*(J) ≤ 75.4/L + 2; asymptotically degree ≈ 4π/L (centred) to 6π/L (at an end); quadrature certificates show M* ≳ 2π/L, so A(J) ≍ 1/L is the true order.
- **Proposition T_m for every m (§2).** The twisted prime sums over every class are ≤ η₀ν_DH/h for almost all D once A > 𝒜(m). Route F (A > (m+7)/3) is RSFam's argument written for general m; routes H, H′ are new.
- **Propositions F′, F″ (§3, new).** Family zero density for Sym^m f ⊗ θ_χ over all |D| ~ X: X^{1+3(1−σ)} + X^{3/4+(11/2)(1−σ)} + X^{(m+6)(1−σ)/(3−2σ)} (RSFam: X^{(m+10)(1−σ)/(3−2σ)}). Per-D mollifiers make class II cost the family size X^{3/2}; class I uses the per-D sieve, Pascadi–Thorner's Hölder sieve and Humphries–Thorner, with D-independent coefficients from a squarefree detector (RCApps3 §2, no pole for m ≥ 1).
- **Theorem 1 (§4):** every interval, explicit A(J); every primitive form of discriminant D represents primes p ≤ 2|D|^A with a_E(p)/(2√p) ∈ J; p = a² ± nb² for almost all n. Examples: θ_p ∈ [π/4, π/2]: A > 5; |θ_p − π/2| ≤ 0.05: A > 127/3.
- **Theorem 2 (§5):** discrepancy ≤ λ_n(α) + λ_n(β) + 1/(n+1) (≤ 5/(n+1); half-lines 3/(n+1)) for A > 𝒜(2n); Corollary 2.1: genuine equidistribution with A = A(|D|) → ∞, o(X) exceptions, no explicit rate.
- **Theorem 3 (§6):** D > 0 (narrow classes) is included above; sectors for D < 0 (Hecke angle in arcs of length ≥ |D|^{−κ}) with the same thresholds (8/3 for degree-one minorants).

## Thresholds (minorant degree M)
| M | 1 | 2 | 3 | 4 | 5–8 | ≥ 8 |
|---|---|---|---|---|---|---|
| this note | 9/4 | 3 | 10/3 | 11/3 | 11/3 | (M+3)/3 |
| RSFam D_M | 9/4 | 3 | 10/3 | 11/3 | (M+7)/3 | (M+7)/3 |

## Uniformity in m (NOTE §7)
Explicit in m: conductors, coefficient bounds, the per-D large sieve, mollifier lengths, Thorner's zero-free region (IMRN 2025, arXiv 2404.06482, eq. (5.1): no Landau–Siegel zeros for Sym^m f ⊗ χ, every Dirichlet χ, every m; re-read in the TeX 2026-10-21, [FIX.md](FIX.md)), c_m = 10⁷(m+1)⁴. Not explicit: Brumley–Thorner–Zaman's constant B_m (genus twists; enters δ logarithmically), Humphries–Thorner's constant for π′ = Sym^m f (route F only; routes H, H′ use only GL₂ family sieves, whose constants do not depend on m), Lemma B's C(ε, m). Hence no explicit rate when M grows with X.

## Relation to RSFam and significance (for the reviewer)
RSFam §7 already stated Theorem D_M (any degree-M minorant, A > (M+7)/3) with a short proof; its DISCOVERIES line covers only degree ≤ 2. This note proves Proposition T_m for every m with all m-dependence, adds explicit minorants for every interval, the equidistribution form, sectors, and lowers the thresholds for M ≥ 5 by up to 4/3 (new family density). The step from T_m to intervals and equidistribution is classical; the reviewer may prefer to merge this into the RSFam entry rather than list it separately.

## Limits (NOGO with present tools)
- Slope 1/3 per degree: the second moment Σ_χ|L(½, Sym^m f ⊗ θ_χ)|² ≍ |D|^{(m+1)/2} sits at the large-sieve diagonal; Lindelöf on average would remove it.
- The constant 11/3 (5 without PT): class I lengths X^{3/2}…X^{7/2} below the family threshold; a family sieve "beyond the conductor" (RCApps3 §6) would remove it.
- Theorem C (per-D density, all m, already in RSFam) never helps for m ≥ 2. Per-D (every D) statements need a log-free per-D density.

## Novelty
arXiv (34 phrase queries; 4 bulk listings, 378 records scanned), the 400 arXiv ids 2610.00001–00400, Hexagon (all 342 works), Zenodo (20 queries), OpenAlex (8), Semantic Scholar, openai/math (749), TeX of Humphries–Thorner, Pascadi–Thorner, Thorner 2021 (arXiv 2002.10450) and Thorner 2025 (arXiv 2404.06482, the source of (J5)). Nothing states a Sato–Tate condition on primes in ideal classes of varying quadratic fields with any power of |D|. GXL evidence API unavailable (HTTP 429).
