REPAIRED: no mathematical gap in ChebST or RSFam. (J5) is in the paper both notes cite, Thorner, IMRN 2025, no. 20, rnaf307 = arXiv 2404.06482v2 (eq. (5.1), i.e. Theorem 1.1 with n = 0; read in the TeX). The finding checked a different Thorner paper (Res. Math. Sci. 2021, arXiv 2002.10450). Every statement and constant stands; citations are now pinned; one missing reference (Proposition F″ in ChebST Theorem 3) is added.

# FIX: the (J5) citation (2026-10-21, brief `.scratch/orca/chebstfix-brief.md`)

## Source, read in the TeX of arXiv 2404.06482v2 ("Exceptional zeros of Rankin–Selberg L-functions and joint Sato–Tate distributions")
- Theorem 1.1 (thm:main, l. 97–118): F totally real, π, π′ ∈ RA_F (non-dihedral, regular algebraic), m, n ≥ 0, M = max{m, n} ≥ 1, m ≠ n: L(σ, Sym^m π × (Sym^n π′ ⊗ χ)) ≠ 0 for σ ≥ 1 − c/(M²[M^{max{y_π,y_π′}}log(MC(π)C(π′)) + log C(χ)]), c absolute and effective. Corollary 1.2 (l. 123–129): the same at height t. Remark (l. 130–132): for n = 0 put C(π′) = 1, y_π′ = 0.
- Proof, eq. (5.1) (eqn:SymL, l. 600–604), case n = 0, no ray-class restriction: L(σ, Sym^m(π) ⊗ ξ) ≠ 0 for σ ≥ 1 − c/(m^{2+y_π}log(mC(π)) + m² log C(ξ)), every ξ ∈ 𝔉₁. With ξ = χ∣·∣^{it}, C(ξ) ≤ (4/3)q_χ(3+∣t∣), this is (J5) as we state it: every m ≥ 1, uniform in χ, no exceptional zero.
- Checked for our twists (ω_f = 1, χ = χ_d ≠ 1): the auxiliary series 𝒟 (eqn:D_def) has its only pole at s = 1, of order 3 (L(s, Sym^m × Sym^m) and ζ²), non-negative coefficients at unramified places, and the factor L(s, Sym^m f ⊗ χ_d)⁴. So a real zero near 1 would be a zero of 𝒟 of order ≥ 4, more than its Goldfeld–Hoffstein–Lieman bound of 3.
- arXiv 2002.10450: its Cor. ZFR (self-dual π × π′) allows one real simple zero, and its region for Sym^m f has no twist. The finding is right about that paper, but neither note cites it for (J5). ChebST wrote "Thorner, IMRN 2025" without an arXiv number, probably the cause of the mix-up; the number is now given.

## Uses of (J5), and what each needs
| place | needs | met by |
|---|---|---|
| ChebST §2.1: c₁′ := c/(2m²), δ_e | the shape c/(m²[m^y log(mC(f)) + log(C(χ)(3+∣t∣))]), uniform in χ | (5.1); in (e) the bracket is ≤ 2L for ∣D∣ ≥ D₀(f, m, δ) |
| ChebST §2.4(e), end of §2, §6 (D > 0); RSFam §5.4(e), §5.6, §6, §7 (Theorem D_M) | genus twists L(s, Sym^m f ⊗ χ_d), ∣d∣ ≤ ∣D∣^δ, ∣γ∣ ≤ T₁: no zero with β > 1 − c₁′/L (L ≤ 1.3δ log∣D∣), including the real zero a self-dual quadratic twist could have | (5.1); for m ≤ 2 also Hoffstein–Ramakrishnan and Banks |
| RSFam §5.1: c₁′ := c₁/(m+3), m ≤ 2 | the Iwaniec–Kowalski region with no exceptional zero | Hoffstein–Ramakrishnan, Banks, and (5.1) |
| class-group twists χ ∈ Ĝ_c (RSFam §5.6, ChebST §2) | nothing: (J5) is not used | every zero with β ≥ σ₀, ∣γ∣ ≤ T_X (real ones near 1 included) is counted by Proposition F (RSFam) or F′/F″ (ChebST); for D ∉ E₁ none has β ≥ σ_J |
| ChebST Theorem 3 (sectors, ∞-type ≠ 0) | nothing: no genus factors | E₁′ near σ = 1 needs F″ (with F′ alone the term X^{1+3(1−σ)} exceeds the Markov room as σ → 1); F″ is now named, with the same changes as F′ |
| Verdicts, §7, §9 and REPORTs of both notes | bookkeeping | arXiv numbers added |

## Verdict
- DISCOVERIES "Primes with prescribed Frobenius-trace sign in every ideal class" (RSFam Theorem D, m ≤ 2): survives unchanged (A > 3; one-sided conditions A > 9/4; same constants).
- ChebST candidate (Theorems 1–3, thresholds 𝒜(M), Corollary 2.1): survives unchanged.
- (J5) is a real input: without it, one genus factor d whose twist has a real zero near 1 would spoil about X/∣d∣ discriminants, and the exceptional set could not be shown to be O(X^{1−c}). It holds, by (5.1).

## Edits (marked 2026-10-21 / FIX.md)
ChebST NOTE: header line; §2.1 source pinned and c₁′ checked; (e), end of §2, Verdicts, §7 table, §9: arXiv numbers; Theorem 3 (iii)(a): F″ for σ_j ≥ σ₂ (V := X^{7/2}M_X^5, class I by (LS2) only). ChebST REPORT: numbers. RSFam NOTE: (J5) pinned (TeX lines, (5.1), note on 2002.10450); §7 Theorem D_M; §10 note. RSFam REPORT: number.
