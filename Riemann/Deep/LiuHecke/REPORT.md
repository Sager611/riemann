NO REFUTATION: in Liu, arXiv 2610.26835, the step from cyclotomic fields to every number field (Part III: run Part II over any F ⊇ μ12, then descend through K(ζ12)) is sound as written and adds no unproved step. What is still unchecked is Part II's own analytic derivation, which uses nothing about the field beyond μ12 ⊂ F.

# Audit of the general-field step in Liu's Hecke half-plane paper (2026-10-21)
Source: arXiv e-print TeX `number_fields_v2.tex` (11,410 lines). Details and the full ingredient table: [NOTE.md](NOTE.md).

## The claim that Part II "applies verbatim" (Thm vi:prop:transfer, TeX l. 10740–10862) is literally true
- In Part II (l. 1585–10423), "cyclotomic" and Q(ζ_m) occur only in theorem statements and in "fix a cyclotomic F" set-up lines. Searches for Galois, automorphism, Frobenius and embedding find only local uses. The only cyclotomic-specific steps, the Part II descent and the class 𝒞, are replaced in Part III.
- Every proof works over one fixed F with S-arithmetic: S-integers R (S chosen so R has trivial class group), allowed generators via Hilbert reciprocity, Tate–Weil Poisson summation on R ⊂ F_S, and Hecke functional equations.
- Each input that depends on the field exists for every F ⊇ μ12:
  - Theta function: Kazhdan–Patterson (with the 1985 correction) works over any global field containing μ_n. Tame uniqueness is KP Cor. I.3.6. Part II keeps the S-part as a finite sum of local products, so uniqueness at places above 3 is never needed.
  - Gauss-sum corrections Λ = Λ_{Q(ω)}∘N_{F/Q(ω)} and G_0: Davenport–Hasse prime by prime plus Weil's index product formula. F/Q(ω) need not be normal.
  - Quadratic large sieve: Goldmakher–Louvel Thm 1.1 (any number field). Cubic large sieve: Blomer–Goldmakher–Louvel Thm 1.3 (any K ⊇ μ_n, n ≥ 3), with bound (MN)^ε(M+N+(MN)^{2/3}). Both statements were checked in the arXiv TeX.
  - OpenAI's primitive-roots preprint (§§3–8, App. A) is written for "a fixed cyclotomic field" but uses only μ12 ⊂ F. Cureton–Thorner (Hexagon 2610.00200v2, Prop. "field-scope") re-checked it for every F ⊇ μ12 and formalized the KP cubic constructions over such F in Lean.
- The degree of F enters only through constants, log powers (1+log Z)^{r_F+|S_f|−1}, height orders and AFE lengths U^{1/2}H^{r_F}:
  - H ≤ Z^τ with τ < d_min/10⁴.
  - Lengths above 1/2 use M2 without prime sums, and the coefficient D^{−θ} ≤ U^{−θ/2} covers them.
  - Exponents and mesh are fixed before F (Lemma v:lem:arithmetic-independent-mesh).
  All of this was checked, including the balance identity −θm − gv(1−2m)/2 = −θ/2.

## Descent (Prop vi:prop:descent, l. 10875): correct
- χ of order d, M/K its cyclic field, E = M(ζ12): ζ_E = ∏ L_F(s,ψ) over ψ ∈ Gal(E/F)^ (Gal(E/F) embeds in Gal(M/K)); ζ_E = ζ_M·∏_{ω≠1} L_M(s,ω) with entire factors; ζ_M = ∏_j L_K(s,χ^j).
- A ratio of residues settles s = 1. Ramified characters are included, also those ramified at real places.
- There is no height threshold on zeros. Constants that depend on K enter only the Z-thresholds of the Mellin argument, and its Gaussian weight excludes zeros at every height. Hu (Lemma "descent") and Cureton–Thorner (Lean) have the same lemma.

## Ingredient status
- (b) Already in the literature for any F ⊇ μ12: generators and reciprocity, Gauss corrections, theta datum, both large sieves, Hecke functional equation and AFE, prime ideal theorem in ray classes.
- (a) Proved in Part II, with nothing specific to F: orbit Poisson summation and its kernels, the radial functional equation, the inverse moment (which gives the 11/12 seed), the product moment M2 (13/18 ≤ κ ≤ 1, including arbitrary lengths without prime sums), the detector and row count, the averaged numerator, and the scalar system.
- (c) Missing or false: none found.

## Comparison
- Cureton–Thorner get 31/32 for every K, Lean-checked, by the same route. They avoided the 7/8 moments as "further requirements".
- Hu (Hexagon 2610.00216, 13 Oct) gets 7/8 for every unitary Hecke character of every number field by the same route. He transfers OpenAI §§17–18 by citing their locations, and his list of field-sensitive inputs (App. A) matches Liu's. His paper is on Hexagon, not arXiv, and is not limited to special fields.
- Rubinstein–Restrepo Borrero (Hexagon 2610.00220) treat only Q(√−3), via the CM curves y² = x³ + d.

## Remaining risk, outside the field step
- Not machine-checked anywhere: the general-field moments M1, M2 and row count; M2 for κ < 3/4, needed since κ♭ = 0.74977 (Liu 2610.12234 Prop 28.1 goes below OpenAI's printed 3/4 ≤ κ ≤ 1); and the Xu–Sang averaged numerator, even over Q(√−3).
- Only OpenAI's 7/8 (ζ and Dirichlet L-functions) is Lean-checked; our port builds with no sorry. Part II's ~5,000 lines of moment derivations were not checked line by line here.
