NO RESULT: this attempt proves H for no θ ≥ 0 (its θ = 0 case is the open Hughes–Pearce-Crump conjecture Σ_γ ζ′(ρ)² ~ (T/2π)L³/6, arXiv 2509.07792); every route tried meets the same wall: up to approximate-functional-equation errors, S_B is a first moment at the zeros with a twist of length T^{1+2θ} ≥ T (Bui–Heath-Brown reach T^{1/2}), and the part beyond their reach (normalized coefficient ≥ 1/24 for ∫P ≥ 0) is a formally related χ²-dual sum that is not small in the data. No counterexample found: the recipe constant is re-derived and approximately supported by a new numerical check at zero index 10¹².

# hmean — proof attempt of the holomorphic mean value H (node N8, RH branch)

## Statement
- Target (explore55): H: S_B = Σ_{0<γ≤T} ζ′(ρ)²B(ρ)² = (TL³/2π)(u² − 1/12) + o(TL³), u = ½ + θ∫P, 0 < θ < ½; H at θ_H > ¼ would give RH ⇒ κ* > 19/27.
- Proved: nothing new about H. Exact identity (no RH): ζ′(ρ) = −χ(ρ)ζ′(1−ρ) at every zero, so S_B = −Σ_γ χ(ρ)|ζ′(ρ)|²B(ρ)² under RH.
- Found (recipe, leading order; numerically confirmed): modulus localization. With Ẽ_K(s) = Σ_{n≤K}(ℓ − log n)n^{-s}, ℓ = log(γ/2π), τ = γ/2π,
  F(τ^κ) := Σ_γ ζ′(ρ)Ẽ_{τ^κ}(ρ)B(ρ)² ~ (T/2π)ℓ³ f(κ), f(κ) = 1/6 + μ + μ² + μ(1−κ)² + (1−κ)³/3, μ = θ∫P (0 < κ < 1);
  f(1⁻) = u² − 1/12 (H); f(0⁺) = ½ + 2μ + μ² is the normalized coefficient of ℓΣζ′(ρ)B(ρ)². Dropping the moduli n > τ^κ shifts the coefficient:
  S_B − F(τ^κ) ~ −(T/2π)ℓ³[μ(1−κ)² + (1−κ)³/3].
- Obstruction (heuristic, located; not a non-provability claim):
  1. χ² form: on the left edge ζ′(s)² = χ(s)²[…]², stationary phase gives phases e(2√(abm/h)) with Λ(m); for m ≥ T^{1+ε} the explicit formula
     formally returns the χ²-dual part of the same zero sum (stationary-phase heuristic; absolute values over zeros give no saving there).
  2. χ¹ form: S_B ≈ Σ_γ ζ′(ρ)A(ρ), A = Ẽ_τB² of length T^{1+2θ} (numerical reading, C1): Bui–Heath-Brown's S₁ problem with y → T^{1+2θ}.
     Formal substitution makes their large-moduli bound (Lemma 2, y^{1/3}T^{5/6+ε}) and O(yT^{1/2+ε}) term T^{7/6+2θ/3+ε}, T^{3/2+2θ+ε} ≫ TL³,
     already at θ = 0 (their hypotheses were not checked for this twist): the large-moduli estimate is where the route stops.
  3. Not only an error term: BHB's reach for the twist Ẽ_KB² is K = τ^κ with κ < ½ − 2θ (θ < ¼; empty for θ ≥ ¼). The remainder is
     S_B − F(τ^{½−2θ}) ~ −(T/2π)ℓ³c, c = μ(½+2θ)² + (½+2θ)³/3 ≥ 1/24 for ∫P ≥ 0 (P = x: c is 25 % of u² − 1/12 at θ = 0, 75 % at θ = ⅛),
     and up to AFE errors S_B − F(K) = −Σ_γ χ(ρ)²ζ′(1−ρ)G_{τ/K}(1−ρ)B(ρ)² (G_M = Σ_{n≤M}(log n)n^{-w}): a χ²-object again.
  4. The χ² piece is not small: in the Riemann–Siegel split ζ′(ρ) = Ẽ_{√τ}(ρ) + χ(ρ)G_{√τ}(1−ρ), the pure χ² part Σχ(ρ)²G(1−ρ)²B(ρ)²
     is .313 ℓ² per zero at zero index 10¹² (θ = ½; S_B = .485), up from .160 and .178 at ℓ = 9.1, 12.1: a proof through this split needs more than an upper
     bound for it (control including its cancellation against the χ⁰ and χ¹ pieces).
  5. Avoiding H by a lower bound for Σ(Im ζ′(ρ)B(ρ))² through short Dirichlet polynomials gains nothing near BHB's optimum (projection
     coefficients ∝ P′(1)/θ − P(1) → 0; quadratic gain vs linear loss; heuristic).

## Proof
[DERIVATION.md](DERIVATION.md): A identities; B recipe main term of H (re-derived, sympy); C localization formula; D the obstruction; E numerics.

## Evidence (all runs single-threaded under memcap 290 MB)
- [logs/recipe_check.log](logs/recipe_check.log): sympy, ∂_α∂_β D|₀ = ℓ³/6 + ℓ²λ + ℓλ²; Φ entire; f(κ) as stated; S₁ recipe = BHB's (½ + θ∫P).
- [logs/kappa_low.log](logs/kappa_low.log), zeros 50001–10⁵ (ℓ_rms 9.11; 59 s, peak 70 MB): F(1) = exact S_B to 4 digits; tails F(κ) − F(1) within 0.007 of
  f(κ) − f(1) for κ ≥ 0.5 (θ = 0, ¼, ½).
- [logs/kappa_hi2.log](logs/kappa_hi2.log), zero indices 10¹²+1…+10⁴ (height 2.7·10¹¹, ℓ = 24.475; 151 s, peak 242 MB; first run [logs/kappa_hi.log](logs/kappa_hi.log), 145 s, 210 MB, same numbers):
  θ = ½: F − S_B at κ = .1…​.5: .443 .341 .244 .162 .100 vs .446 .331 .237 .162 .104; S_B − F at κ = .5…​.9 (χ²-tail): −.100 −.067 −.028 −.012 −.002
  vs −.104 −.061 −.032 −.013 −.003; θ = ¼ same quality (all κ-differences within .010); exact S_B/ℓ² = .4847 vs u² − 1/12 = .4792 (θ = ½,
  the boundary of H's interval, an extrapolation), .3161 vs .3073 (θ = ¼).
- [logs/kappa_low2.log](logs/kappa_low2.log) (57 s, 72 MB), [logs/kappa_hi3.log](logs/kappa_hi3.log) (160 s, 243 MB): same F as above plus the three-piece split
  (DERIVATION E3); P₀ + 2X₁ + X₂ = S_B to 4–5 digits.
- [logs/kappa_z6.log](logs/kappa_z6.log), zeros 1981053–2001052 (ℓ_rms 12.10; 320 s, peak 161 MB): tails within .009 of f(κ) − f(1) for κ ≥ .3 (θ = ¼, ½);
  pure χ² piece X₂ = .178 (θ = ½), between .160 (ℓ 9.1) and .313 (ℓ 24.5).
- Tiny runs: [logs/tiny.log](logs/tiny.log), [logs/kappa_hi_tiny.log](logs/kappa_hi_tiny.log), [logs/kappa_hi2_tiny.log](logs/kappa_hi2_tiny.log) (identity of D6 per sample to 4 digits at K = √τ).

## Checks
- Recipe vs literature: θ = 0 constant 1/6 = Hughes–Pearce-Crump's 1/Γ(k+2) at k = 2 (random-matrix + ratios, independent of explore55); the same bookkeeping
  reproduces BHB's S₁ and Fujii's Σζ′(ρ) ~ (T/4π)L².
- Falsification of the recipe steps (brief item 2): arithmetic factor is 1 at zero shifts (balanced ratio; its derivatives are O(1/L)); only the swaps
  α↔δ, β↔δ exist (one (1−s)-numerator); no extra main term found; sampled values of f(κ) are approximately supported at two (three) heights.
  H stands as a conjecture.
- Pipeline: F(κ = ½) direct and via the χ²-tail agree to 4 digits; Im parts ≤ 0.013 for θ = ¼, ½ (θ = 0 noisy, |Im| ≤ 0.1: heavy tails of ζ′).

## Review (reviewer subagents, formulas pasted inline)
- Main-term block (B, C): GAP; ROUTINE. Algebra 1–5, 7, 9 OK; the numbers "support the polynomial quantitatively, particularly the interior-cutoff tails".
  Gaps fixed in DERIVATION: ζ′(ρ) ≈ Ẽ_τ(ρ) is not justified by the sharp AFE (now: AFE of ℓζ + ζ′ with lengths K, τ/K, accurate for K ≍ √τ, checked
  to 4 digits; K = τ is a numerical reading only); f(κ) stated for 0 < κ < 1 (κ = 0 as the limit, n = 1 at full weight); keep ℓ = log(γ/2π) until the end.
- Obstruction block (D): GAP; ROUTINE. Item 1 (identity) rigorous; items 2–7 heuristic (contour remainders, stationary-phase reassembly, transfer of BHB's
  bounds and the ≥ 1/24 claim need proofs; the last needs ∫P ≥ 0); item 8 overclaimed. Wording fixed: the supported conclusion is "these routes do not
  prove H", with the failure located; θ = 0 open is a literature fact.
- Final assembly: GAP; ROUTINE. NO RESULT confirmed as the category; 18 wording items, all applied: "≥ T" (equality at θ = 0), "no counterexample found",
  sign of S_B − F, "coefficient ≥ 1/24" not "fraction", the κ < ½ − 2θ boundary only for θ < ¼, AFE errors stated, "formally related χ²-dual sum"
  instead of "self-referential", "formal substitution" for BHB's bounds, agreement "within .010" (not .007), "zero index 10¹²", θ = ½ as an extrapolation,
  the κ = 0 rate and the ratios-conjecture implication marked unverified, novelty only as "new numerical check in this report".

## Novelty
- Unmollified case: proved nowhere; conjectured (RH) by Hughes–Pearce-Crump, arXiv 2509.07792 (conj:secmom, conj:n1), with numerics on 10⁶ zeros; 2509.07788
  (same authors) proves only a twisted first moment. Durkan–Hughes–Pearce-Crump 2601.18025 treat Σ χ(ρ)X^ρ (χ, not χ²). Benli–Elma–Ng 2311.13554,
  Hughes–Pearce-Crump 2106.03005, Hughes–Martin–Pearce-Crump 2305.14253, Ng 0706.1763, Bui–Florea–Milinovich 2310.03949: first, absolute or negative moments.
  ([.scratch/hmean/scout_k2.md], [scout_lit.md], [scout_tools.md] at the time of writing.)
- Listings (python3 .scratch/lst.py 'simple zeros|discrete moment|mean value.*zeros|derivative of the Riemann zeta' math.NT --zenodo '"simple zeros"', 20 Nov):
  arXiv 2611.21551 (Najnudel–Wei, complex moments of ζ′ on vertical lines via random Euler products, not at zeros); Zenodo 21975237 (Yang–Yang, 79.62 %
  claim via Weil-matrix trace moments, "certified-candidate"), 24227401 (correction, auxiliary ratio only), others unrelated; none uses Σζ′(ρ)²M(ρ)² or Re Y.
- Repo: rg 'simple zero|19/27|Bui.Heath' TREE.md NOGO.md Riemann/Deep/*/REPORT.md: only explore55 and N8 record lines. The localization formula f(κ) is a routine
  recipe consequence (reviewer); not a discovery. Nothing here passes the DISCOVERIES gate.

## Caveat
- The obstruction is located, not proved: it says where the contour/Gonek/BHB route fails and that the missing part is a fixed fraction of the main term
  (recipe + numerics), not that H is unprovable from RH.
- The localization numerics test the recipe, not provability; the κ = 0 discrepancy (θ = ½: 1.009 vs 1.0625 at ℓ = 24.5) is presumably a lower-order
  term (not derived).
- H still has the support listed in explore55 (CUE, Odlyzko at 10¹², Davenport–Heilbronn); the θ = 0 constant now also matches Hughes–Pearce-Crump.
- Unverified conditional: RH + the CFZ ratios conjecture for Σ_γ ζ(ρ+α)ζ(ρ+β)/(ζ(ρ+z₁)ζ(ρ+z₂)) should give H for 0 < θ < ½ (B1–B4), hence
  κ* ≥ 7/8 − ε; uniform error control, the z-integral tails and the proportion calculation are not written: COROLLARIES material at most.

RH relevance: N8

## Suggested NOGO line
- hmean (Orca, proof attempt of explore55's H) NO RESULT: Σ_γ ζ′(ρ)²B(ρ)² is unproved even for B ≡ 1 (Hughes–Pearce-Crump's open k = 2 conjecture,
  2509.07792); up to AFE errors S_B = first moment Σζ′(ρ)A(ρ) with a twist of length T^{1+2θ} ≥ T (BHB: < T^{1/2}); recipe localization
  f(κ) = u² − 1/12 + μ(1−κ)² + (1−κ)³/3 (approximately supported at zero index 10¹², within 0.010) puts a coefficient ≥ 1/24 at moduli beyond BHB,
  a χ²-dual sum that is not small in the data; recipe constant intact
  [HMean](Riemann/Deep/HMean/REPORT.md)

## Files
[DERIVATION.md](DERIVATION.md), [NOTES.md](NOTES.md) (working notes), [kappa.py](kappa.py) (localization test; modes low/high), [recipe_check.py](recipe_check.py) (sympy),
[rs.py](rs.py) (Riemann–Siegel Z′, copied from explore55), [logs/](logs/). Data: Odlyzko's zeros1, zeros3 (downloaded to .scratch/hmean/data, deleted after the runs).
