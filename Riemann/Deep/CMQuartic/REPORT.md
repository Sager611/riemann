DISCOVERY: class group characters of CM fields over a fixed totally real field F have a thin-family large sieve with constant ≈ N + √d_K (new relative lattice lemma). Consequences: for each such field Σ_{χ∈Ĉl_K} N(σ;t,χ) ≪ (d_K Tⁿ)^{2(1−σ)/(3−2σ)+ε}; and for every A > 3/2, almost all CM fields K/F with prime (or δ₀-rough) relative discriminant have prime ideals of degree one in every ideal class in every interval (x, x + x·d_K^{−κ}], x ≥ d_K^A. The previous bound, valid for every field, was d_K^{694} (Thorner–Zaman 2017).

# CMQuartic report (2026-10-21)
Proofs: [NOTE.md](NOTE.md). Paper proofs, not refereed, no Lean. F totally real of degree n, fixed (quartic CM: n = 2); K/F CM; D = N𝔡_{K/F}; d_K = d_F²D.

## Answer to the brief
- **The large sieve exists (Lemma 3).** Weight each ideal by its O_F-content (the largest ideal of F dividing it) and use the characters trivial on
  ideals from F; ≤ h_F cosets of them cover Ĉl_K. Constant C_F ε⁻¹log^{2n−1}(3d_K)(N + √d_K); both terms are attained.
  - Lemma 1: for β, β′ in an ideal 𝔠 with β′/β ∉ F, ββ′ᶜ − βᶜβ′ ∈ N_{K/F}(𝔠)𝔇_{K/F}, so its norm is ≥ (N𝔠)²d_K/d_F².
  - Lemma 2 (angle boxes + Lemma 1): a class holds ≤ 2^{n−1}(1 + (4π)ⁿe^{nr_F}d_F N/√d_K) ideals of norm ≤ N not divisible by an ideal from F.
- **Why Density2 stopped.** With Z-primitivity the method fails for CM fields as in degree 3: the ideals 𝔭O_K (𝔭 of degree one in F) crowd into
  ≤ h_F classes, a loss of d_K^{1/4−o(1)} (Remark 2.2). The O_F-content weight removes them; Density2's L² pair bound then holds (Remark 2.4).
- **Zero density.** Theorem A: Σ_{χ∈Ĉl_K} N(σ;t,χ) ≤ C(F,ε)(d_K Tⁿ)^{2(1−σ)/(3−2σ)+ε}, 0.6 ≤ σ ≤ 1 (Pascadi–Thorner Thm 1.3: d_K^{3(1−σ)/(3−2σ)}).
  Proposition C: summed over the fields with D ≈ X, X^{10(1−σ)/(3−2σ)+ε} (Humphries–Thorner over F).
- **Application (Theorem B, Corollary B1).** Family: CM K/F with odd 𝔡_{K/F} having no prime factor of norm ≤ D^{δ₀} (includes prime 𝔡_{K/F},
  e.g. K = F(√−π)). For every A > 3/2, all but O(X^{1−c}) such K with D ≤ X have ≥ H/(16h_K log 2x) primes of degree one in every class with norm
  in (x, x+H], for all x ≥ D^A, H ≥ xD^{−κ}. So some β ∈ O_K has prime norm ≤ 2d_K^A; every β ∈ O_K∖F has norm ≥ d_K/(4ⁿd_F²), so exponent 1 is
  necessary. For n = 2, A ≥ 1.55: κ = min(10⁻³, δ₀/4, κ₀(F)), c = min(4·10⁻³, δ₀/3). No Siegel-zero loss: genus factors are fixed (handled by the
  Thorner–Zaman log-free density for one fixed Dedekind zeta function) or have conductor > D^{δ₀}.
- **Conditional (routine; COROLLARIES.md).** Under Hu's half-plane 7/8 (Liu's B♭): every class of every CM K/F has a prime ≪ d_K^{8/5+ε} (d_K^{1.5997}).

## Where it stops (NOTE §7)
- All CM K/F: small genus conductors need a Q-aspect log-free density with Deuring–Heilbronn for Hecke characters of F (known over Q only).
- F varying: Lemma 2 and Theorem A are polynomial in d_F (sketch), but the Humphries–Thorner and Thorner–Zaman constants are per F.
- A < 3/2: same three-way coincidence at σ = 5/6 as RCApps2 §4. Fields with no quadratic subfield: no involution; Density2's counterexample stands.

## Checks (pari/gp, .scratch/cmquartic/, all pass)
- Lemma 1: Δ ∈ 𝔇_{K/F} for 1800 random pairs in 6 fields. Lemma 2/3: class counts ≤ 1/175 of the bound; h_K·max/(N + √d_K) ≤ 0.61 (10 fields,
  d_K ≤ 1.5·10⁸, h_K ≤ 408). Least prime per class: worst is the principal class, 1.04–1.97 × d_K/(16d_F²). Parameters: params_cm.py.

## Novelty
- About 45 arXiv API queries and 3 GXL claims: nothing implies Lemmas 2–3, Theorem A, Proposition C or Theorem B. Closest: Habegger–Pazuki 1506.02485
  (average class counts via subconvexity), Sardari 1802.06193 (GRH), Thorner–Zaman 1604.01750 (d_K^{694}), almost-all Chebotarev for fixed groups.
- New idea: Lemma 1 (the CM involution gives an O_F-determinant controlled by the relative different). The rest reuses Combos, Density2, RCApps(2).

## Leads
- Q-aspect log-free density + Deuring–Heilbronn over F (all CM K/F); F-uniform Humphries–Thorner (all quartic CM fields); infinite-type Hecke
  characters of K (CM Hilbert forms of higher weight) via angular boxes in each embedding.
