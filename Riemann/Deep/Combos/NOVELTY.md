NEW: no published result implies Theorem 1 (Σ_χ N(σ;t) ≪ |D|^{2σ(1−σ)+ε}T³ over the ring class characters of D). Every published density bound for this family is weaker in |D| for every fixed σ in (1/2, 1). The method is classical (Carlson's argument plus a weighted large sieve), so what is new is the statement for this family, not the technique.

# Novelty of the ring class zero-density theorem (2026-10-20)

Searched: about 105 arXiv API queries (class group, ring class, Hecke, dihedral or quadratic form, combined with zero density, large sieve, least prime or Linnik; hits in `.scratch/combosv/ax/all.json`); the arXiv math.NT listings for Aug–Oct 2026 (1938 titles); the Hexagon site (16 terms); the GXL evidence API (4 claims). Full texts searched with rg: 21 arXiv sources (`.scratch/combosv/src/`), Kowalski–Michel 2002, Blomer 2004, and the theses of Zaman (Toronto 2017) and Gaudet (Rutgers 2023). Closest published bounds for Σ_χ N(σ,T,χ), with χ over the ring class characters of D = f²d_K:

| Source | Bound | Compared with Theorem 1 |
|---|---|---|
| Pascadi–Thorner, arXiv 2508.14888 (2025), Thm 1.3 (any finite set of GL₁ forms over K) | (\|D\|⁵\|d_K\|^{−2}T⁸)^{(1−σ)/(3−2σ)+o(1)}; for fundamental D, \|D\|^{3(1−σ)/(3−2σ)} | weaker for every σ, since 3/(3−2σ) > 2σ (because 4σ²−6σ+3 > 0); better in T |
| Kowalski–Michel, Pacific J. Math. 207 (2002), Thm 2 (their Cor 7 treats x²+ny², exponent unspecified) | T^B\|D\|^{c(1−σ)/(2σ−1)}, c > 5.5, σ ≥ 3/4 | weaker for fixed σ |
| Weiss 1983, Cor 4.4; Thorner–Zaman, arXiv 1510.08086, Thm 1.1 | (e^{O(1)}d_K²QT²)^{81(1−σ)}, with d_K²Q ≥ \|D\| | weaker for fixed σ < 1 |
| An, arXiv 2004.14410; Humphries–Thorner, arXiv 2103.05634; Thorner–Zaman, arXiv 1906.07717 | exponents ≥ 6.5(1−σ) | weaker |

Kowalski–Michel, Weiss/Thorner–Zaman and An are log-free, so within O(1/log|D|) of σ = 1 they beat Theorem 1, which loses |D|^ε. Theorem 2 does not use that range.

- Large sieve. Specialising Huxley 1968, Schaal 1970 and Schumer 1986 to this family gives constant ≍ N + |D| (as quoted by Ditchen, arXiv 1312.1502, and Gaudet, arXiv 2408.01610 §6). Pascadi–Thorner Thm 1.1 also gives ≍ N + |D|. Lemma 2 (N + |D|^{1/2}, with weight m(𝔞)^{1+2ε}) was not found anywhere. Its key step, writing 𝔞 = m·(primitive ideal) with χ(m𝔞) = χ(𝔞), is classical: Blomer (Ann. Inst. Fourier 54 (2004), §3) uses it for moments.
- Related work with no zero-density estimate: Duke–Friedlander–Iwaniec 1995 (first moment, subconvexity); Blomer 2004 (non-vanishing); Fouvry–Iwaniec 2003 (low-lying zeros); Ditchen 2013 (large sieve averaged over D); Gaudet 2024 and his thesis (sieve, explicit L = 1.3·10^10); Hexagon 2610.00178 §54 (Thorner–Zaman 2019, A unspecified); FSW 2610.18491; Liu 2610.26835.
- Records for Theorem 2, which is conditional on (L). The unconditional record is |D|^455 (Zaman's thesis, Cor 1.4.1, which Gaudet calls the record), not the |D|^694 the report quoted. (L) combined with Pascadi–Thorner would give |D|^{3/(3−2B♭)+ε} = |D|^{2.40} for fundamental D (up to |D|^{4.0} in general). (L) combined with Kowalski–Michel gives about |D|^{11}. Neither bound is stated anywhere. Theorem 2's exponent 2B♭+ε (1.7498) remains new.
- Gaps closed in NOTE.md:
  - Lemma 2 is now called weighted, with an example showing why the weight is needed.
  - Theorem 1's height dependence is stated: T³ per window, T⁴ up to height T.
  - The exponent reads 2B♭+ε throughout.
  - The factor y^{−0.4} is kept in the low-σ bounds of Theorems 2 and 8.
  - The small-|D| case of Theorem 2's interval claim is proved via Theorem 1(b).
  - (L) is marked on Theorem 2, Corollary 4 and the records.
  - The record is corrected from 694 to 455.
  - The ε bookkeeping of Lemma 4 is fixed.
  - REPORT.md was updated to match.
- Significance: Theorem 1 is new but routine. Theorem 2 is a new conditional record: 455 unconditionally, 1.75 under Liu's unrefereed half-plane.
