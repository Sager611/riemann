NEW: no source states Theorem C, Proposition F or Theorems D, D⁺, D_M, or any bound with an explicit power of |D| (per D or for almost all D) for primes in ideal classes, or represented by a form of discriminant D, with a sign or Sato–Tate condition on a_E(p) or λ_f(p). Caveat: a per-D bound |D|^C with an inexplicit C of order 10^10 looks derivable from published inputs (see below) but is stated nowhere.

# RSFam novelty audit 2 (2026-10-21, after review)
Scripts and raw outputs: `.scratch/rsfamfix/` (`axq.py` → `axq.out`, `axq2.out`, `other.py` → `other.out`, `oa_auth.out`, `oa_more.out`, `zenrec.py`, `gxl.out`, TeX in `src/`). First audit: [NOTE](NOTE.md) §8.

## What was searched
- arXiv API: 65 new queries (Sato–Tate × arithmetic progressions / Linnik / imaginary quadratic / ideal class / Hecke characters / binary quadratic / dihedral / Chebotarev / short intervals / almost all / base change; sign changes × quadratic forms; least prime × elliptic curve / Hecke eigenvalue; log-free × Rankin–Selberg / Hecke characters; author lists of Zaman and Zywina; 2 failed) and the author lists of Lemke Oliver and Wong (`axq2.out`). The remaining 18 queries (Thorner's author list and phrase variants) hit the API rate limit and then timed out; they were replaced by 15 OpenAlex searches (`oa_more.out`) and Thorner's OpenAlex work list.
- OpenAlex (journals outside arXiv): 10 + 15 searches (`other.out`, `oa_more.out`), and the full work lists of Thorner, Wong, Zaman and Lemke Oliver (`oa_auth.out`).
- Zenodo: 15 queries (2025–2026) and the 6 closest records read (Lee 23452649, 22168461, 21743631, 22907575, 23118170, 22112048).
- Hexagon (hexagonmath.org API): all 323 works (none new since the RCApps2 audit), keyword scan of titles and abstracts.
- openai/math preprints: 749 (none new), keyword scan.
- GXL evidence API: 5 claims (least prime represented by a form with a_E(p) > 0 polynomial in |D|; effective joint Chebotarev–Sato–Tate polynomial in the discriminant; Sato–Tate Linnik in progressions; almost-all imaginary quadratic fields with positive-λ_f prime ideals in every class; zero density for f ⊗ θ_χ over class group characters).
- TeX read: Thorner 2404.06482 (IMRN 2025), Thorner–Zhao 2601.04189, Akbary–Wong 2103.15713; first audit: BTZ 1804.06402, Lemke Oliver–Thorner 1505.03122, Humphries–Thorner 2103.05634, Pascadi–Thorner 2508.14888.

## Closest results
| Work | Proves | Why it is not Theorem C/F/D |
|---|---|---|
| Lemke Oliver–Thorner, IMRN 2019 (1505.03122), least-prime Sato–Tate theorem and Remark 5 | least prime 𝔭 of a totally real K with cos θ_𝔭 ∈ I, polynomial in 𝔮(π) and D_K; remark: the same for p ≡ a mod q | no ideal-class condition; the family log-free density they need exists only over Q with Dirichlet characters |
| Brumley–Thorner–Zaman, JEMS 2022 (1804.06402) | log-free density for families π × π₀ (Thm 1.2), Landau–Siegel repulsion (Thm 1.4) | input (J4) of the NOTE; no Sato–Tate or class application |
| Thorner, IMRN 2025 (2404.06482) | no Landau–Siegel zero for Sym^m π × (Sym^n π′ ⊗ χ), F totally real; effective (joint) Sato–Tate rates | no class or progression condition; used as input (J5) for m ≥ 3 |
| Thorner–Zhao 2026 (2601.04189) | no exceptional zeros for new Rankin–Selberg and triple-product families; self-twist theorem | zero-free regions only; their Prop. 1.1 (4), (6) (Ramakrishnan–Wang, Luo) and Thm 1.7 exclude exceptional zeros of Sym^m f ⊗ θ_χ (NOTE §5.6, not used) |
| Akbary–Wong 2021 (2103.15713) | Siegel–Walfisz (N𝔮 ≤ (log x)^γ) and short-interval (N𝔮 ≤ x^θ) sums of Rankin–Selberg coefficients over primes in ray classes | fixed base field (constants depend on F, suppressed); coefficient sums, no sign condition; class groups of varying K not covered |
| Pascadi–Thorner 2025 (2508.14888) | general family zero density | weaker for this family (NOTE §3, §4); no application to classes |
| Thorner 2021 (2002.10450); Chen–Park–Swaminathan (1507.07122, 1506.09170); Hammonds et al. 2021 (1906.07903, GRH) | least prime with θ_p ∈ I: level aspect, CM curves, progressions | no ideal classes |
| Wong, JNT 2019; Amri (2203.11498); Malik–Prabhu (2012.12534); Bucur–Kedlaya (1301.0139); Lee (2607.18832); Yang (2609.07080) | Chebotarev–Sato–Tate for a fixed Galois extension | fixed field; for ring class fields every bound is exponential in h(D) |
| Thorner–Zaman, ANT 2017 (1604.01750); Zaman's thesis | least prime represented by a form, \|D\|^694, \|D\|^455; elliptic-curve applications via division fields | no Sato–Tate or sign condition |
| Vaishya et al. (2204.08310, 2401.18055, 2601.16352); Hua (2023–2026) | λ_f at integers represented by one fixed form | fixed form, not primes, no uniformity in D |
| Kumar–Kumari–Mishra (2604.17532); Hamdar–Wang (2604.24753) | effective joint Sato–Tate of two forms; vertical laws | no classes |
| Zenodo, Lee 23452649 / 22168461 / 21743631 (unrefereed) | local Sato–Tate in progressions with moduli up to exp((log X)^β); joint Chebotarev–Sato–Tate for a fixed field | progressions only, moduli far below any power of X |

## Caveat: a per-D bound is assembled from known pieces, not stated
Per D, Linnik's method seems to give "every class has a prime p ≤ |D|^C with λ_f(p) > 0" from: BTZ Theorem 1.2 (log-free density, family {π_χ : χ ∈ Ĝ_c} of dihedral GL₂ forms of conductor ≤ 12|D|, π₀ = Sym^m f, exponent 10^7(2(m+1))^4(1−σ)), BTZ Theorem 1.4 (repulsion by a Siegel zero of a genus factor), exceptional zeros excluded by Ramakrishnan–Wang (m = 1), Luo (m = 2) and Thorner–Zhao for the non-genus twists and by Hoffstein–Ramakrishnan, Banks, Thorner 2025 for the genus twists, and Thorner–Zaman for m = 0. C would be at least of order 10^10. We did not write this out and found it stated nowhere. Theorem D's content is the explicit exponent 3 (9/4 for one-sided conditions) for almost all D, in short intervals; Theorem C and Proposition F are new statements by a classical method.

## Residual risk
arXiv full text was read for about ten sources only; 18 arXiv API queries could not run (rate limit, then timeouts) and were replaced by OpenAlex searches; journals outside arXiv were seen only through OpenAlex metadata; Zenodo indexes metadata only.
