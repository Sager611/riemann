# SerreNum: novelty check (2026-11-06)

Claim checked: an unconditional bound for Serre's open-image constant C_E that is polylogarithmic in the conductor
(C_E ≤ 2(11.5 log rad(2N_E) + 68)^16 for every non-CM E/ℚ), and an unconditional effective isogeny theorem with a
polylogarithmic prime (p ≤ (17 log rad(2N_E N_E') + 100)^32), both with numbers.

## Searches
- **arXiv monthly listings** (`.scratch/lst.py 'Serre|open image|surjectiv|isogen|Galois representation' math.NT`),
  2026-10 and 2026-11: 13 title hits; none bounds C_E or the effective isogeny prime (Drinfeld modules, p-parity,
  7-adic images via Kummer descent 2610.06521, Serre weights, lifting, modularity). Zenodo `"open image"`: 31 records,
  none relevant.
- **arXiv API by author** (all papers, newest first): Zywina (latest 2602.10865; open-image work 2206.14959,
  2403.16147: algorithms, no conductor bound), Lemos (2002.02714, non-split Cartan), Le Fourn (2505.00680 etc.),
  Furio–Lombardo (2305.17780; 2603.04021 with Bisatt: p-adic images given a non-split Cartan mod p; no C_E bound),
  Mayle (2109.08656 GRH; 2212.11472 GRH, abelian varieties; 2501.02345 smallest surjective prime ≤ 7, a different
  question), Chen–Swidinsky (2501.00202, GRH).
- **arXiv API by phrase**: "open image" + "Riemann hypothesis" → 2507.06913, 2501.00202, 2109.08656; "isogeny theorem"
  + effective → 2501.00202; ti:"open image" → nothing new; Cureton + Thorner (all papers citing the Lean 31/32
  theorem) → 2611.02280, 2610.37179, 2610.18491 (primality, sparse polynomials, class groups): none on elliptic curves.
- **2507.06913 (Wang–Wei, Jul 2025)**: unconditional, but for 100% of pairs (zero-density, families ordered by
  conductor) and "comparable to GRH" only for semistable families; not every curve. Not a pre-emption.
- **Zenodo** phrase queries: "Serre constant" 0 hits; "effective isogeny" 0; "open image theorem" 2 (unrelated);
  "surjectivity" "quasi-Riemann" and "Cureton" "Thorner" return unrelated OR-matches.
- **HAL** (2025–2026): "open image" AND (Serre OR elliptic): 0; "isogeny" AND "Chebotarev": 0.
- **Our files**: `rg -n -i 'open image|Serre constant|Mayle|Swidinsky|2501.00202|2109.08656'` over TREE.md,
  COROLLARIES.md, NOGO.md, DISCOVERIES.md, REFUTATIONS.md, Riemann/Deep: only Harvest4/HARVEST.md (candidate 3, the
  brief's source). ChebExp and ArtinHalfPlane give the Chebotarev input only.

## Prior unconditional bounds (what this improves on, in shape)
- Masser–Wüstholz 1993 (C_E ≤ c·h^γ, height), explicit: Pellarin 2001, Lombardo 2015, Le Fourn 2016.
- Kraus 1995: C_E ≤ 68 rad(N_E)(1 + log log rad N_E)^{1/2}; Cojocaru 2005: (4√6/3) N_E ∏(1 + 1/p)^{1/2} + 1.
- Zywina 1508.07661 (Prop. 1.8, Thm 1.10): C_E ≤ max(17, (p+1)/2) for each p with v_p(j_E) < 0; ≤ max(17, log den j_E);
  ≤ max(37, (2√3/3)N_E^{1/2}∏(½ + 1/(2p))^{1/2}) if E has no multiplicative prime.
- Effective isogeny: linear in the level (Sturm bound for the two newforms; Kraus, Cojocaru, Zywina 1508.07661
  §4, Riemann–Roch lemma), or Chebotarev with power-of-discriminant bounds (Kadiri–Wong Δ^{310}: a power of R).
- GRH: Serre 1981 Thm 21–22; Mayle–Wang (482 log R + 2880)², C_E ≤ 964 log rad(2N_E) + 5760; Chen–Swidinsky
  (124 log R + 561)², C_E ≤ 446 log rad(2N_E) + 2254 (two of their steps have gaps, DERIVATION §6).

## Assessment
No unconditional bound polylogarithmic in the conductor was found, for C_E or for the effective isogeny prime.
The new input is Cureton–Thorner's Lean-checked half-plane (Oct 2026) through chebexp; the substitution into the
Mayle–Wang / Chen–Swidinsky argument is direct, so a referee may call it routine. New beyond the substitution: the
excluded-primes Chebotarev lemma (no √m doubling, degree 192 instead of 384), the 2-group bound 128 for reducible
residual pairs, the case split (a)–(d) with degrees 36/48/128/192, and the closing of four gaps in the sources.
Numerically the bounds beat the unconditional linear ones only for astronomically large conductors
(rad N_E > 10^48 against Kraus; single bad prime p > 10^51 against Zywina), see REPORT.md.
