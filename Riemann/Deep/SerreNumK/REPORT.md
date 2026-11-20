RESULT: Serre's open image step over Galois number fields K without GRH, in numbers. (A) Non-isogenous non-CM E₁, E₂/K have a good degree-one prime 𝔭 with a_𝔭(E₁) ≠ a_𝔭(E₂) and N𝔭 ≤ (aρ + b(K))^32, a = 41.256 + 41.472n_K, b(K) = (41.472 + 10.368n_K) log d_K + 41.472n_K log(384n_K) + 0.1187, ρ = log N_{K/ℚ} rad(2O_K𝔑_{E₁}𝔑_{E₂}). (B) For non-CM E/K and ℓ > 2(aρ_E + b(K))^16, G_E(ℓ) ⊂ C_ns^+(ℓ) forces G_E(ℓ) ⊂ C_ns(ℓ). (C) For K = ℚ(√D), 1 < D < 50, and image not in C_s^+(ℓ), ρ̄_{E,ℓ} is surjective for ℓ > max(10^70.36, 2(124.201ρ_E + b(K))^16), b(K) ≤ 876.93. (C) needed a second GRH input removed: Banwait–Derickx's Type 2 primes. Certified in Arb.

# SerreNumK: explicit Serre open image over Galois number fields without GRH (node N11)

Agent serrenumk, 2026-11-19, 23:33 UTC – 20 Nov 00:05 UTC (all compute single-threaded, ≤ 35 MB). Brief: `.scratch/orca/serrenumk-brief.md`.
Full proof: `DERIVATION.md`. Source: arXiv 2611.19556 (14 Nov 2026), whose every main statement assumes GRH.

## Claim
Input: chebexp (A). For every Galois L/K with L ≠ ℚ and every class C, some unramified degree-one 𝔭 has Frob_𝔭 = C and
N𝔭 ≤ (c·log 3Δ_L)^32, with c = 0.484, 0.29, 0.151, 0.108, 0.094 for log Δ_L ≥ 0, 10, 100, 1000, 10⁴. It rests on
Cureton–Thorner's Lean-checked ζ_L ≠ 0 for Re s > 31/32. Write ρ := log N_{K/ℚ} 𝔪.

**Theorem A (effective isogeny theorem over K; replaces the paper's Thm 'Eff iso for K').** Let K be Galois, and let
E₁, E₂/K be non-CM and not K-isogenous. Put 𝔪 = rad(2O_K𝔑_{E₁}𝔑_{E₂}). Then some prime 𝔭 ∤ 𝔪 with N𝔭 = p prime has
a_𝔭(E₁) ≠ a_𝔭(E₂) and

  N𝔭 ≤ (aρ + b(K))^32,  a = 0.108(382 + 384n_K),  b(K) = 0.108((384 + 96n_K) log d_K + 384n_K log(384n_K) + log 3).

Once (382 + 384n_K)ρ + (384 + 96n_K) log d_K + 384n_K log(384n_K) ≥ 11490, the factor 0.094 may replace 0.108.

**Theorem B (non-split Cartan step; replaces Thm 5.4).** Let E/K be non-CM and put ρ_E = log N rad(2O_K𝔑_E). For every
prime ℓ > 2(aρ_E + b(K))^16, G_E(ℓ) ⊂ C_ns^+(ℓ) implies G_E(ℓ) ⊂ C_ns(ℓ). If K is formally real, G_E(ℓ) ⊄ C_ns^+(ℓ).

**Lemma 2 (Momose Type 2 isogeny primes, new; replaces the GRH step of Banwait–Derickx).** Let K be real quadratic with
d_K ≤ 188. No K-rational p-isogeny of Momose Type 2 has log p ≥ 162.0, that is, p ≥ 10^70.356.

**Corollary C (replaces Cor. 'total theorem').** Let K = ℚ(√D) with D squarefree and 1 < D < 50, and let E/K be non-CM
with ρ̄_{E,ℓ}(G_K) ⊄ C_s^+(ℓ). Then ρ̄_{E,ℓ} is surjective for every ℓ > max(e^{162}, 2(124.201ρ_E + b(K))^16).

Constants from logs/cert_output.txt, rounded up. a₁ and b₁ are the 0.094 constants, valid once ρ ≥ ρ₁:

| K | d_K | a | b(K) | a₁ | b₁(K) | ρ₁ |
|---|---|---|---|---|---|---|
| ℚ | 1 | 82.729 | 246.904 | 72.005 | 214.898 | 12.017 |
| ℚ(√5) | 5 | 124.201 | 651.302 | 108.101 | 566.874 | 4.749 |
| ℚ(√2) | 8 | 124.201 | 680.540 | 108.101 | 592.322 | 4.513 |
| ℚ(√13) | 13 | 124.201 | 710.742 | 108.101 | 618.609 | 4.270 |
| ℚ(√47) | 188 | 124.201 | 876.930 | 108.101 | 763.254 | 2.932 |

The other 25 fields ℚ(√D), D < 50, are in the log: b(K) runs from 651.302 to 876.930.

| K | N rad | Thm A: log₁₀ N𝔭 ≤ | GRH (paper) | Thm B: log₁₀ ℓ ≤ | GRH | Cor. C: log₁₀ ℓ ≤ |
|---|---|---|---|---|---|---|
| ℚ | 10 / 10⁶ / 10¹⁰⁰ | 84.51 / 98.65 / 135.21 | 7.66 / 8.54 / 10.76 | 42.56 / 49.62 / 67.90 | 4.13 / 4.57 / 5.68 | (ℚ: serrenum) |
| ℚ(√2) | 10 / 10⁶ / 10¹⁰⁰ | 95.53 / 106.22 / 141.00 | 8.33 / 9.02 / 11.12 | 48.06 / 53.41 / 70.80 | 4.47 / 4.81 / 5.86 | 70.36 / 70.36 / 70.80 |
| ℚ(√5) | 10 / 10⁶ / 10¹⁰⁰ | 95.10 / 106.05 / 140.99 | 8.31 / 9.01 / 11.12 | 47.85 / 53.32 / 70.79 | 4.46 / 4.81 / 5.86 | 70.36 / 70.36 / 70.79 |
| ℚ(√47) | 10 / 10⁶ / 10¹⁰⁰ | 97.07 / 107.31 / 141.09 | 8.47 / 9.08 / 11.13 | 48.84 / 53.96 / 70.85 | 4.53 / 4.84 / 5.87 | 70.36 / 70.36 / 70.85 |

The table uses the tight envelope G(R)^32 (DERIVATION §3). The GRH columns use the paper's C₂ as its proof gives it.
In Corollary C the Type 2 threshold 10^70.36 dominates until N rad(2O_K𝔑_E) > about 10^93.

## Before
- **GRH, the source.** arXiv 2611.19556 (14 Nov 2026) assumes GRH throughout. It proves
  N𝔭 ≤ (C₁ρ + C₂)², C₁ = 510.352 + 513.024n_K, and ℓ ≤ 2C₁ρ + 2C₂.
  - The GRH enters only through Bach–Sorenson Thm 5.1 and Table 1 (l.368–375, 1026–1056), and through the
    Banwait–Derickx Type 2 step behind the paper's Thm 'bound theorem'. That second use is not stated in the paper (§5 below).
- **Unconditional, over ℚ.** serrenum (COROLLARIES): p ≤ (17 log R + 100)^32 and C(E) = 2(11.5 log rad(2N_E) + 68)^16.
- **Unconditional, over number fields.**
  - Masser–Wüstholz, made explicit by Lombardo (2015), give bounds polynomial in the Faltings height and [K:ℚ], not
    in the conductor alone.
  - Unconditional Chebotarev bounds are powers of Δ_L (Kadiri–Wong Δ_L^310), so they give powers of N rad 𝔑.
  - Larson–Vaintrob (GRH): c(E,K) ≪_K log N_E.
  - No bound polylogarithmic in N rad(𝔑_E) over K ≠ ℚ without GRH was found (Novelty).
- **Isogeny primes over ℚ(√D), 1 < D < 50.** Banwait–Derickx Thm 1.9 (ANT 2025, arXiv 2203.06009v3) lists them, all
  ≤ 163. Its correctness needs GRH for Momose Type 2 primes (TeX l.1457: "the only result in our work which requires
  the Generalised Riemann Hypothesis"). The statement itself does not say so.

## Method
- **Kept algebra.** The 2-adic deviation group with [L:K] ≤ 192 (Mayle–Wang), the paper's L̃ = L·K(√α) with
  v_𝔮(α) = 1 for every 𝔮 | 𝔪, the paper's discriminant bound log d_{L̃} ≤ R(K, ρ), and Hensel–Serre. All re-derived in
  DERIVATION §2; R's coefficients are rechecked exactly in cert.py.
- **Replaced step.** Bach–Sorenson on L̃/K is replaced by chebexp on the same L̃/K and the same lifted class.
  - chebexp's c(Λ) jumps down at 10, 100, 1000 and 10⁴, so it is used through the envelope
    G(R) = sup_{Λ≤R} c(Λ)(Λ + log 3) ≤ 0.108(R + log 3) for R ≥ 1400. Here R ≥ 2285 always.
  - Bach–Sorenson's 2.5n_L term has no counterpart.
- **Theorem B.** The paper's twist argument, with two repairs.
  - ε is unramified at 𝔭 | ℓ over K (Raynaud/Serre inertia with e ≤ n_K).
  - 𝔭₀ | ℓ is impossible, because ordinary reduction at a degree-one 𝔭 | ℓ puts (u *; 0 1) with u ≠ ±1 in the image.
    The paper skips this step, as Mayle–Wang do over ℚ.
- **Lemma 2 (new).** It uses Banwait–Derickx's unconditional Condition CC (their Prop. 'cond:CC') with two chebexp primes.
  - q₁ splits completely in F = ℚ(√D, √−p), with log Δ_F ≤ 2 log p + 2 log d_K + log 16.
  - q₂ ≠ q₁ splits completely in F(√q₁).
  - Both are < p/4 once log p ≥ 162. They cannot both be primitive cube roots of 1 mod p, since p ∤ q₁ + q₂ + 1.
    So CC fails for p.
  - Arb chain: 7352 intervals of width ¼ on [162, 2000], then a monotone tail.

## Certificate (one command)
```
cd Riemann/Deep/SerreNumK && uv run --no-project --with python-flint python3 cert.py > logs/cert_output.txt
```
- Ends with `ALL CERTIFIED (0.6 s)`, one core, max RSS 35 MB (logs/cert_time.txt).
- Every decimal enters as an Arb ball built from a string, and every comparison is a certain Arb comparison.
- It checks:
  - the exact coefficients of R;
  - the envelope thresholds 1400 and 11490, and ĝ(R) = 0.108R + 152.3 ≥ G(R);
  - a, b, a₁, b₁ and ρ₁ for ℚ and the 30 fields ℚ(√D), against the identity 0.108(R + log 3) = aρ + b;
  - the Type 2 chain and its tail.
- It also prints the table above.

## Checks
- **GRH audit.** rg 'GRH|Riemann|Bach|Lagarias|Odlyzko|conditional|assum' over the TeX; every hit was read
  (DERIVATION §1).
  - GRH enters only through Bach–Sorenson. That covers Lemma 'bach norm prime' and Thm 'Eff iso for K', and through
    them Thm 5.4 and the corollary.
  - It also enters through the cited Banwait–Derickx Thm 1.9 (Borel step, l.290). That use is closed by Lemma 2.
  - The other cited results are algebraic or computational: Mayle–Wang Prop. 12 and Lemma 15, Serre–Tate, Dickson,
    and Bilu–Parent–Rebolledo (not used for K ≠ ℚ). Banwait–Derickx's Type 1, generic and modular-curve steps are
    unconditional by their own statement.
- **Sanity (PARI, logs/sanity_output.txt, 9 s, 14 MB).** Over ℚ(√5) and ℚ(√2), 24 curves each: 8 base changes plus
  16 twists by non-rational d. That gives 276 non-isogenous pairs per field.
  - The least separating good degree-one prime is ≤ 31 (ℚ(√5)) and ≤ 41 (ℚ(√2)).
  - log p₀ − 32 log(aρ + b) ≤ −224.
  - Twist step: the least good degree-one 𝔭 with χ_D(𝔭) = −1 and a_𝔭 ≠ 0 is ≤ 47.
  - The 5-, 17- and 37-isogeny curves (11a1, j = −17·373³/2¹⁷, j = −7·11³) stay non-surjective over K at
    ℓ = 5, 17, 37. Theorem B's bound there is 10^50.4–10^52.7.
  - Hensel–Serre (4.4) holds on the 8 fields K(E[2]), with equality where L = K.
- **LMFDB** (API, JSON, logs/lmfdb_check.txt). 400 curves over ℚ(√5) with conductor norm ≤ 256: 378 are non-CM
  with a non-maximal prime, the largest being 7 (2.2.5.1-41.2-a2, not a base change). That is far below every bound.
  The ℚ(√2) pages returned no rows.
- **K = ℚ.** Theorem A gives (82.729ρ + 246.904)^32, serrenum (17ρ + 100)^32: the same shape, as required.
  The constant is about 5× worse because the paper's L̃ doubles the degree and pays log|Nα| ≤ 2n_Kρ + ….
- **Paper slips found (minor, for REFUTATIONS).**
  1. The stated C₂ (l.153–156, 939–942) has 513.024 where its proof (l.1094–1100) gives 513.024 log d_K. The stated
     constant is too small by 513.024(log d_K − 1), e.g. 554 for ℚ(√47).
  2. Thm 5.4 claims "not in C_ns^+", but its proof covers only "in C_ns^+, not in C_ns". The claim is fine for
     formally real K.
  3. 𝔭₀ | ℓ is not excluded, and neither is ε being ramified at 𝔭 | ℓ over K (Lemma 15 is only cited). Both are closed
     in DERIVATION §4.
  4. Thm 'bound theorem' and Cor. 'total theorem' silently depend on GRH through Banwait–Derickx Thm 1.9.

## Review
One `reviewer` subagent saw R1–R6 pasted with exact formulas (logs/review.txt). Verdicts:
- R1 envelope, R2 kept algebra and R, R4 Theorem B with the 𝔭₀ ∤ ℓ step, R5 Lemma 2 and R6 Corollary C: CORRECT.
- R3 Theorem A: GAP, only because the paste rounded b(ℚ) = 246.9037… to 246.9. cert.py rounds up and certifies 246.904,
  which is the value used here.
- Overall: ROUTINE. "The argument transfers the supplied inputs rather than establishes a new arithmetic input."

## Novelty
Searched 2026-11-19 (logs/novelty_lst.txt, logs/novelty_zenodo_serre.txt).
- **arXiv listings.** math.NT, 2026-10 and 2026-11: `lst.py 'Serre|open image|surjectiv|isogen|Galois representation'`
  gives 17 title hits.
  - Only 2611.19556 bounds Serre's constant over number fields, and it assumes GRH.
  - The rest are Drinfeld modules, lifting, Serre weights, modularity and 7-adic images over ℚ (2610.06521).
- **Zenodo.** "open image": 34 records, none on Galois representations. "Serre constant": 0.
- **Our files.** rg 'Serre|open image|Larson|Momose|Banwait' over TREE, COROLLARIES, NOGO, DISCOVERIES, RECORDS,
  REFUTATIONS and Riemann/Deep. Hits:
  - SerreNum (over ℚ only);
  - R492 (the candidate that set this task);
  - Harvest4/HARVEST.md item 5 ("Larson–Vaintrob without GRH", planned, not done; it proposed LV's own field
    K(√−ℓ, √−3), giving ℓ ≲ (log)^64; Lemma 2 goes through Condition CC instead, giving (log)^32 and 10^70.4);
  - ChebExpAudit (the input).
- **Assessment.**
  - Theorems A, B and Corollary C are a routine corollary of chebexp and 2611.19556, as serrenum was of Mayle–Wang.
  - New small items:
    - the second, unstated GRH use (Banwait–Derickx Type 2), and its unconditional closure for real quadratic
      fields (Lemma 2: two chebexp primes plus the cube-root trick);
    - the 𝔭₀ ∤ ℓ step and ε unramified at 𝔭 | ℓ over K;
    - the C₂ typo.
  - Candidate for COROLLARIES (the orchestrator decides).

## Caveat
- **Inputs.**
  - chebexp: Cureton–Thorner's 31/32, Lean-checked, paper unrefereed; chebexp itself audited, PASSED.
  - The paper's algebra, re-derived in DERIVATION §2. The integral-basis claim |σw_j| ≤ n_K√d_K is checked only
    through Minkowski II plus a reduced basis.
  - Mayle–Wang's [L:K] ≤ 192 is taken over K unchanged (pure group theory).
  - The paper's exceptional-image and C_ns cases for K = ℚ(√D) (l.292–330) were read, not re-proved.
  - Banwait–Derickx's unconditional parts are cited, not re-run: Condition CC, the generic and Type 1 bounds, and the
    Magma determinations behind Thm 1.9.
- **Size.**
  - The exponents 32 and 16 make the numbers huge: Theorem B is 10^48–10^71 against 10^4.1–10^5.9 under GRH.
  - In Corollary C the Type 2 threshold 10^70.36 dominates until N rad > 10^93.
- **Not optimal.**
  - The paper's L̃ = L(√α) doubles the degree, and its bound log|Nα| ≤ 2n_Kρ + … is crude.
  - serrenum's excluded-prime term inside chebexp's master inequality would avoid L̃ and cut a from 124 to about 20
    for quadratic K. That needs new Arb chains and was not done.
  - The Type 2 threshold drops to about e^117 for p ≡ 2 mod 3 (one prime suffices). This is a float estimate, not used.
- **Scope.** Theorem A needs K Galois, as in the paper. Corollary C needs K real quadratic with D < 50.
- **Routine?** Yes relative to chebexp. The substitution is mechanical. The new content is that the corollary had
  a hidden GRH input, now removed.

RH relevance: N11

## Files
- `REPORT.md` (this file) and `DERIVATION.md` (every replaced step and constant: §1 GRH audit, §2 kept algebra,
  §3 Theorem A, §4 Theorem B, §5 Lemma 2 and Corollary C, §6 remarks).
- `cert.py`, the Arb certificate.
- `sanity.gp`, PARI sanity checks.
- `logs/`: cert_output.txt, cert_time.txt, sanity_output.txt, sanity_time.txt, novelty_lst.txt, novelty_zenodo_serre.txt,
  lmfdb_check.txt, review.txt.
