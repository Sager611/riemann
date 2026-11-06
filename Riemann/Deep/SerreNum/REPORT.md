RECORD: first unconditional bounds polylogarithmic in the conductor for Serre's open image theorem and the effective isogeny theorem over ℚ, with numbers. For every non-CM E/ℚ, ρ̄_{E,ℓ} is surjective for every prime ℓ > 2(11.5·log rad(2N_E) + 68)^16 (> 2(4.5·log rad(2N_E) + 23.5)^16 if ρ̄_{E,2} is surjective). For non-isogenous E, E'/ℚ, some prime p ∤ 2N_EN_E' with p ≤ (17·log rad(2N_EN_E') + 100)^32 has a_p(E) ≠ a_p(E'). The only analytic input is Cureton–Thorner's Lean-checked half-plane, through chebexp. Certified in Arb.

# SerreNum: Serre's open image theorem with polylogarithmic numeric bounds and no GRH (node N11)

Agent serrenum, 2026-11-06, 10:54–12:15 UTC (all compute single-threaded, < 120 MB). Brief: `.scratch/orca/serrenum-brief.md`. Full proof: `DERIVATION.md`.

## Claim
R := rad(2N_EN_E'), a_p(E) = p + 1 − #E(𝔽_p). Input: ζ_L(s) ≠ 0 for Re s > 31/32 for every number field L.

**Theorem 1 (effective isogeny theorem).** If E, E'/ℚ are not ℚ-isogenous, some prime p ∤ 2N_EN_E' has
a_p(E) ≠ a_p(E') and p ≤ (a log R + b)^32, with (a, b) by the mod-2 representations:

| case | ρ̄_{E,2}, ρ̄_{E',2} | auxiliary field degree | (a, b) |
|---|---|---|---|
| (a) | traces differ somewhere | ≤ 36 | (3.5, 17.5) |
| (b) | isomorphic, image GL_2(𝔽_2) | ≤ 48 | (4.5, 23.5) |
| (c) | both reducible (rational 2-torsion on both) | ≤ 128 | (11.5, 68) |
| (d) | every pair (needed only for isomorphic images of order 3) | ≤ 192 | (17, 100) |

**Theorem 2 (Serre's constant).** For every non-CM E/ℚ and every prime ℓ > C(E), ρ̄_{E,ℓ} (and ρ_{E,ℓ} into GL_2(ℤ_ℓ))
is surjective, where

  C(E) = 2(11.5 log rad(2N_E) + 68)^16, and C(E) = 2(4.5 log rad(2N_E) + 23.5)^16 if ρ̄_{E,2} is surjective.

| N_E | rad(2N_E) ≤ | C(E) ≤ (all non-CM) | C(E) ≤ (ρ̄_{E,2} surjective) | GRH (Mayle–Wang) |
|---|---|---|---|---|
| 11 | 22 | 3.50·10^32 | 2.95·10^25 | 8740 |
| ≤ 10^6 | 2·10^6 | 1.72·10^38 | 2.99·10^31 | 19,747 |
| ≤ 10^100 | 2·10^100 | 1.84·10^55 | 5.30·10^48 | 228,398 |

(Rounded up from logs/cert_output.txt: log₁₀ C ≤ 32.5433, 38.2337, 55.2643 and 25.4688, 31.4748, 48.7242.)
Theorem 1 at the same R: p ≤ 7.40·10^69, 1.90·10^81, 2.28·10^115 (case (d)); 3.06·10^64, 7.34·10^75, 8.45·10^109
(cases (a)–(c)).

**Conditional on Liu's B_2** (ζ_L ≠ 0 for Re s > B_2 = 0.873382 for every L; Hexagon 2611.00082, unrefereed; our
audit Riemann/Deep/Hexagon2): the same proof with exponent 1/(1 − B_2) = 7.8978 gives p ≤ (67 log R + 360)^7.8978
(all pairs), (45 log R + 225)^7.8978 (cases (a)–(c)), and C(E) = 2(45 log rad(2N_E) + 225)^3.9489,
2(17 log rad(2N_E) + 74)^3.9489 if ρ̄_{E,2} is surjective: C(E) ≤ 2.61·10^10, 8.41·10^11, 1.59·10^16 (4.01·10^8,
1.58·10^10, 3.36·10^14) for the three rows above. This beats Kraus's bound once rad(N_E) > 10^10.5.

## Before
- **Unconditional, C_E.** Linear in the conductor or its radical: Kraus 1995, 68 rad(N_E)(1 + log log rad N_E)^{1/2};
  Cojocaru 2005, (4√6/3)N_E∏(1 + 1/p)^{1/2} + 1; Zywina 1508.07661: max(17, (p+1)/2) for every p with v_p(j_E) < 0, max(17, log of the denominator of j_E), and
  max(37, (2√3/3)N_E^{1/2}∏(½ + 1/(2p))^{1/2}) without multiplicative primes. Height bounds: Masser–Wüstholz 1993
  (c·h^γ), explicit by Pellarin, Lombardo, Le Fourn. None is polylogarithmic in N_E. Wang–Wei 2507.06913: small for
  100% of pairs ordered by conductor, not every curve.
- **Unconditional, isogeny prime.** Linear in the level (Sturm bound for the two newforms), or a power of R via
  unconditional Chebotarev (Kadiri–Wong Δ_L^310).
- **GRH.** Serre 1981 Thms 21–22; Mayle–Wang 2109.08656: p ≤ (482 log R + 2880)², C_E ≤ 964 log rad(2N_E) + 5760;
  Chen–Swidinsky 2501.00202: (124 log R + 561)², C_E ≤ 446 log rad(2N_E) + 2254 (two steps have gaps, below).

## Method
- **Algebra (Serre, Mayle–Wang, Chen–Swidinsky, reused).** The 2-adic deviation group of (ρ_{E,2}, ρ_{E',2}) is a
  finite quotient Gal(L/ℚ), L ⊂ ℚ(A[2^∞]) unramified outside R, with a conjugation-closed set C on which the traces
  differ. |G| ≤ 255 and divides 2^a·9, so ≤ 192. New: if both residual representations are reducible, the kernel of
  reduction M/2M → M_2(𝔽_2)² is nil, so G is a 2-group and ≤ 128. For isomorphic surjective residual representations
  Carayol's theorem gives the φ-quotient of order ≤ 48 (Chen–Swidinsky); if the mod-2 traces differ, L = ℚ(E[2], E'[2]),
  degree ≤ 36. Hensel–Serre: log Δ_L ≤ (n − 1) log R + n log n.
- **Analysis (replaces Bach–Sorenson's GRH Chebotarev).** chebexp's certified master inequality plus one new term:
  the primes p | m that are unramified in L cost at most n_L‖k‖ log m/x (Lemma 3.2, via the exact count
  |C_G(g)|/|H| of degree-one primes of L^{⟨g⟩} above p). So p ∤ R comes for free, without Mayle–Wang's L(√m), which
  would double the degree and the discriminant.
- **Certification.** Chains of intervals in ρ = log R ∈ [log 2, ∞), with x = (aρ_a + b)^32 fixed per interval,
  Λ ≤ Λ_max(ρ_b) and n_L ∈ {2, n_max} (affine), and a monotone tail; 33 + 51 + 26 + 80 intervals (and 53 + 57 + 47
  for the Liu-conditional chains).
- **Serre's constant.** Mayle–Wang's argument: for ℓ > 37 non-surjective, image ⊂ C_ns^+(ℓ), ε_ℓ = χ_D with D | N_E;
  Theorem 1 for (E, E_D), rad(2N_EN_{E_D}) = rad(2N_E), gives p_0 with χ_D(p_0) = −1, a_{p_0} ≠ 0, ℓ | a_{p_0}, so
  ℓ ≤ 2√p_0. The twist pair is in case (b) (ρ̄_{E,2} surjective), (c) (reducible) or a degree-6 case (image C_3).
  New: p_0 = ℓ is excluded (ordinary reduction at ℓ is incompatible with C_ns^+), a step Mayle–Wang skip.

## Certificate (one command)
```
cd Riemann/Deep/SerreNum && uv run --no-project --with python-flint --with mpmath --with numpy python3 cert.py > logs/cert_output.txt
```
- Ends with `ALL CERTIFIED (84 s)`; one core, peak RSS 49 MB (logs/cert_time.txt). Seven chains: U36, U48, U128, U192
  (unconditional) and L48, L128, L192 (Liu's B_2).
- Imports the Arb kernel code of ../ChebExp/cert.py (unchanged). New code: `lhs_interval_R`, `lhs_tail_R`,
  `algebra_checks` (GL_2(𝔽_2), |GL_2(ℤ/2^k)|, the caps 192/128/48/36, C_ns^+(ℓ) facts by enumeration).
- Smallest margin LHS/Φ(1): 0.00198 (U36, U48, U128), 0.00064 (U192), 0.00188–0.00200 (L chains). `design_R.py`
  (floats) made `chains.json`; not part of the proof.

## Checks
- **Sanity (PARI, logs/sanity_output.txt).** 318 curves of conductor ≤ 300 (small coefficients), 262 classes with
  distinct a_p, 34,191 certainly non-isogenous pairs: the least separating good prime is ≤ 31; log₁₀(p/(17 log R + 100)^32) ≤ −64.8. Twists of
  11a1, 37a1 and the 17- and 37-isogeny curves of Zywina's examples (conductors 14450, 1225): least p with χ_D(p) = −1,
  a_p ≠ 0 is ≤ 19 for every D | 2N_E. Their non-surjective primes 5, 17, 37 lie far below C(E) (10^32.5–10^34).
  Hensel–Serre holds on 765 fields ℚ(E[2], E'[2]), with equality for ℚ(√±2) (the bound is sharp).
- **GRH audit.** rg 'GRH|Riemann|conditional|assum' in both TeX sources, every hit read: GRH enters only through
  Bach–Sorenson (Mayle–Wang Thm 9/Cor. 10; Chen–Swidinsky Thm 2.4/Cor. 2.5/Prop. 2.6) (DERIVATION §6).
- **Reviewer** (pasted excerpts): Lemmas 3.1–3.2 correct; 2-group bound correct after naming the group kernel
  G ∩ (1 + I) (fixed); p_0 ≠ ℓ correct; the two Chen–Swidinsky gaps confirmed (the second with their definitions
  spelled out, done in DERIVATION §6).
- **Float cross-checks.** design_R.py (numpy, ChebExp/design.py kernels) predicted margin ≥ 0.003 on every interval;
  the Arb certificate finds ≥ 0.00064 (U192) and ≥ 0.00188 elsewhere. `xcheck.py` (mpmath, 40 digits, written
  independently of cert.py and ChebExp) re-evaluates the first interval and the tail of all seven chains: all positive,
  smallest 0.0030 (logs/xcheck_output.txt).

## Gaps found in the sources
1. Chen–Swidinsky Thm 1.3, first case: non-isomorphic mod-2 representations can have equal traces (one rational
   2-torsion point vs full 2-torsion), so their class C_2 is empty; the bound is unproved for those pairs.
2. Chen–Swidinsky Thm 1.4 (used for their C_E ≤ 446 log rad(2N_E) + 2254): for √d ∉ ℚ(E[2^∞]) and reducible ρ̄_{E,2}
   they claim α = β = 1, but α = 2 > β = 1 and their φ-set is empty.
3. Chen–Swidinsky take m = rad(N_EN_E'), which allows p = 2 where ρ_{E,2} is ramified.
4. Mayle–Wang's proof of Thm 1 does not exclude p_0 = ℓ (closed here).
Candidates for REFUTATIONS.md (minor, fixable gaps in GRH-conditional papers).

## Novelty
NOVELTY.md. arXiv listings Oct–Nov 2026, arXiv API (Zywina, Lemos, Le Fourn, Furio–Lombardo, Mayle, Chen–Swidinsky,
phrases, all papers citing Cureton–Thorner), Zenodo, HAL, our files: no unconditional polylogarithmic bound.
2507.06913 is unconditional only for 100% of pairs.

## Caveat
- **Inputs.** Cureton–Thorner's 31/32 (Lean-checked, paper unrefereed) via chebexp; Faltings, Serre–Tate, Carayol,
  Serre, Mazur, Bilu–Parent–Rebolledo, Balakrishnan et al. (qualitative, cited).
- **Size.** Exponents 32 and 16 make the numbers huge. C(E) beats Kraus only for rad(N_E) > 10^48.2 (10^40.5 if
  ρ̄_{E,2} is surjective) and Zywina's bounds for curves with one bad prime p only for p > 10^51; Theorem 1 beats a
  linear-in-the-level bound only for R > 10^119. The gain is the shape, polylogarithmic in the conductor without GRH.
- **Routine?** The proof is Mayle–Wang's with the GRH step replaced; a referee may call it a corollary of chebexp.
  The excluded-primes lemma, the 2-group bound and the case split are small but new.
- **Case (d)** (isomorphic residual images of order 3) still needs degree 192.

RH relevance: N11

## Files
- `REPORT.md` (this file), `DERIVATION.md` (every step and constant), `NOVELTY.md`.
- `cert.py` (certificate), `chains.json` (intervals and kernels), `design_R.py` (how the chains were found).
- `sanity.gp` (PARI sanity checks), `xcheck.py` (independent mpmath re-evaluation).
- `logs/`: cert_output.txt, cert_time.txt, sanity_output.txt, xcheck_output.txt, design_*.txt (float design runs, not part of the proof).
