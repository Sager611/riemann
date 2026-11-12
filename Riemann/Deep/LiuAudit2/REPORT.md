AUDIT INCOMPLETE: no false or missing number-field step in G1, G3, G4, G5; the field-specific coefficient algebra of the returned families and theta coefficients (Lemmas 3.4, 3.5, 4.2, M1 unit return) holds at all 418 tested exterior primes of norm ≤ 3000 of Q(ζ12) (2 of the 420 had no small generator), with explicit allowed generators; open: the per-row analytic bounds and the bookkeeping of the M1/M2 inductions (Lemmas 10.6, 12.13), the KP theta input (Prop 4.1, cited) and the sieve chain of Lemma 13.3. The T12 flag is withdrawn. It came from a phase bug in my own T12 code: every CRT pair phase was conjugated. Corrected, and with Liu's smoothing and the full c·n³ coefficients up to norm 3·10^5, a_0 = Λ̄γ_2 shows no residue (max z² over 81 classes 5.45, family-wise p = 0.29; no growth; random phases). The γ_2 positive control shows the predicted X^{5/6} term in 77 of the 81 classes.

# Liu 2610.26835v2, Part II open items G1, G3, G4, G5 against OpenAI's Lean Q(√−3) model (liuaudit2, 12 Nov 2026)

Brief: .scratch/orca/liuaudit2-brief.md. Predecessor: ../LiuAudit/REPORT.md (audit 1). "l." = line of number_fields_v4.tex (arXiv e-print, 13,294 l.). "oai l." = line of openai/math preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/build/paper.tex (16,677 l.). Lean paths are relative to Riemann/Deep/QRH/OAI/NumberTheory/DirichletL/; they were read as text and not built. "CT Lean" = AxiomMath/DedekindZetaZeroFree, main branch, read through GitHub raw files and not built.

## Short answer
- **Nothing false found.** Every change from Q(√−3) to F ⊇ μ12 is handled in the text (§2). Each handling is either exact algebra, tested in §3, or standard analysis, re-derived in §4.
- **New tests over F = Q(ζ12)** (smallest case: degree 4, r_F = 2). Setup: S = {two complex places, p2 = (1+i), p3 = (√−3)}, U_S = ⟨ζ12, 1+ζ12, 1+i, √−3⟩, U_S/U_S^6 ≅ (Z/6)^4. B_S is an explicit isotropic complement, and no counterexample was found in 20,699 checks:
  - **Lemma 3.4** (allowed generators):
    - Allowed generators exist.
    - Sextic reciprocity is exact on 11,175 pairs. Control: 1,441 of 1,770 pairs of plain O_F generators fail.
    - χ_a(ε) = H_S(ε,a)^{−1}.
    - E_S is injective and Lagrangian.
  - **Lemma 3.5** (unit projection):
    - Orthogonality holds over all 1,296 unit classes.
    - 𝒢(U,εh) = ξ̄_ε(U)𝒢(U,h).
  - **Lemma 4.2** (theta coefficient characters):
    - γ_2^3 = −Λ with Λ = Λ_{Q(ω)}∘N_{F/Q(ω)}. Control: the conjugate fails at 404 of 418 primes.
    - Sextic duplication holds.
    - On allowed generators, γ_3 is a quadratic character of the B_S class. On plain generators it is no character.
    - G_0 is a character of order dividing 6.
    - a_0γ_1 = μG_0 and μγ_{−1} = a_0Ḡ_0 hold on 2,262 squarefree allowed n.
- **Consequences.** Given these inputs, two steps follow by algebra:
  - the M1 coefficient return ν ↦ νG_0ξ̄_{ε1} ↦ νξ̄_{ε1}ξ_{ε2} (l.3674–3723);
  - the closure identity of Lemma 10.6 (l.5007–5011).

  Audit 1's Lemma 12.10 test needed an inserted reciprocity factor. That factor is now 1 for allowed generators.
- **No-residue claim (T12/T13): consistent with Liu; the flag is withdrawn.**
  - Bug. The first T12 run conjugated every CRT pair phase: it used exponent 2k for χ_P(a)² = ζ3^k. That broke the Galois symmetry: class-0 sums of a_0 must be real, but only 11 of 26 norms were. Its numbers are void.
    - The fix restores 212 of 212 norms. The saved data were corrected exactly via γ_2^fix(c) = γ_2^bug(c)²/∏γ_2(p).
    - T1–T8 were not affected; they use sextic exponents directly.
  - T13, with Liu's own smoothing (W ∈ C_c^∞(1/2,1), l.2899–2903) and the full coefficient set γ_2(c)ϑ(cn³)/(√q_c q_n) of Lemma 5.1 (l.2634):
    - Liu's a_0 (ϑ = Λ̄): no growth. Mean slope of log(|S|/σ) is −0.03 (a pole would give +1/3). Max z² over 81 classes 5.45, family-wise p = 0.29. Phases random.
    - γ_2 positive control: the predicted Patterson term X^{5/6}. 77/81 classes have z² > 9, mean slope 0.27, best exponent 5/6, fixed phase.
    - The other Λ^k twists are clean, except one leak: the γ_2 term reaches k = +2 in class 7 through degree-2 primes, where Λ² = 1. It disappears when those primes are removed.
- **Prop 5.2.** Confirmed: the local Mellin factor 2^{2t+1}c^{−2t−2}Γ(5/6+t)Γ(7/6+t), the growth (1+|τ|)^{4r_Fσ}, and holomorphy on Re t > −5/6.
- **Lean.**
  - OpenAI's port covers Q(√−3) only.
  - CT Lean (Cureton–Thorner, every F ⊇ μ12) contains counterparts of Lemma 3.4 and of the Gauss correction of Lemma 4.2, but no M1 or M2.
  - Over F, M1 and M2 have no Lean model.
- **Open.** These are bounds with implied constants and Z^ε, so they cannot be falsified at computable sizes:
  - the per-row raw-volume and centered Mellin bounds;
  - the returned-family measures (δ_edge ≤ C_*(η+ξ+ε_0)) and the allowances behind the closure inequalities;
  - Lemma 13.3's sieve exponents;
  - Prop 4.1 (Kazhdan–Patterson, cited).

## 1. Correspondence table (open items)
Verdict column: ✓ correct (checked); ✓n numerically confirmed; P prose, not refereed line by line.

| # | Liu item (l.) | OpenAI 7/8 paper (oai l.) | Lean, Q(√−3) (file:line) | field-specific content | verdict |
|---|---|---|---|---|---|
| G3.1 | Prop 10.2, canonical inverse family bound (3940) | Lemma 17.2, canonical marked estimate (9296) | Descent/CanonicalRankExistence.lean:34 `canonical_marked_energy_exists`; Descent/CanonicalRankExistenceInduction.lean:25, 45 | allowed generators; unit classes; norm-sheet kernels (§2 rows 1–2) | P |
| G3.2 | Lemma 10.6, closure of inverse families (4977–5061) | §17.5, canonical reduction and its induction (9783–11601) | Inversion/InitialEnergyCallerFamily.lean:57 `normalized_child_uniform` (child = unit·primary generator, frequency h ↦ u^5h, u ∈ μ6), :87 `fresh_initial_split_fixed_family` | unit twist ξ̄_ε replaces h ↦ u^5h (3708–3723); coefficient identity (5007–5011) | ✓n coefficients; ✓n closure inequalities from the stated allowances (T11); P allowances (14η row, l.4746; ν(Ω) ≪ Z^{11η/2+ε}, l.5022–5027) |
| G3.3 | Thm 10.7, inverse moment M1 (5062–5078) | Lemma 17.1, marked inverse moment, unmarked case (9250) | Descent/InitializationEnergy.lean:73 `initial_energy_projection`, :120 `initial_weighted_count` | rows over R/Γ including unit classes; finite group ⊇ ⟨ν, G_0, Υ_S⟩ (5076–5078) | P (by Rem 10.9, l.5254–5266, it implies ζ_F ≠ 0 on Re s > 11/12) |
| G3.4 | Thm 11.2, marked M1 (5657) | Lemma 17.1 with slots (9250) | as G3.3; Inversion/InitialProfile.lean:256 `sourceDensity_moment_seminorm` | prime-polynomial factors | P |
| G1.1 | §12.7, two Poisson transforms (6555–7016) | §18.4, first transform (13114); §18.6, second transform (13596) | Moments/FirstPoisson.lean:116 `coprime_pair_poisson`; Moments/Extraction.lean:14–69 (Lemma 18.2); Moments/Cancellation.lean:54 `masked_rectangle_saving` (Lemma 18.3) | Gauss bridges for allowed generators | ✓n (audit 1; Lemma 12.10 now needs no inserted factor) |
| G1.2 | positive enlargements: Jensen over the pool, row ball (6776–6830, 6933, 6980, 7064, 7197) | §18.5, enlarging the Gauss-row norm (13411) | Moments/AmplificationFamilyEnergy.lean:21 `amplification_family_energy` (energy ≤ 16·Σ over the ray-character family) | pool from the prime ideal theorem over F (6787); the size of the finite character group replaces 16; ball without unit restrictions (7197) | ✓ positivity; P constants |
| G1.3 | Lemma 12.13, returned families (7407–7514); Rem 12.15 (7540–7573); App C.4 (12741–13100) | §18.8, completion of the induction (14779) | Energy/WidthInduction.lean:89 `certified_terminal`; Energy/CappedWidthInduction.lean:96 | weighted measures, δ_edge | P |
| G1.4 | per-row raw-volume bound Z^{a_0−b_2−r+2θ_N+ε_0} (7258–7327); Prop 12.4 (5997); Lemma 12.3 (5917) | §18.7, rows with inducing characters in Θ (14312) | Moments/ExceptionalPair.lean:54 `fixedInducingRow_of_product`; Moments/FiniteProfileExceptionalVolume.lean:18 | Kummer count of exceptional rows; finitely many unit classes (7180–7199) | ✓ exponents; ✓ Lemma 12.3 (standard); P volumes |
| G1.5 | Thm 12.16, M2 (7575–7753) | Lemma 18.1 (12531) | Energy/CertifiedExistence.lean:42 `terminal_certificate` (3/4 ≤ κ, 51/100 ≤ β, 2β−1 ≤ κ) | 13/18 ≤ κ (../Kappa); Θ ⊃ Υ_S (7596–7597), add-back Lemma 12.14 (7516–7539) | ✓ skeleton (audit 1); P |
| G4.1 | Prop 4.1, exceptional theta (2050–2351) | cubic theta from Eisenstein residues (§5, 1648ff) | Eisenstein/SquarefreeCoefficients.lean:122 `infinityCoefficient_squarefree`, :137 `…_squarefree_cube`; Eisenstein/PrimeGaussFactors.lean:504, 542 | KP over F (cited); S_f cutoff vectors; r_F factors K_{1/3} | citation; its W_p are used consistently in T7 |
| G4.2 | Lemma 4.2, coefficient characters (2352–2515) | coefficient formula (§5) | Eisenstein/PrimeGaussFactors.lean:627 `normalizedCubicUnitGauss_eq_gaussTwo`; SquarefreeCoefficients.lean:122 (factor `gaussTwo`) | Λ through N_{F/Q(ω)}; G' a Weil-index character on B_S; G_0 = G'ξ̄_4 | ✓n (T1, T2, T6, T7) |
| G4.3 | Prop 5.2, radial FE (2912–2943); multiplier (2812–2866) | completed reflection (Prop 5.1, 1676) | GaussianTheta.lean:269 `theta_functional_equation`; Hecke/Theta.lean:223 `completed_functional_equation` | Γ-ratio^{r_F}; finite factors e_v holomorphic on Re t > −δ, δ < 1/3 | ✓n (T9) + proof read |
| G4.4 | Lemma 13.3, mean squares (13.6)–(13.7) (8238–8502) | Lemma 15.1 (8111); Prop 15.2, Gram bound (8340) | MeanSquare/PolynomialHeight.lean:441 `completed_mean_square_norm_height`; MeanSquare/UniformPolynomialHeight.lean:460; CubicSieve/ (91 files), QuadraticSieve/ (36 files) | large sieves over F cited (GL Thm 1.1, BGL Thm 1.3); unit classes finite (8304); unit factor ξ_ε in the Gram correlation (8441–8455) | ✓n final exponent steps (T10); P sieve input |
| G4.5 | Prop 13.4, direct estimate (8503–8540) | Prop 15.3 (8564) | MeanSquare/PolynomialHeight.lean:456 `short_completed_average_norm_height` | — | ✓ given G4.4 (audit 1) |
| G5.1 | Lemma 13.9, detector Sobolev (9149–9313) | §8, detector (4201–4687) | Hecke/DetectorRowwise.lean:35 `rowwise_energy`; Hecke/DetectorFourier.lean:13 `uniform_fourier_moment`; Detector.lean:69 `simultaneous_saturation` | none (parameter box, O((log U)^2) dyads) | ✓ |
| G5.2 | §13.4–13.5, count and assembly (9226–9743) | §9, sextic row count (4687); §19 (15095) | Hecke/DetectorRowCount.lean:27 `inverse_count`, :42 `plain_count` | detector over F (Lemma 13.7, 8925) | ✓n formulas (audit 1, ../Hexagon2); P derivation |

The OpenAI paper cites no Lean names. Lean matches are by content: four read-only scouts proposed them, and all 18 declarations above were re-located by rg.

## 2. What changes from Q(√−3) to F ⊇ μ12, and where Liu handles it
| change | Q(√−3) device (Lean) | Liu over F (l.) | status |
|---|---|---|---|
| 1. unit group μ6 → U_S of rank \|S\|−1 (3 for Q(ζ12)) | primary generators, Descent/GeneratorTransport.lean:50 `generator_eq_unit_primary`; child frequency u^5h | allowed generators [a] ∈ B_S, unique mod Γ = U_S^6, h = εa_h (Lemma 3.4, 1908–1963; 1928–1931); unit twists ξ_ε go into one finite coefficient group (3708–3732, 5074–5078, 7596–7597); Θ ⊃ Υ_S is twist-invariant, while principal-only exclusion is not (2027–2033). M1: row unit sectors k = εk_0 before an unrestricted orbit norm (3892–3915). M2: unrestricted Γ-orbit row sums with no sector projector (6110–6117, 6688–6691); dual unit normalization (3537); the twist τ_1 ↦ τ_1ξ̄_ε on both columns (6342–6350, 6394–6397); the complete-support unit factor is 1 by reciprocity (6542–6549); the exceptional split is made after the twist, with the add-back of Lemma 12.14 (7608–7614) | ✓n T3–T8 |
| 2. lattice Z[ω] ⊂ C → R ⊂ F_S, sums over unit orbits | sectors and annuli (oai 2464–2638) | Poisson over unit orbits on G/Γ, G = F_S^1 (Prop 8.1, 3494–3600); compact dual support (Lemma 9.1, 3614); folded kernel K(T) = ∫_{F_S^1}Φ̂(tℓ_T)dt, \|𝒟^jK\| ≪ (1+\|log T\|)^{d_F}, d_F = [F:Q]/2+\|S_f\|−1 (3752–3770); origin kept without a factor N_S (2023–2025) | ✓ (audit 1; d_F re-derived: a simplex of r_F−1 complex logs plus \|S_f\| valuation shells) |
| 3. class number 1 → Cl(R) = 1 by enlarging S | Z[ω] is a PID | injectivity of U_S/U_S^6 → V_S; uniqueness of allowed generators (1947–1953) | ✓n T3 (rank 4 = \|S_f\|+[F:Q]/2) |
| 4. more roots of unity | primary sign conventions | χ_p(−1) = 1 (2498–2499, 3699); the Hilbert pairing is alternating (1940–1941) | ✓ |
| 5. degree 2 → 2r_F | one complex place | Γ-ratio^{r_F} (2821–2866); log power d_F (3753); Γ(s)^{r_F} in Lemma 12.6 (audit 1) | ✓n T9 |
| 6. different | (√−3) | S contains its support, so O_w is self-dual at exterior w (1808, 2416–2418) | ✓ |
| 7. archimedean places | one complex place | all complex, no sign conditions (2453–2454); K_{1/3} at each (2073–2076, 2806–2812); Λ has angular degree ±1 at each, which removes Patterson-type residues (2397–2400; Rem 10.9 item 2, 5290–5296) | ✓n T1 (this also fixes the embedding convention: −Λ, not −Λ̄); T12/T13 (corrected): no residue for a_0, Patterson term for γ_2 |
| 8. theta function | Eisenstein-residue construction, Eisenstein/ (82 files) | KP exceptional representation (Prop 4.1), cited: KP I.2.9, I.3.6, II.2.1, II.2.5; KPerr | citation (§4 G4) |
| 9. large sieves | CubicSieve/, QuadraticSieve/ | GL Thm 1.1, BGL Thm 1.3 (n = 3), published (Lemma 5.3, 2944–2988) | cited |
| 10. amplifier pool | prime ideal theorem in Q(ω) (oai 13450–13460) | prime ideal theorem over F (6787) | ✓ |

Candidate gaps (changes the text does not handle explicitly): none found.
- Only the weighted-measure and allowance bookkeeping (§4 G1, G3) is asserted rather than displayed. Its F-specific inputs are rows 1–2, and they are verified.
- Row 7's residue-free claim for a_0 is consistent with the corrected T12/T13 sums.

## 3. Tests (code/, logs/)
Every run is single-threaded and under perl alarm; time and memory are from /usr/bin/time -l.
- theta_chars_q12.py 3000: 29 s, 82 MB.
- theta_chars_q12.py 600: 12 s, 29 MB (second run).
- theta_chars_q12.py 1500 with seeds 1 and 2: 41 MB each. B_S is replaced by two other isotropic complements, f″_j + Σ_i s_ij e_i with s symmetric. T4–T8 pass again, and λ changes as it should.
- mellin_q12.py: 61 s, 35 MB.
- lemma133_exponents.py: 8 s, 61 MB.
- m1_closure.py: 4 s, 59 MB.
- theta_bias_q12.py 3e5: run twice, 631 s and 624 s, ≤ 241 MB maximum RSS (footprint 20 MB). Memory was watched with top during both runs. The data are saved to .scratch/liuaudit2/data/bias3e5.npz, which is not committed.
- theta_bias_analysis.py and theta_bias_twists.py: < 1 s, 43 MB.
- Follow-up runs:
  - fix_pairs.py: 2 s;
  - theta_smooth_q12.py and theta_smooth_fits.py: 1–3 s, ≤ 47 MB;
  - theta_bias_q12c.py at 3·10^4: 13 s; the C and numpy Gauss sums agree to 7.6e−15;
  - theta_bias_q12c.py at 10^6: stopped at the time limit after 63% of the C Gauss sums. Python ≤ 93 MB, C ≤ 82 MB.

Local square and cube classes at S are computed mod 8 (= p2^6) and mod 9 (= p3^4). This is valid because 1+p2^5 ⊂ F_{p2}^{×2}, since v_{p2}(4) = 4, and 1+p3^4 ⊂ F_{p3}^{×3}, since 3e/2 = 3 < 4. The dimensions found (5, 1, 1, 5) match the theory. H_S is computed on classes from tame symbols by Hilbert reciprocity.

| test | statement (l.) | instances | result |
|---|---|---|---|
| T1 | γ_2(p)^3 = −Λ(p), Λ = Λ_{Q(ω)}∘N_{F/Q(ω)}, π_+ ≡ 1 mod 3 (2364–2394) | 418 exterior primes of norm ≤ 3000 (2 of 420 had no small generator) | 0 failures, error ≤ 1.7e−14 |
| T2 | γ_1γ_2 = χ̄_p(4)γ_3γ_2^3 (2482–2486) | 418 | 0 failures |
| T3 | E_S injective, isotropic; H_S nondegenerate; H_S is a function of classes (1913–1953) | ranks 4/4 and 8/8; 2,070 further pairs | 0 mismatches |
| T3b | N_S^{−1}Σ_ε ξ_ε(b)ξ̄_ε(c) = 1_{b=c} (1975–1979) | all 1,296 unit classes, 576 pairs | 0 failures |
| T4 | χ_a(b) = χ_b(a) for allowed a, b (1921–1925) | 11,175 pairs | 0 failures (plain generators: 1,441/1,770 fail) |
| T5 | χ_a(ε) = H_S(ε,a)^{−1} (1921–1925) | ε ∈ {ζ12, 1+ζ12, 1+i, √−3, 4}; 2,090 checks | 0 failures |
| T6 | G'(𝔭) = γ_3(a) = λ([a]), λ a character of B_S of order 2; G_0 = χ̄(4)γ_3 of order dividing 6 (2437–2481) | fit on 4 primes, tested on 418 | 0 failures (plain generators: no character, rank 6 → 7) |
| T7 | a_0γ_1 = μG_0, μγ_{−1} = a_0Ḡ_0 on squarefree n (2356–2360, 2488–2499); CRT factor χ_{c1}(c2)^jχ_{c2}(c1)^j | 2,262 n with 2–3 primes; 84 CRT checks | 0 failures |
| T8 | 𝒢(U,εh) = ξ̄_ε(U)𝒢(U,h) (Lemma 3.5, 1990) | 352 | 0 failures |
| T9 | ∫_0^∞xK_{1/3}(c√x)x^t dx/x = 2^{2t+1}c^{−2t−2}Γ(5/6+t)Γ(7/6+t); \|ratio^{r_F}\| ~ \|τ\|^{4r_Fσ}; no pole on Re t > −5/6 (2846–2866) | 18 values of (c,t); σ ∈ {−0.3, 0.25, 1, 3}, τ up to 10^5 | PASS (quadrature error ≤ 3e−11 near the endpoint) |
| T10 | Lemma 13.3: identity (8355–8359); E(e) ≤ E(0)−e/6; E(0) ≤ M′+(d+5ℓ−1)_+/4, giving (13.6); the three branches of (13.7): zero frequency Q′/Y′, nonprincipal (Q′/Y′)P_a²/Y′, principal (Q′/Y′)P_a^{1/6} (8376–8500) | symbolic, plus 400,000 random admissible points | PASS (margin ≥ 0.021) |
| T11 | Lemma 10.6 / §10.4: the two exact identities (4903–4908) from the child centres (4612–4618), F = N+V, r = N−3ℓ; then F′−M′−Q′−z′ ≥ F−M−Q−z_0−15η and 4F′−3M′−6z′ ≥ 4F−3M−6z_0−44η (4918–4919) from the allowances M′ ≤ M̂+14η, Q′ ≤ Q+ℓ̂+t̂+ĝ−θ̂+η, F′ ≥ F̂, ĝ−θ̂ ≥ −η | symbolic, plus 200,000 random points | PASS (15η = 14η+η; 44η = 3·14η+2η) |
| T12 | Rem 10.9 item 2 (5290–5296), Prop 5.2 proof ("the theta integral is entire", 2892): sharp class sums over squarefree c, 81 classes of B_S^{(3)}, for c = γ_2 (type 0), a_0 = Λ̄γ_2, Λγ_2. Corrected run (the first run's CRT pair phases were conjugated: void) | 26,022 prime ideals, 69,796 squarefree c ≤ 3·10^5, 12 checkpoints | γ_2: energy 86.3, coherence 0.99, 77/81 classes > 3σ, 99.6% explained by c_bX^{5/6}. a_0: energy 0.61, coherence 0.20 (random 0.32), class 0 z = 1.0, none > 3σ. Λγ_2: same |
| T12b | the same data: Λ^kγ_2 and Λ^k for k = −3…3 | same | only k = 0 deviates; Λ^{±2,±3}γ_2 energy 0.9–1.2; pure Hecke sums Λ^k energy 0.13–0.24 |
| T13 | Liu's smoothed theta sums with the full coefficients b(cn³) = q_n^{1/2}γ_2(c)Λ(c)^kΛ(n)^{3k} (Lemma 5.1, l.2634; Prop 5.2, l.2899–2903); W a C_c^∞ bump on (1/2,1); 8 disjoint dyadic windows D ∈ [2,343, 299,977]; WLS fits S ≈ CD^θ against the exact random-walk variance; family-wise p over 81 classes. Pole prediction: constant-term power \|y\|^{1/3} gives a pole of D_m at s = 1/3, i.e. θ = 5/6 for b-sums, with fixed phase | 69,796 squarefree c × 17 cube cofactors n (q_n ≤ 61) | k = 0: z²_max 214.6, 77/81 classes > 9, mean slope 0.27 (pole: 1/3), Rbar = 1.00 in class 0, θ = 5/6 explains 96%. k = −1 (a_0): z²_max 5.45, p_fw 0.29, mean z² 0.66, slope −0.03, class 0 Rbar 0.25 (null fraction 0.61). k = ±1, ±3: p_fw ≥ 0.22. k = +2 class 7: 10.98 (p_fw 0.0014, 0.008 over 486 tests); it disappears without degree-2 primes (Λ² = 1 on them leaks the γ_2 term) |

The per-row raw-volume bounds and the theta mean squares (Lemma 13.3) were not run. They are bounds with implied constants and Z^ε, and their smallest cases are dominated by the diagonal. At computable sizes (rows of norm ≤ 3·10^3, lengths ≤ 3·10^2) no outcome could contradict them. The count behind the raw volume is uniform for a structural reason:
- a nonzero x ∈ cR has |x|_S ≥ q_c (product formula);
- U_S-invariance makes boxes of fixed volume a compact family;
- so #{x ∈ (a+cR)/Γ : q_x ≤ X} ≪ X/q_c + log X.

## 4. Referee verdicts on the open items
- **G3, M1.**
  - Unit return (l.3674–3723): correct, given Lemmas 3.4 and 4.2.
    - CRT exponent 3 per orientation: χ_{c1}(c2)^3χ_{c2}(c1)^3 = χ_{c2}(c1)^6 = 1.
    - ν ↦ νG_0 ↦ νG_0Ḡ_0 = ν.
    - ψ̄_1(ε_1) = ξ̄_{ε1}(u_1)ξ_{ε1}(u_2)·(frozen scalar), and likewise for ψ_2.
    - The frozen scalar ∏χ̄_p(ε_1)^{t_p} is constant on each unit-class slice. There are N_S slices, recombined by finite Cauchy–Schwarz (l.2014–2021).
  - Unit sectors (l.3892–3915): rows k = εk_0, ν(n)χ_n(k) = ν_ε(n)χ_n(k_0) with ν_ε = νξ_ε. Positive enlargement to the unrestricted Γ-orbit norm comes before orbit Poisson. Correct.
  - The two Gauss identities the transforms use are exactly those tested in T7:
    - first transform: a_0γ_1 = μG_0, giving c(u;y) with ν_1 = νG_0ξ̄_{ε1} (l.4459–4471);
    - second transform: μγ_{−1} = a_0Ḡ_0, giving L_2 with ν_2 = ν_1Ḡ_0ξ_{ε2} (l.4545–4551).
  - Lemma 10.6, coefficient identity (l.5007–5011): exact. χ_n(k′)χ_n(f′)^4 = χ_n(d_kk″)χ̄_n(d_2)χ_n(CJv′)^4, using χ^5 = χ̄, with zeros retained.
  - §10.1 kernel bound: correct (re-derived).
  - Closure inequalities (l.4918–4924 = 4985–4988): they follow exactly from the two identities (l.4903–4908) and the stated allowances (T11). M′ ≤ M−2(V+ℓ)+32η follows from M̂ ≤ M−2(V+ℓ̂)+3η, with ĵ ≤ 2ℓ̂+η, plus 14η + 2η.
  - Gap: the allowances themselves, the projection mass ν(Ω) ≪ Z^{11η/2+ε_1} (l.5022–5027) and the weighted integrals (App. C.2, l.12566) are asserted, not displayed.
    - Fixable: they are exponent bookkeeping on the same scheme as OpenAI §17.5.
    - The Lean argument for that scheme transfers once the unit twists are absorbed. Over Q(√−3) it is `normalized_child_uniform`; over F the twist is (l.3708–3723).
    - Not refereed here.
- **G1, M2.**
  - Jensen step (l.6805–6816; App. C l.12850–12860): correct.
    - #𝒜_h ≥ #𝒜/2 after removing O(log Z) primes.
    - |avg|² ≤ avg|·|²; each H = hp^6 has O(1) representations.
    - 𝒢(u,hp^6) = 𝒢(u,h), and 𝒢(p^iu,h) = χ_p(u)^{2i}𝒢(p^i,h)𝒢(u,h) needs allowed reciprocity (T4).
  - Row-ball enlargement: correct (all summands nonnegative; l.7197).
  - Unit handling in M2: correct, given reciprocity for allowed generators (T4).
    - No sector projector enters a row Poisson step (l.6110–6117).
    - The twist ξ̄_ε acts on both columns, so the common presentation is preserved (l.6342–6350).
    - The complete-support unit factor (χ_E(a)/χ_a(E))(χ_b(D)/χ_D(b))(χ_b(a)/χ_a(b)) = 1 (l.6542–6549).
  - Exceptional-row count (l.7180–7199): correct (Kummer).
  - Per-row exponent: correct. sup_{v≥0}{−2v/3 − (L−v)_+} = −2L/3 = −M/6 (l.7313–7316) and the F_2 table (l.7240–7250), both recomputed.
  - Lemma 12.3: correct. L^S(s,ϑ)·∏_{p|R_*}(1−ϑ(p)q_p^{−s}) is entire (nonprincipal ϑ) or has its only pole at s = 1. Shifting to Re s = 0^+ costs 2^{ω(R_*)}·Z^ε.
  - Gap: the returned-family measures and δ_edge ≤ C_*(η+ξ+ε_0) (Lemma 12.13 proof, l.7434–7514), which Remark 12.15 item 2 lists as "described in the text".
    - Fixable in principle: the Q(√−3) scheme is Lean-checked (Energy/WidthInduction.lean:89).
    - Its F-specific inputs (units, kernels, orbit Poisson) are verified here.
- **G4, theta.**
  - Prop 4.1 is a citation of Kazhdan–Patterson, not a proof. It was not refereed, and our Lean has no model over F.
    - Its exterior Whittaker values W_p(a(p^eu)) = χ_p(u)^{2e}·{q^{−k}, q^{−k−1/2}γ_2(p), 0} are consistent with the cocycle C_p(ϖ,u) = χ_p(u)^{−2} (l.2101–2104) and with the CRT form of γ_2 tested in T7.
  - No-residue mechanism (Rem 10.9 item 2; Lemma 5.1 proof l.2675–2684; Prop 5.2 proof l.2892): consistent with Liu at every tested size; the T12 flag is withdrawn.
    - Theory. Liu's constant-term Jacquet characters are pure norm characters |y_v|^{1/3} (the |y|^{2/3} constituent is absent). The weight L_υ(y)^{−1} has angular degree ±1, so ∫e^{−im_vθ}dθ = 0 kills them. Without the angular weight, the same constant term gives a pole of D_m(s) = Σa_m q^{−s} at s = 1/3. In the unit normalization b = q^{1/2}a this is a Patterson term C_bŴ(5/6)X^{5/6} with fixed phase, versus √X with drifting phase for noise.
    - Cause of the flag. The first T12 code conjugated every CRT pair phase: it used 2k instead of k for χ_P(a)² = ζ3^k.
      - Detected because class-0 sums of a_0 must be real: σ_11 maps class 0 to itself and a_0(σ_11𝔞) = conj(a_0(𝔞)). They were real for only 11 of 26 norms.
      - Corrected code: 26/26 norms at 3·10^4. Corrected data (exact identity γ_2^fix = γ_2^bug²/∏γ_2(p)): 212/212 norms.
      - The class-0 "stable phase" was this broken symmetry, not a residue.
    - Corrected sharp sums (T12): a_0 energy 0.61 at 3·10^5, below the random-walk value 1. γ_2 energy 86.3 with 77/81 classes above 3σ.
    - Liu's smoothing, full coefficients (T13): W ∈ C_c^∞(1/2,1); coefficients γ_2(c)ϑ(cn³)/(√q_c q_n) of Lemma 5.1 (l.2634), c squarefree, n arbitrary, with ϑ = Λ̄ (his ϑ = Λ̄Ḡ_0η with η = G_0; class restrictions are twists by Υ_S, which lies in his finite group).
      - a_0: max z²(θ = 5/6) over 81 classes is 5.45 (class 16), family-wise p = 0.29, mean z² 0.66.
      - Growth slope of log(|S|/σ) is −0.03 ± 0.14; a pole gives +1/3.
      - Class 0: z² = 0.10, Rbar = 0.25.
      - γ_2 control: slope 0.27 ± 0.07, 77/81 classes, Rbar = 1.00.
    - Not tested: the super-polynomial decay that Liu's identity implies for exact theta sums. With the degree-8 multiplier |m(t)| ~ |τ|^{8σ} and a C_c^∞ bump, the dual side decays only like exp(−cD^{1/16}), so random-walk-sized sums (energy ≈ 1) are expected at D ≤ 10^6. Correspondingly the γ_2 residuals stay at noise level.
    - Range. The C path (gauss_cubic.c, agreement 7.6e−15 with numpy) was started at X = 10^6 and stopped at the time limit after 49,052 of 78,256 degree-1 Gauss sums (63%). The decision rests on X ≤ 3·10^5: 69,796 squarefree c, 8 dyadic windows. Extending the run needs about 25 min of C time, or the cubic character from a precomputed discrete log.
  - Lemma 4.2: correct, proof checked.
    - Davenport–Hasse: J = −π_+^f.
    - Weil index: γ_3 = 𝔴_S(p^{−1})^{−1}; λ is a character on the isotropic B_S.
    - G' is trivial on α ≡ 1 mod^×𝔪.
    - All confirmed by T1, T2, T6, T7.
  - Prop 5.2: correct. Mellin factor and Stirling exponent by T9; the shift to −a_0 with a_0 < 1/3 is allowed by the finite factors (l.2829–2845).
  - Lemma 13.3, final exponent steps: correct (T10, plus the case analysis re-derived).
    - u = z_a gives E ≤ O/2+H ≤ M′.
    - Otherwise u ≥ v/2, so u+ℓ_b+½(T_d−v−3ℓ_b)_+ ≥ T_d/4, as stated (the true minimum is T_d/3).
    - Then E ≤ M′/2+z_a/2+(1+ℓ′)/4 ≤ M′+(d+5ℓ−1)_+/4, using z_a ≤ ℓ′ and Δ_H ≥ 0. Equality holds at z_a = ℓ′, so (13.6) is tight.
  - (13.7): correct. Zero frequency, nonprincipal and principal (Kummer count (P_a/q_C)^{1/6}) branches recomputed, with the unit factor ξ_ε retained (l.8441–8455).
  - T_d ≤ H−3d (l.8321–8322) also follows: 2A_0 ≤ O, O+H ≤ M′ = 1−ℓ−2d and z_a ≤ ℓ−d give H+2A_0+2z_a ≤ 1+ℓ−4d = N_*−3d.
  - Lemma 13.3, gap: the inputs to E(e_λ) are not refereed — the hybrid sieve (Lemma 11.1 identification, GL and BGL) and the powerful-part count (l.8323–8343). Audit 1 checked the consequence, Prop 13.4.
- **G5.**
  - Lemma 13.9: correct. Unit-interval Sobolev on a box of fixed dimension (eq. l.9163–9168). Derivatives of x^{−σ}, x^{it} add powers of log x on fixed annuli; the cost is polynomial in T_1.
  - §13.4–13.5: the formulas are confirmed (audit 1; ../Hexagon2). The prose derivation was not re-refereed.
- **No step is an error.**

## 5. Literature (12 Nov 2026; read-only scout plus own fetches)
- **Versions.**
  - arXiv 2610.26835: v1 19 Oct, v2 29 Oct; no v3.
  - Hexagon 2611.00082: v1 4 Nov only; v2 returns 404.
- **No erratum and no public doubt.**
  - QuasiRiemannTracker (GTimG) items updated since 10 Nov: #128, #129, #131, all closed. None concerns Part II.
- **CT Lean** (Cureton–Thorner, Hexagon 2610.00200v2; 5,528 .lean files) is the only Lean code over general F.
  - It has counterparts of Lemma 3.4:
    - Arithmetic/Ideals/AllowedGenerators.lean: `exists_allowed_generator` l.75, `allowed_generator_unique_mod_sixth` l.99, `allowed_generator_mul` l.140;
    - Arithmetic/Units/SUnitLagrangianComplements.lean:51 `exists_lagrangian_complement`.
  - It has counterparts of the Lemma 4.2 Gauss correction:
    - CharacterSums/Gauss/Sextic/WeilCorrection.lean:65 `gammaTwo_vartheta_tau_one` (γ_2·ϑ·τ_1·q^{−1/2} = −η(p));
    - CharacterSums/DavenportHasse/PrimeGaussIdentities.lean:185.
  - It also has unit-orbit norms (Theta/Unfolding/UnitOrbitNormSummability.lean) and Metaplectic/ (427 files) and Theta/ (510 files) trees.
  - It has no M1 or M2: their route avoids the general-field Möbius mean square (../LiuHecke/NOTE.md l.67).
  - Not checked: whether the `correction` field of its `GeneralArithmeticContext` is constructed for every F, and whether the trees build. GitHub code search needs authentication.

## 6. Dependencies
- ../LiuAudit/REPORT.md §5 is unchanged:
  - least prime represented by a form of discriminant D: |D|^{1.7498+ε};
  - Artin L zero-free on Re s > B_2;
  - Chebotarev (B);
  - the verdicts "sound given 2610.26835v2".

  All rest on all of Part II and stay "human, unrefereed", with the Lean-backed fallbacks listed there.
- The new tests remove one condition: Lemma 12.10 no longer depends on an untested Lemma 3.4. They do not change any entry.
- The withdrawn T12 flag does not change anything.

## 7. Reviewers
- Reviewer 1, on T1–T8 and the M1 unit return: GAP, on the wording only.
  - Sampled floating-point agreement confirms the identities; it does not prove them.
  - So item (ii) must read "the algebra of l.3674–3723 is correct, and its inputs are confirmed on every tested prime".
  - INCOMPLETE is justified.
- Reviewer 2, on the first (buggy) T12 data: INCONCLUSIVE. This is superseded: the excess was my phase bug. No reviewer was run on the corrected data, since nothing is reported FAILED.
  - "Open numerical flag, not decisive" is the right wording.
  - The random-walk null is heuristic, and nested checkpoints weaken the fits.
  - A nonzero angular type does not by itself prove the cancellation; the projection must be checked for the actual kernel.
  - Next test: the exact residue of Liu's theta Mellin transform.
- This report uses both reviewers' wording.

RH relevance: N7 (zero-free regions for Hecke families), N11 (load-bearing input of Liu's B♭ and B_2 for every finite-order Hecke L-function over every number field, and of our Artin, Chebotarev (B) and binary-form entries).

Files: code/theta_chars_q12.py (T1–T8, controls; optional second argument = shear seed), code/mellin_q12.py (T9), code/lemma133_exponents.py (T10), code/m1_closure.py (T11), code/theta_bias_q12.py (T12; second argument = npz path), code/theta_bias_analysis.py, code/theta_bias_twists.py (T12b), code/fix_pairs.py (phase correction), code/theta_smooth_q12.py and code/theta_smooth_fits.py (T13), code/theta_bias_q12c.py with code/gauss_cubic.c (C Gauss sums); logs/theta_chars_q12.log (bound 3000), _600.log, _shear1.log, _shear2.log, mellin_q12.log, lemma133_exponents.log, m1_closure.log, theta_bias_q12.log, theta_bias_q12_run2.log, theta_bias_q12_1e4.log, *_fixed.log (corrected T12, T12b, T13), fix_pairs.log, theta_bias_q12c_1e6.log (stopped run); first (buggy-phase) versions, void: theta_bias_q12.log, theta_bias_q12_run2.log, theta_bias_q12_1e4.log, theta_bias_analysis.log, theta_bias_twists.log, theta_smooth_q12_3e5.log, theta_smooth_fits.log (+ .err files with /usr/bin/time). Run: `uv run --no-project python3 code/theta_chars_q12.py 3000 [seed]`; `uv run --no-project --with mpmath python3 code/mellin_q12.py`; `uv run --no-project --with sympy python3 code/lemma133_exponents.py`.
