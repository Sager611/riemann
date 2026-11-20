# PROBES — finite consequences of theorems on the cone, tested independently (alpmath)

Conventions: file paths relative to `.scratch/alp/lean/src/`; code in `code/`, logs in `logs/`; every run under `.scratch/memcap.sh 290`. "On cone": reached by name from the root (OUTLINE.md), unless marked "Main build only" (imported by Main.lean but outside the import closure of MainTheorems.lean).

## G1–G8 Green's parametrisation of cuspidal characters of GL_d(F_g)
Statement (`Development/AbelianVarietyTheory/PublishedTheorems/GaloisCohomology₂/GreenRegularCharacters.lean:13`, GreenRegularCharacters): for g prime, d > 0 and K a field with finrank_{F_g} K = d (so K = F_{g^d}) there is τ from characters θ of K^× to class functions on GL_d(F_g) with: θ regular (θ^{g^i} distinct, i < d) ⇒ τθ cuspidal (IsCusp, `Development/AbelianVarietyTheory/Definitions/NonGaloisDescent/LabelDatum.lean:47-65`: d > 0, irreducible, and the normalised restriction res a (d−a) χ 1 1 = g^{−a(d−a)}·Σ_{u ∈ U_a} χ(u) vanishes for 0 < a < d; zero sums are equivalent), τθ(1) = ∏_{i=1}^{d−1}(g^i − 1), and on y with irreducible characteristic polynomial and eigenvalue z ∈ K: τθ(y) = (−1)^{d−1} Σ_{i<d} θ(z^{g^i}); for regular θ, θ′: τθ = τθ′ iff θ′ = θ^{g^i}; every cuspidal irreducible character is some τθ.
Proof in archive: `Development/AbelianVarietyTheory/Proofs/InertiaEigenvalues/GreenCuspidalCharsByRegularCharsProof.lean:137` greenRegularCharacters_holds (general proof, not decide). On cone: `…/EveryCharacterOfLiteratureOfVWin.lean:399` → field greenRegularCharacters of P (`:694`) → GlobalTheory.GreenBridge.green_of_regular in every_character_of_theta_integrals_local_matching (`…/SiegelEisensteinSecondTermHypotheses.lean:69-72`) → Green (`LabelDatum.lean:108`). Literature: Green 1955; Deligne–Lusztig 1976 (KNOWN).
Test (`code/green_probe.py`, log `logs/green_probe.log`): character table of GL_d(F_g) by Burnside's algorithm (class-multiplication constants, simultaneous eigenvectors), orthonormality error ≤ 5·10⁻¹³; cuspidal = vanishing sums over every U_a; F_{g^d} built from an irreducible polynomial; elliptic classes = charpoly with d distinct roots in K, none in a proper subfield. Checked: number of cuspidal irreducibles = number of Frobenius orbits of regular θ; every cuspidal degree = ∏(g^i − 1); the multiset of predicted elliptic value vectors (−1)^{d−1}Σθ(z^{g^i}) = that of the cuspidal characters.

| # | group | order | classes | cuspidal = regular orbits | degree | elliptic values | result |
|---|---|---|---|---|---|---|---|
| G1 | GL₂(F₂) | 6 | 3 | 1 = 1 | 1 | match | PASS |
| G2 | GL₂(F₃) | 48 | 8 | 3 = 3 | 2 | match | PASS |
| G3 | GL₂(F₅) | 480 | 24 | 10 = 10 | 4 | match | PASS |
| G4 | GL₂(F₇) | 2016 | 48 | 21 = 21 | 6 | match | PASS |
| G5 | GL₂(F₁₁) | 13200 | 120 | 55 = 55 | 10 | match | PASS |
| G6 | GL₃(F₂) | 168 | 6 | 2 = 2 | 3 | match | PASS |
| G7 | GL₃(F₃) | 11232 | 24 | 8 = 8 | 16 | match | PASS |
| G8 | GL₄(F₂) | 20160 | 14 | 3 = 3 | 21 | match (sign (−1)³ = −1 tested) | PASS |

## Cone membership of the probe sources
`code/lit_reach.py` (log `logs/lit_reach.log`): all 196 declarations of `…/ThirdPlaceSeries/EveryCharacterOfLiteratureOfVWin.lean` are reachable by name from every_character_of_literature_of_vWin, i.e. every "Published" statement listed there is in the proof term of riemannHypothesis. Probes G, H, R, Q, W, L, S test such statements (name-level cone). Probes A, C, B test GroundTheory/Development lemmas in the file-level cone only: `code/uses.py` (`logs/uses.log`) finds no importer that names KroneckerLimit.full, ClassMean.full or the quadratic class-number formulas, so these may be dead code for RH; B1 and B4 are used one level up (ExistsIntModLevelPrimePowEqAddMul.lean; Primlevel.lean, ZeNorm.lean).

## H1–H4 Hilbert class field and CM (on cone)
HilbertClassFieldJ (`Development/AbelianVarietyTheory/PublishedTheorems/Divisors₂/ExistsEqOfNatCastMulMem.lean:66`): K imaginary quadratic, L₀ a lattice of O_K: K(j(L₀)) is the Hilbert class field, [K(j):K] = h_K, j(𝔞) = j(𝔟) ⟺ same class, the Galois conjugates of j(L₀) are the j-invariants of the classes. ComplexMultiplicationJSt (`…/PublishedFoundations/ComplexMultiplicationJPremises.lean:38`): Artin map acts by Shimura reciprocity. On cone: `…EveryCharacterOfLiteratureOfVWin.lean:225` hilbertClassFieldJ' := GlobalTheory.Laid.hcfj_holds (`Development/AbelianVarietyTheory/Proofs/DivisorProofs/HcfjHolds.lean:7`) := CompanionTheory.Hypotheses₂.hilbertClassFieldJ (`Development/AutomorphicSpectrum/Proofs/GlobalProofs/HilbertClassFieldJ.lean:7`, a theorem, derived from complexMultiplicationJ, `…/ComplexMultiplicationJacquetShalika.lean:91`); `:331` complex_multiplication_j_st_of_literature. KNOWN (Weber, Hasse, Deuring; Shimura).
Test `code/cm_probes.py`, `logs/cm_probes.log` (mpmath kleinj at 150 digits, reduced forms, sympy factorisation):
| # | consequence | range | result |
|---|---|---|---|
| H1 | H_d = ∏_C (X − j(τ_C)) ∈ Z[X] | 122 fundamental d ∈ [−399, −3], h ≤ 19 | PASS |
| H2 | j distinct on distinct classes | same 122 | PASS |
| H3 | H_d irreducible over Q (all 122); over K for h ∈ 2..4, abs(d) < 300 (43 fields); known values j(ρ)=0, j(i)=1728, j((1+√−7)/2) = −3375, j(√−2) = 8000, j((1+√−163)/2) = −640320³, H_−23 | | PASS |
| H4 | split p ∤ d·disc(H_d): H_d mod p has h roots in F_p ⟺ p represented by the principal form | 14 discriminants, 2957 pairs (d, p), p < 3000 | PASS |
First runs failed H1 (Re τ passed as a float, 1e-16 error; then an absolute tolerance at coefficients ~10¹²⁵) and H4 (my representation search used x ≥ 0 only): own bugs, fixed, logged in code comments.

## R1, Q1, W1, L1, S1 further on-cone statements
`code/cone_probes2.py`, `logs/cone_probes2.log`.
| # | statement (file:line; cone line in EveryCharacterOfLiteratureOfVWin.lean) | finite consequence tested | result |
|---|---|---|---|
| R1 | RadfordShufflePolynomial (`…/TwistedTraceFormula/RadfordShufflePolynomialHolds.lean:675`; :479): MvPolynomial on Lyndon words → shuffle algebra bijective | graded dimensions: ∏_{Lyndon w}(1 − t^{abs(w)})⁻¹ = 1/(1 − abs(A)·t), abs(A) = 2, 3, length ≤ 10 | PASS |
| Q1 | QuadraticClassFieldTheory.ProductFormula (`…/TwistedTraceFormula/QuadraticClassFieldTheory.lean:224,241`; :477): ∏_v η_v(t) = 1 | K = Q, L = Q(√a): Hilbert reciprocity ∏_v (a, b)_v = 1, 3000 random pairs abs(a), abs(b) ≤ 10⁵ | PASS |
| W1 | WeylBranching (`…/PublishedFoundations/FiniteIntertwining.lean:117,129`; :543): U(n+1) → U(n) restriction multiplicity 1 exactly on interlacing μ | s_λ(x₁..x_n, 1) = Σ_{μ interlacing λ} s_μ(x), 65 dominant λ, N = 2..4, entries ≤ 3 | PASS |
| L1 | HeckeLEntireAndNonzeroAtOne (`…/TwistedTraceFormula/HeckeLEntireAndNonzeroAtOne.lean:28`; :427): L(1, χ) ≠ 0, χ ≠ 1 unitary | K = Q: all 1042 nontrivial χ mod q ≤ 60, min abs L(1,χ) = 0.332 | PASS |
| S1 | SiegelKlingen (`…/Divisors₂/SiegelKlingenDeligneRibet.lean:40`; :231): L(1 − k) ∈ Q for ray functions on totally real F | ζ_F(1−k) for F = Q(√d), 9 fields, k = 2, 4, 6: Hurwitz at 1 − k = Bernoulli rationals to 2·10⁻²⁸ (formula consistency); independent S1′ (`code/review_fixes.py`): ζ(1−k) and L(1−k,χ) from the functional equations and values at s = k, 1.8·10⁻²⁷; ζ_F(−1) = 1/30, 1/12, 1/6 for d = 5, 8, 13 | PASS |

## M1, P1, U1 more on-cone statements
`code/cone_probes3.py`, `logs/cone_probes3.log`.
| # | statement (file:line; cone line in EveryCharacterOfLiteratureOfVWin.lean) | finite consequence tested | result |
|---|---|---|---|
| M1 | AbMulNSt (`Development/AbelianVarietyTheory/PublishedTheorems/Landherr₂/DimEqZero.lean:10`; :237): abelian B/k̄ of dimension g, char k ∤ n ⇒ #B[n] = n^{2g} | g = 1: y² = x³ + ax + b over F_p, p < 60, n = 3, 5, 7: ψ_n mod p separable of degree (n² − 1)/2 (so #E[n] = 1 + 2·(n²−1)/2 = n²), 172 cases | PASS |
| P1 | PartialZetaFunctionalEquation (`…/Divisors₂/PartialZetaFunctionalEquation.lean:11`; :227): d^{s/2}Γ_ℝ(s)^n Z_C(s) = d^{(1−s)/2}Γ_ℝ(1−s)^n Z_{C′}(1−s), C·C′ = [different] | F = Q and Q(√d), d = 5, 8, 13, 17, 29 (h = 1), 24 random complex s, rel err 1·10⁻³⁰ | PASS |
| U1 | ChevalleyUnits (`…/GaloisCohomology₂/ChevalleyUnits.lean:7`; :295): ∃ a coprime to b with u ≡ 1 mod a ⇒ u = v^m | Q(√D), D = 5, 8, 12, 13, 21; m = 2, 3, 5; b = 1, 2, 3, 5: such a ≤ 400 found in all 60 cases (all residues of ±ε^k mod a checked; first version missed k = 0, fixed) | PASS |

## N1, K1, J1 more on-cone statements
`code/cone_probes4.py`, `logs/cone_probes4.log`.
| # | statement (file:line; cone line in EveryCharacterOfLiteratureOfVWin.lean) | finite consequence tested | result |
|---|---|---|---|
| N1 | HasseNormSt (`Development/AbelianVarietyTheory/PublishedTheorems/Landherr₂/ExistsNormTensor.lean:12`; :417): x a local norm from F(√δ) at every place ⇒ x = a² − δb² in F | F = Q, all 893 pairs (δ, x) ∈ [−30, 30]² (δ non-square) that are local norms everywhere (Hilbert symbols): a rational solution found for each | PASS |
| K1 | ChebotarevDegreeOne (`…/GaloisCohomology₂/ChebotarevDegreeOneIffCheb.lean:12`; :323): every g ∈ Gal(L/F) is the Frobenius of a degree-one prime outside any finite S | L = Q(ζ_m), m ≤ 30, every a mod m: a prime p ≡ a, p > 500; L = splitting field of x³ − 2: primes 601, 503, 523 for the classes id, transposition, 3-cycle | PASS (weak: instances) |
| J1 | BrauerDegreeZero (`…/GaloisCohomology₂/BrauerDegreeZero.lean:12`, indFun `:9`; :275): χ_ρ − dim ρ = Σ c_i (Ind_{H_i} λ_i − Ind_{H_i} 1), c_i ∈ Z | S₃, D₄, Q₈, A₄, S₄: for every irreducible χ, χ − χ(1)·1 lies in the Z-span of Ind_H λ − Ind_H 1 over all subgroups and linear λ (exact lattice membership after Burnside character tables) | PASS |

## A1–A8 Kronecker limit / Chowla–Selberg folder (file-level cone)
`GroundTheory/HeightsAndTheta/Heights/ChowlaSelberg/`: each Statement is proved in a P*.lean file (e.g. `PCotangentSum.lean:17`, `PPowerIntegral.lean:28`, `PGammaZetaLimit.lean:33`, `PLogEta.lean:12`) and assembled in `Assembly.lean:32` KroneckerLimit.full, `:36` PhTau.full, `:39` ClassMean.full. Timeline item 9 (Chowla–Selberg). Test `code/analytic_probes.py`, `logs/analytic_probes.log` (mpmath 30–60 digits), A7 rerun `code/a7_rerun.py`, `logs/a7_rerun.log`.
| # | statement | test | result |
|---|---|---|---|
| A1 | CotangentSum (`CotangentSumStatement.lean:12`): Σ_n ((n+u)² + a²)⁻¹ = (π/a)(1 + 2 Re Σ_k e^{2πi(u+ia)k}) | 20 random (a, u), rel err 3·10⁻³¹ | PASS |
| A2 | PowerIntegral (`PowerIntegralStatement.lean:12`): ∫(t² + a²)^{−s} = √π Γ(s−½)/Γ(s)·a^{1−2s} | 12 random (a, s > ½), 1·10⁻³¹ | PASS |
| A3 | GammaZetaLimit (`GammaZetaLimitStatement.lean:12`): lim_{s→1⁺} 2√πΓ(s−½)/Γ(s) y^{1−s} ζ(2s−1) − π/(s−1) = 2π(γ − log 2) − π log y | 4 values of y at s − 1 = 10⁻²⁵, 1·10⁻²⁴ | PASS |
| A4 | LogEta (`LogEtaStatement.lean:12`): Σ_m (1/m) Re Σ_k q^{mk} = −π Im τ/12 − log abs(η(τ)) | 10 random τ, 2·10⁻³¹ | PASS |
| A5 | DeltaEta (`Statements.lean:47`): g₂³ − 27g₃² of Z + Zτ = (2π)¹² η²⁴ | 8 random τ with g₂, g₃ from the classical Eisenstein q-series, 5·10⁻³⁰; A5′ (`code/review_fixes.py`): direct lattice sums (box 150) for both G₄ and G₆ at τ = 0.31 + 1.07i agree with the q-series (6·10⁻⁶, 2·10⁻¹¹) and give Δ = (2π)¹²η²⁴ to 7·10⁻⁶ (box truncation) | PASS |
| A6 | KroneckerLimit (`Statements.lean:53`): lim_{s→1⁺} E(τ,s) − π/(s−1) = 2π(γ − log 2 − log(√y·abs(η)²)) | E by its Fourier–Bessel expansion, 6 random τ, s − 1 = 10⁻²⁰, 1·10⁻²⁰ | PASS |
| A7 | ClassMean (`Statements.lean:137`; classMean `:84`, logPeriod `:27`): mean over classes of log(4π Im τ_C abs(η(τ_C))⁴) = γ + log(2π/√p) − L′(1,χ)/L(1,χ) | prime discriminants −p, p ≡ 3 mod 4, p < 200 (24 fields): authoritative `code/a7_rerun.py` at 100 digits, max err 4.5·10⁻³⁰ | PASS (the run in analytic_probes.log used a float Re τ: 8·10⁻²⁰ against a 10⁻²⁰ threshold; script patched afterwards) |
| A8 | KroneckerSplit (`KroneckerSplitStatement.lean:12`): E(τ,s) = 2y^s ζ(2s) + 2Σ_{m≥1} y^s Σ_n ((n + m x)² + (m y)²)^{−s} | s = 2, τ = 0.3 + 1.1i: split = Fourier to 8·10⁻³¹; direct sum (box 300) to 3·10⁻⁵ | PASS |

## C1–C4 quadratic class number formulas, B1–B4 Bernoulli (file-level cone)
| # | statement | test | result |
|---|---|---|---|
| C1 | `GroundTheory/GlobalFields/AnalyticClassNumber/QuadraticConductor.lean:188`: h = −(w/(2abs(d)))Σ_{a mod abs(d)} χ(a)a, d < 0 | 911 fundamental d ∈ [−3000, −3], h by reduced forms, exact rationals | PASS |
| C2 | `…/AnalyticClassNumber/Quadratic.lean:254`: h = w√abs(d) L(1,χ)/(2π) | same 911, L(1,χ) by digamma, 9·10⁻³⁰ | PASS |
| C3 | `…/Quadratic.lean:267`: h·R = √d L(1,χ)/2, d > 0 | authoritative C3′ (`code/review_fixes.py`, `logs/review_fixes.log`): all 244 fundamental d ∈ [5, 800], ε by continued fractions (no skips), h independently = number of cycles of reduced indefinite forms (halved when N(ε) = +1): max abs(h_formula − h_cycles) = 2·10⁻³⁰. (First version: integrality only and a unit search that could skip fields silently; reviewer finding 1–2) | PASS |
| C4 | `Development/AbelianVarietyTheory/Proofs/CovolumeBounds/ExistsIdealInClassOfNormLeSqrtDiscr.lean:263`: h ≤ √abs(d)(1 + log √abs(d))^{n−1} | n = 2, the 911 imaginary fields | PASS |
| B1 | `Development/AbelianVarietyTheory/Proofs/FoundationProofs/PrimeMulBernoulliMem.lean:88`: p·B_{k+2} is p-integral | k ≤ 200, p < 100, exact | PASS |
| B2 | `Development/AutomorphicArithmetic/AutomorphicVarieties/GaloisLibrary/BernoulliNumbers/BLval.lean:8`: ζ(1−k, x) = −B_k(x)/k, k ≥ 1, 0 < x ≤ 1 | k = 1..12, 5 values of x, mpmath at 1 − k (may use Bernoulli internally), 3·10⁻³³; independent B2′ (`code/b2_fourier.py`): Hurwitz's formula ζ(1−k,x) = 2Γ(k)(2π)^{−k}Σ cos(πk/2 − 2πnx)/n^k summed exactly over residues with ζ(k, r/b) at positive k, k = 2..12: 4·10⁻⁴² (an nsum version in review_fixes.log failed numerically on the odd-k sine series) | PASS |
| B3 | `…/BernoulliNumbers/BernoulliPolyReflection.lean:6`: B_k(1−x) = (−1)^k B_k(x) | k ≤ 40, exact | PASS |
| B4 | `…/BernoulliNumbers/GeneralizedBernoulliChangeLevel.lean:22`: genB_M(χ↑M) = genB_N(χ)·∏_{p∣M}(1 − χ(p)p^{k−1}) (genB `MapGenBChar.lean:11`) | 775 cases: N ∈ {1,3,4,5,7,9,11}, all χ, M = Nm with m ∈ {1,2,3,4,6}, k = 1..5; 2·10⁻⁹ | PASS |

## Tally
39 probes, 0 counterexamples: 23 on the name-level cone (G1–G8, H1–H4, R1, Q1, W1, L1, S1, M1, P1, U1, N1, K1, J1), 16 in the file-level cone only (A1–A8, C1–C4, B1–B4). All test consequences of general Lean proofs (no decide/norm_num on these numbers). Largest run 124 MB RSS (`/usr/bin/time -l`, analytic_probes); every run under memcap 290.

## Review
Independent reviewer (subagent aca2af55, 20 Nov 18:28 UTC) on this file: no counterexample, no bug in the Green logic (cuspidality, elliptic classes, regular orbits, sign, matching), readings of A7 and of the positivity hypotheses confirmed; 1 blocking + 6 minor findings (C3 silent skips and integrality-only; B2/S1 independence; A5 only G₄ at τ = i; B2/B4 ranges; A7 script; Green wording). All fixed above: C3′, B2′, S1′, A5′ added (`code/review_fixes.py`, `code/b2_fourier.py`), ranges and wording corrected, A7 rerun designated authoritative.
