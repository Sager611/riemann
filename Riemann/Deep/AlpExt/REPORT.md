BLOCKED: the archive's host step needs every kept character to be a character of Gal(E/ℚ) for a CM host field E built inside ℚ(ζ_N) (KerChi, CharacterBoxes/TwistForm.lean:9-12; KeptInXp, TotientGaloisClasses/Analytic.lean:138; HostDer.cm, TotientGaloisClasses/ConjEMemHostFIff.lean:89), and the field cut out by the cubic class-group character ψ of ℚ(√−23), its Hilbert class field, embeds in no CM field, because its subfield ℚ(α), α³ = α + 1, has one real and one complex place.

# AlpExt: which steps of Alpöge's claimed RH/GRH proof are specific to Dirichlet characters, and where does it stop for ψ on ℚ(√−23)? (alpext, 20 Nov 2026, 21:13–21:44 UTC)

Archive .scratch/alp/lean/src/ (READ-ONLY; no Lean, no lake); paths relative to it. CB = Development/AutomorphicArithmetic/CharacterBoxes, TGC = Formalism/SharpRestriction/RamificationFrobenius/TotientGaloisClasses. Step table: SCOPE.md. Probes: PROBES.md. Target: L(s, ψ) = L(s, ρ) = ζ_K/ζ, K = ℚ(α) the cubic field of discriminant −23, ρ = Ind ψ the 2-dimensional representation of S₃ = Gal(H/ℚ), H = ℚ(√−23, α).

## Short answer
- Before the host nothing breaks: the level and the windowed Landau endgame have KNOWN analogues for L(s,ρ) (Hecke's continuation; P3: functional equation with root number +1, L(1,ρ) = 2π log α/√23 to 15 digits, real coefficients Λ(n)·tr ρ(Frob_n), positive combinations for 0 < w ≤ 1).
- The first step that fails is the host. The target must be a character of the Galois group of the CM field E on which the theta chain's hermitian spaces live (D7.host, Formalism/LocalRingCardinality/RhoDetPrimeEqMuTildeMulChiPow.lean:80-82). In the archive E is abelian, cut out in ℚ(ζ_N) by a subgroup of Dirichlet characters (fieldChi → AnyGrp.host, SCOPE §3). For ψ no CM E exists over any base: subfields of CM fields are totally real or CM (Lemma A), and ℚ(α) ⊂ H is neither (P0).
- Exact scope of the host step (Lemma B): a finite-order Hecke character ψ of an imaginary quadratic F passes iff ψ^c = ψ iff ψ = χ∘N_{F/ℚ}, whose L(s,ψ) = L(s,χ)L(s,χχ_F) the archive already covers. Every ring class character of order ≥ 3 fails; the first is ψ on ℚ(√−23) (−23 is the first fundamental discriminant whose class group has exponent > 2; P0, P0b).
- After the host, the steps transfer locally but all run on the CM host: the inert-prime divisibility keeps its exponent a/p + a/p³ + ⋯ at degree-1 and degree-2 supersingular primes of ℚ(α) (P2, finite check), Bost's slope inequality is general, the Galois representations are built for every CM field.
- Not shown: that no mechanism outside the archive can carry ψ. Three natural rewrites fail (§4); the third fails for the archive's mechanism only (host-class weights have no ρ-component, P1 C; twists by roots of unity change no valuation).

## 1. Where the target enters (SCOPE §0, §3)
- PChar = primitive DirichletCharacter ℂ N (TGC/PrimePairSizeCondition.lean:7-11). It occurs in 54 of 175,187 files, none in the theta-lift, slope, lattice or Galois-representation trees. Four doors: the level Level.cyc ψ (TGC/ExplicitFormulaErrorOfHostConditions.lean:105-121), the host (KerChi, TwistForm, KeptInXp; produced by fieldChi, CB/FieldChi.lean:20-26), the pair choice (PairGood ψ.N, CB/TwistForm.lean:22), VWin (IsChiField, CB/ExistsPageN.lean:19-25; same field by isChiField, CB/IsChiField.lean:5-6).
- The core sees the target only through the host it is handed: divisibility places are rational p with host.cls p inert (InertK, Formalism/LocalRingCardinality/DataDIVDefs.lean:45-46); its prime sums are weighted by Frobenius classes in GG = Gal(E/ℚ) (ArithData, Formalism/SizeWindowBrackets/FirstPrimeNeSecondPrime.lean:38-42).

## 2. Step verdicts for ψ (details, file:line and types (i)–(iv): SCOPE.md §1)
| Step | Uses | Verdict for ψ |
|---|---|---|
| S1 level, S2 windowed Landau endgame | (i) in the Lean; proof needs continuation, explicit formula, positive weights | NEEDS-NEW-INPUT, KNOWN: Hecke L-functions of F (entire, FE); not in Mathlib (P3) |
| S3–S5 host: E ⊂ ℚ(ζ_N) Galois, CM, target ∈ Gal(E/ℚ)^ | (i), (iii) | FAILS: ψ is no Dirichlet character, and its field lies in no CM field (P0) |
| S6 orthogonality ASK − MAIN = De | identity (iv) on abelian E0 × H; characters (i) | identity: KNOWN with Irr(G) incl. ρ (P1 B); ASK as written FAILS: partner = cconj + g + topPl t needs central complex conjugation for all class weights (P1 A) |
| S7 m-fold theta chain, inert-prime divisibility | (iv) on a CM host; auxiliary K = ℚ(√−d) (ii) | divisibility TRANSFERS locally (P2; finite check); the chain FAILS through S4 |
| S8 Bost slope inequality, Chowla–Selberg | (iv); periods of auxiliary K (ii) | TRANSFERS |
| S9 VWin and counts | IsChiField (i); Ferrero–Washington for abelian K (iii) | FAILS as stated; a non-abelian host would need μ = 0 for its cyclotomic ℤ₂-tower: UNKNOWN (open); moot after S4 |
| S10 PairGood/PairEx twists by (·/q), (·/d) | (i), (ii), (iii) | NEEDS-NEW-INPUT, KNOWN in kind (zero density for quadratic twists of a fixed GL(2) form); pair_share unopened |
| S11 ThinScales, exceptional zeros | (i)/(ii) | KNOWN (ψ·χ∘N is never real, no new Siegel zero) |
| S12 Galois representations (Row0b over every CM E) | (iv) | TRANSFERS for CM E; needs S4 |
| S13 imprimitive characters | (i) | TRANSFERS |

## 3. The decisive step (quotes, 13 lines)
```
CB/TwistForm.lean:9-12      def KerChi (ψ : PChar) {φ : Frame} (H : HostRaw φ) : Prop :=
                              ψ.N ∣ H.QS ∧ ∀ (a : ℕ) (h : a.Coprime H.condN), H.art (ZMod.unitOfCoprime a h) = 1 ↔
                                (a ≡ 1 [MOD 2 ^ (φ.e + 2)] ∨ a + 1 ≡ 0 [MOD 2 ^ (φ.e + 2)]) ∧ ψ.val a = 1 ∧
                                  jacobiSym (a : ℤ) H.q = 1 ∧ jacobiSym (a : ℤ) H.d = 1
TGC/ConjEMemHostFIff.lean:37  art : (ZMod (QS * q * d))ˣ →* φ.GG
TGC/ConjEMemHostFIff.lean:41  art_neg_one : art (-1) = φ.cc
TGC/ConjEMemHostFIff.lean:89  cm : IsCMField H.E
TGC/Analytic.lean:138       def KeptInXp (S : Finset PChar) : Prop := ∀ ψ ∈ S, ∃ θ ∈ φ.Xp, H.IsPrimOf θ ψ.χ
TGC/ExplicitFormulaErrorOfHostConditions.lean:146   kept : H.KeptInXp Λ.KEPT
Formalism/LocalRingCardinality/RhoDetPrimeEqMuTildeMulChiPow.lean:80-82   host : Host R
                                                                          sp : @Spaces R host.F host.E _ _ _ _ host.cmStar _
```
- Reading: the kept characters (all powers of the target) are characters of GG = Gal(E/ℚ), E is CM, and the theta chain's hermitian spaces are over host.E/host.F. For ψ, E would contain H (the field cut out by ψ), and H lies in no CM field.

## 4. Changing the base does not help (SCOPE §2; PROBES P0, P1)
- Lemma A: subfields of CM fields are totally real or CM. Lemma B: ψ's field lies in a CM field ⟺ ψ^c = ψ ⟺ ψ = χ∘N_{F/ℚ} (proofs in PROBES P0).
| Rewrite | First failure | Reason |
|---|---|---|
| R1: base ℚ, target the Artin ρ, non-abelian Galois host E ⊇ H | S4 (and S6) | E ⊇ H is not CM; askP reads w(cconj + i·frob p + t), a function of Frobenius classes for all w only if complex conjugation is central (P1 A: 100/300 non-constant on S₃ × C₂; abelian control 0) |
| R2: base F, class field theory of F for Kronecker–Weber, abelian E/F ⊇ H | S4 | E not CM: no totally definite unitary groups (at a complex place of E⁺, U(n) ≅ GL_n(ℂ)) |
| R3: abelian CM host E ⊇ F (e.g. ℚ(ζ₂₃)), chain twisted by ψ∘N_{E/F} | S5 | ψ∘N is no character of Gal(E/ℚ); the core weights primes by host.cls p, and E ∩ H ⊆ F makes every host-class function ρ-orthogonal (P1 C: exact; ℚ(ζ₂₃) data p ≤ 10⁶, max |mean tr ρ|·√n = 0.77); unit-valued twists change no valuation. A new trace-sensitive construction: NEEDS-NEW-INPUT, UNKNOWN |
- Inert-blindness: on primes inert in F, tr ρ(Frob) = 0 and tr ρ(Frob²) = 2, as for 1 ⊕ χ_F (P1 B: those class indicators lie in the span of linear characters); so an inert-prime sum relative to F itself cannot separate L(s,ρ) from ζ·L(s,χ₋₂₃). P2 shows the supersingular places of ψ's own CM curve are exactly these primes.

## 5. Probes (PROBES.md; all exact or certified, memcap peaks ≤ 14 MB)
- P0: one real root of x³ − x − 1 (signature (1,1)); complex conjugation = transposition, Z(S₃) = 1; first exponent > 2 class group at D = −23; Frobenius types for p ≤ 2·10⁵ (0 mismatches). P0b: blocked ring class orders, e.g. ℚ(√−11) f = 2 (h = 3), ℚ(i) f = 6, 7.
- P1: colour balance with ψ: partner well-definedness fails (A); orthogonality needs ρ (B: rank 4 < 6; exact 9/9 with Irr(G)); CM hosts cannot see ρ (C).
- P2: j₀ = −616750 − 1084125α − 820750α²; v_P(D^a℘(0)) = a·p/(p²−1) + O(log_p a) at every P above p inert in F (degrees 1 and 2; p = 5…43; a ≤ 398), none at split p.
- P3: L(s,ρ) = L(s, η(z)η(23z)): coefficients, FE, L(1,ρ), L'/L(1,ρ) = 1.01706472566, envelope 0.339; |G_ρ|/√y ≤ 0.185 on [10⁵, 10⁶].

## 6. Verdict options (honest tradeoffs)
| Option | For | Against |
|---|---|---|
| TRANSFERS | endgame, orthogonality (with Irr), divisibility, Bost, Galois reps have KNOWN analogues | the host has no analogue: no CM field carries ψ (Lemma A, P0); VWin's Ferrero–Washington input would become an open conjecture |
| BLOCKED (chosen) | exact obstruction at a named step with file:line; robust under the three rewrites; Lemma B gives the exact scope (only ψ = χ∘N pass) | proves failure of the archive's mechanism, not of every conceivable modification (reviewer) |
| UNCLEAR | supplyOfAdmOf/hostOfFieldOf bodies unopened; footprint argument is type-level | the failing requirement is in the types themselves (KerChi, KeptInXp, IsCMField), which no body can weaken |

## 7. Review
- One reviewer, verdict only, evidence inline (logs/reviewer.md): AGREE on the host failure (C2), Lemma A (C3), Lemma B (C4); overall 'BLOCKED at the host step' AGREE for the unchanged archive mechanism; CANNOT JUDGE exhaustiveness of the footprint (C1); DISAGREE with three wordings (partner: centrality needed uniformly over weights; R3 does not exclude new trace-sensitive constructions; P2 is a finite check; positivity needs w ≥ 0; NOGO wording: ψ's values lie in the CM field ℚ(ζ₃)). All fixed in SCOPE/PROBES and here.

## 8. Consequences for the loop (PostRH/OPEN.md item 1)
- If the archive is verified, it proves nothing new for Hecke L-functions of imaginary quadratic fields: ζ_K for disc −23 stays the smallest open Dedekind case. Alpöge expects the Extended RH "with more time" (.scratch/alp/remark.txt l.8-9); by this barrier that needs a host that is not a CM field over which definite unitary groups live, or a new way for the target to enter the divisibility count.
- Next falsification target (rank 1): theta lifts on U(n) relative to F/ℚ itself (F = ℚ(√−23) is CM) twisted by ψ: is any p-adic valuation of a lattice or period datum at a split prime 𝔭 | p of F sensitive to ψ(𝔭)? The archive's divisibility sits at inert places, where ψ is trivial (inert-blindness); unramified twists add no Gauss sums.

## Suggested NOGO line
Alpöge's route: targets must be characters of a CM host's Galois group; ℚ(√−23)'s cubic class-group character cuts out a field in no CM field.

## Files
- Riemann/Deep/AlpExt/REPORT.md (this), SCOPE.md (target entry, step table S1–S13 with file:line and types, rewrites R1–R3, host producer), PROBES.md (P0–P3, Lemmas A and B, limits).
- code/p0_cm_obstruction.py (P0, sympy), code/p0b_ring_class.py (P0b), code/p1_colour_balance_psi.py (P1; idea of AlpChain/code/colour_balance_exact.py, exact Fractions), code/p2_divisibility_d23.py (P2, python-flint, exact in ℚ(α)), code/p3_landau_rho.py (P3, mpmath).
- logs/p0_cm_obstruction.log, p0b_ring_class.log, p1_colour_balance_psi.log, p2_divisibility_d23.log, p2_divisibility_d23_k60.log, p3_landau_rho.log (memcap peaks at the end of each), logs/reviewer.md.
- Scratch (not for the repo): .scratch/alpext/STATE.md, .scratch/alpext/pchar_files.txt (the 54 PChar files).
