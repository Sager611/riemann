# reviewer verdict (alpgal)

## Brief
You are an independent REVIEWER. Read-only. Do NOT run tools unless a single quick 'rg -n' on one named file is essential; you have < 2 minutes. Judge the VERDICT below from the evidence pasted here and reply in <= 25 lines: AGREE/DISAGREE with the verdict, then numbered concrete problems (overclaims, missed degeneracies, wrong inferences), each with a fix.

CONTEXT: Lean 4 archive claiming RH. Its Galois-representation input is Row0b (archive root .scratch/alp/lean/src/, GR = Development/AutomorphicArithmetic/Foundations/GaloisRepresentations). Question: is the formal notion 'Galois representation attached to π' the standard one (Chenevier–Harris, Shin, Caraiani, BLGGT), could a trivial object satisfy it, and does the route use more than it guarantees?

EVIDENCE (verbatim or close paraphrase, file:line):
E1 GR/TypeConditions/WeightsConjBis.lean:69-73: def Row0b : Prop := ∀ (E Ω) [NumberField E] [IsCMField E] [IsAlgClosure E Ω] (l prime) (ι : Qbar l ≃+* ℂ) (a : ℕ), 1 ≤ a → ∀ (π : CuspSub E a) (χ : HChar (maximalRealSubfield E)), (CG.mk π).IsRegAlg → PolarizedPair π χ → ∃ (r : G E Ω →* GL (Fin a) (Qbar l)) (wt : ℤ) (HT : (E →+* Qbar l) → Multiset ℤ), PropsB Rc π ι r wt HT. Here Rc = recT := fun h rc => IsRecT h rc (GR/IntegralStructures/Development.lean:18).
E2 PropsB (:54-67): cont : Continuous r; ss : IsSSRep r (every invariant subspace has invariant complement); away : ∀ w (hl : l ∉ w) A (h : IsPlaceOver A w) rc, Rc h rc → ∀ σ, π.IsLocComp w σ → AwayB (wdOf h hl ι r) (rc (algNorm σ)); weights : ∀ τ, (HT τ) in ℂ = infChar(π, ι∘τ) under z ↦ (a−1)/2 − z; weights_conj : HT(τ∘c) = (HT τ).map (wt − ·); above : ∀ P AI, Nonempty (LocTied l P AI) → ∀ w A h O Bw hj, AboveB …
E3 AwayB x y (:26-31): dim : x.lpDim = y.lpDim; tor : Tor y → y.mono0 → x = y; tame0 : IsQpBig w y.lpDim → TmW y → Sqf y → y.mono0 → x.mono0 ∧ ∀ g ∈ Weil, deg g = 0 → x.charPol g = y.charPol g. IsQpBig w n := E_w ≅ Q_g continuously with (n = 1 ∨ 16 n² < g).
E4 Tor y := ∀ s t ∈ Weil h, Triv y ⁅s,t⁆ (WeightsConjAlt.lean:27); Triv y g := y.charPol g = (X−1)^dim (…/LocalGlobalOfIwahori.lean:21-22). TmW: wild inertia acts unipotently.
E5 AboveB (:39-52): crys : π_w spherical → IsCrysAt ∧ IsdRAt ∧ HT jumps at τ = HT τ; nice : under ((IsQpBig ∧ TmW ∧ Sqf) ∨ Tor) of rc(algNorm σ) and mono0 → dR ∧ pst ∧ jumps ∧ AwayB (wdAt …) (rc (algNorm σ)) (WD of D_pst at l).
E6 G F Ω := Ω ≃ₐ[F] Ω (Formalism/ConjugateSelfDual/Statement.lean:6), Mathlib Krull topology; Qbar l := PadicAlgCl l (:8). No TopologicalSpace instance on PadicAlgCl/GL/AlgEquiv among the 133 TopologicalSpace instance headers in the archive's declaration index.
E7 wdOf (…/AutomorphicMembership/WeilDeligneParamMonodromyEqZeroIff.lean:352-355): if ∃ x, IsWDOf ι r x then choose else trivial WD rep diag 1; IsWDOf = ∃ tame char t, Frobenius φ (deg 1), monodromy N with IsMonodromy t r N and Frobenius-semisimplification a, x = [a] (:340-348).
E8 IsRecT h rc := ∃ a ha, LawsT a ha rc (Foundations/LocalLanglands/IsRecTOfIsRec.lean:76-77); LawsT fields: dim; one : rc(χ) = χ∘a (a = local Artin map); langlands (segment sums); cusp → irreducible; twists; central character; dual; inj; depth; type law. recTExists : ∀ E Ω A w h, ∃ rc, IsRecT h rc is a proved theorem (…/LocalCentralUnits/RecTExists.lean:20-67). IsLocalArtin.ord_eq : ord(a σ) = −deg σ (…/AutomorphicForms₂/HasExponentZeroIff.lean:149); deg 1 = arithmetic Frobenius.
E9 CuspSub E a (…/AutomorphicMembership/ExpoConjOfIsCsd.lean:209-216): nonzero algebraically irreducible subspace of cusp forms (automorphic: left GL_a(E)-invariant, open level, K-finite, Z-finite, moderate growth; cuspidal) stable under GL_a(A_f), K_∞ and Lie algebra. IsLocComp w σ (:375-377): ∃ nonzero GL_a(E_w)-equivariant σ → π. PolarizedPair (…/GaloisGeometry/DualIsConjByOne.lean:53-56): χ_v(−1) = (−1)^a at real v, π^∨ ≅ π^c ⊗ χ∘N. algNorm σ = σ ⊗ |det|^{(1−a)/2}.
E10 LocTied l P AI (Formalism/LocalRingCardinality/AdicCompletionLocal/LocTied.lean:11-37): ring isos P.R.B ≃ Per.Bcris l, P.R.BdR ≃ Per.BdR l compatible with Galois, φ, filtration, t, A_inf; Per.Bcris/BdR/Acris/Ainf are concrete constructions (W(O_C^♭), completions) in …/InstIsAdicCompleteOCSpanSingletonSetCast.lean:23-79.
E11 Consumers: SOfStar.lean (…/ThirdPlaceSeries/SOfStar.lean:28-50): isGalU_of_props_alt projects cont, ss, away, weights, weights_conj, above.nice (dR, jumps), above.crys.1; carU_of_star_alt also uses nice's AwayB at l. These feed the route's literature bundle (EveryCharacterOfLiteratureOfVWin.lean:286, 460, 618) and VWin (vwin_of_row h).
E12 Numerics (code/cm_check.py, logs/cm_check.log): E = Q(i), a = 1, π = ψ of y² = x³ − x, polarized with χ = |·|·χ_{−4}. From E1-E8: at unramified w, ι r(Frob_arith,w) = ψ_w(ϖ_w)^{-1}, i.e. r(Frob_geom,w) = primary generator. Checked for all odd p ≤ 500: a_p (two methods) = α + ᾱ with α ≡ 1 mod (2+2i) at 44 split p; a_p = 0 and tr Frob_(p) = a_p² − 2p = 2(−p) at 50 inert p; #E(F_{p²}) for p ≤ 60 agrees; 0 mismatches.

VERDICT TO JUDGE: WEAKER — the formal notion has the standard core (dimension a ≥ 1, continuity for the Krull and l-adic topologies, semisimplicity, Frobenius = Satake at unramified w ∤ l via AwayB.tor [unramified parameters are abelian with N = 0], full LGC at all w ∤ l when a = 1, HT weights from the infinitesimal character, crystalline + HT at w | l where π_w is spherical, dR/pst/WD at l in the 'nice' cases); absent vs the literature: full LGC at ramified w ∤ l for a ≥ 2 (only abelian N = 0 parameters, and inertia charpolys + N = 0 in tame cases with residue char > 16a²), de Rham at every w | l, the polarization r^c ≅ r^∨ ⊗ ε^{1−a} ⊗ r(χ), an explicit finite ramification set (but implied by tor at the cofinitely many unramified w). Not degenerate: a = 0 excluded; trivial r fails tor whenever some Satake parameter ≠ 1 because Rc is satisfiable (recTExists) and local components exist; junk defaults (wdOf, infChar) would also have to match nontrivial y. Downstream consumes only PropsB fields (E11), so the absent properties are not used.

## Reply
**DISAGREE as stated — INCOMPLETE.** “WEAKER” is plausible; several claimed guarantees remain unproved.

1. **Topology unverified.** E6 neither establishes the intended l-adic topology nor excludes imported instances. Fix: inspect the elaborated topology instances used by `Continuous r`.

2. **At-l guarantees are conditional.** Every `above` conclusion requires `Nonempty (LocTied l P AI)`; concrete period rings alone do not establish that witness. Fix: exhibit an admissible witness. Until then, `HT` is numerical data, not established Hodge–Tate weights of `r`.

3. **`wdOf` needs a realization theorem.** E7 permits a trivial fallback. Equality with that output alone cannot establish actual Frobenius compatibility or finite ramification of `r`. Fix: prove `IsWDOf` exists for each relevant restriction and identify the chosen Frobenius-semisimplification.

4. **Tame scope overstated.** The hypothesis includes `Sqf` and **already assumes** `y.mono0`. Moreover, `IsQpBig` requires \(E_w\cong\mathbf Q_g\), not merely large residue characteristic. Fix: retain the complete antecedent.

5. **Polarization character appears inverted.** From E9, using the rank-one reciprocity convention in E12, the expected formula is
   \(r^c\simeq r^\vee\otimes\varepsilon^{1-a}\otimes r(\chi)^{-1}|_{G_E}\).
   Fix: invert the character, or explicitly document a different convention for `r(χ)`.

6. **Triviality requires qualification.** Rank-one trivial automorphic data can legitimately have trivial attached `r`; that is not degeneracy. Also, `recTExists` supplies reciprocity maps, not local-component existence—E9 only defines the latter. Fix: distinguish legitimate trivial examples from vacuous attachment and provide local-component witnesses.

7. **`Tor` is not definitionally “abelian.”** It asserts unipotent characteristic polynomials on commutators. Fix: supply the finite-inertia/semisimplicity argument converting this into trivial commutator action before identifying its scope with abelian parameters.

8. **Downstream conclusion exceeds E11.** Immediate field projections do not establish that later literature interfaces faithfully express standard mathematical claims. Fix: audit bridge statements, proofs, and load-bearing axioms through `VWin`; limit the present conclusion to the displayed consumers.
