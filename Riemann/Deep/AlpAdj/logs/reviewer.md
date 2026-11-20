# Reviewer run (inline evidence E1-E9)

## Task

You are reviewing a VERDICT about a Lean 4 statement (no files needed; all evidence is inline; do NOT run Lean). Answer in <= 35 lines: verdict AGREE/DISAGREE + numbered concrete objections (each with what to fix). Be adversarial about mathematical faithfulness.

CLAIMED VERDICT: FAITHFUL-KNOWN. The archive's AdjointLiftingTheorem concludes AdZero, which is exactly the vanishing of the Bloch–Kato adjoint Selmer group H^1_f(L+, ad r) of a polarized automorphic Galois representation; its hypotheses include those of Thorne, arXiv:2207.04925, Theorem 5.1, so the statement is a special case of a published theorem (no residual hypothesis there). The brief that commissioned this check had assumed AdZero = an ORDINARY Selmer group H^1_ord and that published results need adequate/enormous RESIDUAL image.

E1 (AutomorphicPointGalois.lean:35-39) AdjointLiftingTheorem: for all primes l, CM fields L with no primitive l-th root of unity, s : HidaSetting l L Ω, V : s.Lev with V.S' = s.S, c with Large L c, 2<=n<=4, residual datum r with r.V0ok (residual characters chib_i : G_{L,S} -> k^x pairwise distinct, all unramified at an auxiliary place v0 and equal on its decomposition group), B universal ordinary determinant ring, x : B -> Qbar_l, rho : G0 -> GL_n(Qbar_l):  AutomorphicPointLocallyGeneric ι x rho -> Par x -> AdZero V rho.  AutomorphicPointLocallyGeneric = {cont; x.IsRep (det of point = det of rho); x.Arith h (psi_{v,j} locally algebraic, weights h strictly increasing); IsAutCsdA ι (rho∘pr); rho = 1 on inertia at v0; HasFlag at each v|l (rho|G_w upper triangular with diagonal psi_{v,1..n}); Free (exists g in G_w with chi_cyc(g)=1 and psi_i(g)^m != psi_j(g)^m for all i != j, m>0); TypeI; N3}.
E2 frame (ContinuousResP.lean:140-162): Γ = GS Ω S = Gal of maximal extension of L+ unramified outside S = {v|l} ∪ St ∪ {v0}; G0 = ker η = image of G_L; c = complex conjugation; ε with ε∘toGS = l-adic cyclotomic character of G_{L+}; G' := ker(ε|G0). HidaSetting also forces L/L+ unramified at all finite places and v|l, St, v0 split in L/L+.
E3 (OrdZero.lean:68-72) AdZero V rho := ∀ P R, IsForm rho P -> IsAdExt rho P R -> ∀ z ∈ Z1(Γ, R), IsBK l (R∘toGS) (z∘toGS) -> ∃ X, ∀ g, z g = R g X - X.  IsForm: rho(g)^T P rho(c g c^-1) = ε(g)^(1-n) P (polarization rho^c ≅ rho^∨ ε^(1-n)). IsAdExt (IsTanZero.lean:31-34): R g X = rho(g) X rho(g)^-1 for g ∈ G0, R(c) X = -(P^-1 X^T P).
E4 (BlochKatoSelmerVanishes.lean:135-140) IsBK l ρ c: away: ∀ finite w ∤ l ∃ place A over w, c restricted to inertia I_A is a coboundary; above: ∀ τ : F -> Om l ∃ j extending τ, c restricted to the decomposition group decJ F j lies in Loc.finiteCocycles. finiteCocycles (FiniteCohomology.lean:112-114) = {c ∈ Z1 : ∃ x ∈ V ⊗ B, ∀ g, c(g) ⊗ 1 = g·x - x}, B = (PP l).R.B an axiomatized period ring (PeriodData, Formalism/WeaklyAdmissibleFiltrations/Statement.lean:75-111: B^{G_K}=K0, injective semilinear φ commuting with G, injective ι : B -> BdR, BdR^{G_K}=K, exhaustive separated multiplicative G-stable filtration, t unit with g t = χ(g) t, φ t = p t, ι t ∈ Fil^1 \ Fil^2, fund: φ b = b and ι b ∈ Fil^0 => b ∈ Q_p); the instance comes from an in-archive existence theorem exists_locPrime_tied.
E5 TypeI V.frame rho := MatIrr G' (mat rho): every G'-stable subspace of Qbar_l^n is 0 or everything.
E6 IsAutCsdA ι r := AutCsd0 ι r ∧ (Continuous r -> WitA ι r); AutCsd0 := ∃ π cuspidal on GL_n(A_L), IsRegAlg π (regular algebraic infinity type), DualIsConj π (π^∨ ≅ π^c), AttC ι r π (r semisimple and, at cofinitely many places w, r unramified and charpoly(r(Frob_w^-1)) = ∏(X - q_w^((n-1)/2) α) over Satake parameters α of π_w).
E7 Thorne, arXiv:2207.04925, Theorem 5.1 (TeX l.1009): 'Let F be a CM number field, let n >= 2, and let (π, χ) be a regular algebraic, cuspidal, polarized automorphic representation of GL_n(A_F). Let ι : Qbar_p -> C be an isomorphism, and suppose that r_{π,ι}|G_{F(ζ_{p^∞})} is irreducible. Then H^1_f(F+, ad r_{π,ι}) = 0.' TeX l.400: 'we impose no condition on the residual representation (which might even be trivial)'. Here ad r is the Clozel–Harris–Taylor adjoint of the extension r : G_{F+} -> G_n (c acts on gl_n by X ↦ -A X^T A^-1).
E8 Argument: Γ-cocycles satisfying IsBK = classes of H^1(G_{L+}, ad~) that are unramified at all w ∤ l and crystalline (B-trivial) at w | l = H^1_f(L+, ad~) (BK classes are unramified outside l, hence inflated from G_{L+,S}; inflation is injective). AttC + Chebotarev + Brauer–Nesbitt => rho∘pr ≅ r_{π,ι}. RACSDC π is polarizable ((π, δ^n_{L/L+}) polarized). TypeI = irreducibility on G_{L(ζ_{l^∞})}. Hence AdZero is Thorne's conclusion.
E9 Side facts: the companion core DenseThresholdMainCase concludes OrdZero (cocycles upper triangular at w|l with diagonal vanishing on inertia, z = 0 on Iu ∩ ker chib at St, no condition at v0) for the same points without Par; the archive derives AdZero from OrdZero plus a local comparison. Test points: weight-one (Bellaïche–Dimitrov) excluded by Arith/IsRegAlg/Free; split-multiplicative exceptional zero excluded by Free (psi_1/psi_2 = chi_cyc); CM-by-L and endoscopic points excluded by TypeI and cuspidality.

Questions: (1) Is AdZero really H^1_f(L+, ad r) as in Thorne (sign of R(c), ad vs ad0, local conditions)? (2) Are Thorne's hypotheses implied? (3) Is any risk (axiomatized B, polarization normalization, Frobenius normalization, Γ = G_{L+,S}) serious enough to change the verdict? (4) Is FAITHFUL-KNOWN the right label, or should it be DISCREPANCY (since the commissioning brief described AdZero as ordinary), or INCOMPLETE?

## Answer

**DISAGREE — INCOMPLETE on the supplied evidence.** The main gap is the period ring, not residual image.

1. **Blocking: “B-trivial = crystalline” is not established.**
   E4 supplies period-ring axioms, not a comparison with Fontaine’s \(B_{\mathrm{cris}}\).
   `exists_locPrime_tied` establishes existence; its name does not establish that comparison.
   Even an embedding into genuine \(B_{\mathrm{dR}}\) alone gives only containment in the de Rham local kernel.
   **Fix:** prove, for the chosen instance and relevant adjoint representations,
   \(H^1_{f,B}\subseteq H^1_{f,B_{\mathrm{cris}}}\).
   This inclusion suffices to apply Thorne; equality is needed for “exactly \(H^1_f\).”

2. **The adjoint formula passes; the identification should include its quantifiers.**
   Taking \(A=P^{-1}\) gives precisely \(-AX^TA^{-1}\).
   The coefficient space is \(\mathfrak{gl}_n\), not \(\mathfrak{sl}_n\); its scalar summand carries the CM quadratic character.
   E3 quantifies over every compatible \(P,R\), whereas Thorne uses the polarized extension.
   **Fix:** identify one automorphic polarization; TypeI and Schur then make all compatible forms scalar multiples, hence give the same adjoint action.

3. **Chebotarev needs matching normalization, not merely matching terminology.**
   E6’s inverse Frobenius and Satake scaling must match the local-global compatibility convention for \(r_{\pi,\iota}\).
   **Fix:** supply that convention and the resulting identification, or directly identify the polarized adjoints.
   Scalar twists do not change the adjoint; normalization differences are not automatically counterexamples.

**Checks that pass:**
- TypeI gives the required characteristic-zero cyclotomic irreducibility: the image of \(\ker\chi_{\mathrm{cyc}}\) is exactly \(G'\). No adequate/enormous residual-image assumption is needed.
- Regular algebraic cuspidal conjugate-self-dual automorphy supplies the standard polarized setting.
- Away from \(l\), E4 gives the usual unramified local condition—not an ordinary condition.
- \(\Gamma=G_{L^+,S}\) is no obstruction: inflation is injective. Conversely, outside \(S\), inertia acts trivially, so locally coboundary cocycles vanish there and descend.
- The commissioning brief misidentified **AdZero**; E9’s **OrdZero** is a different assertion.

**Label:** record the brief’s discrepancy separately. **FAITHFUL-KNOWN** becomes justified after the period-ring bridge and normalization checks; no actual discrepancy with Thorne’s mathematical statement has been demonstrated.

# Second pass (period-ring bridge, Schur, normalization)

## Task

Second, short review (<= 20 lines). Earlier you (a reviewer) returned INCOMPLETE on a FAITHFUL-KNOWN verdict, blocking point: 'B-trivial = crystalline not established' for the Lean condition Loc.finiteCocycles = {c in Z1 : exists x in V ⊗ B, c(g) ⊗ 1 = g·x - x}, with B = (PP l).R.B, PP l := (exists_locPrime_tied l).choose. New evidence: exists_locPrime_tied : ∀ l, ∃ (P : LocPrime l) (AI : LocAinf l P), Nonempty (LocTied l P AI). LocTied (Formalism/LocalRingCardinality/AdicCompletionLocal/LocTied.lean:11-30) has fields: eB : P.R.B ≃+* Per.Bcris l; eD : P.R.BdR ≃+* Per.BdR l; zB : eB (algebraMap Q_l B z) = Per.zB l z for z ∈ Z_l; actB : ∀ g ∈ G_{Q_l}, eB (P.actQ g b) = Per.σB g (eB b); actD similarly; frob : eB ∘ φB = Per.φB ∘ eB; incl : the composite B -> BdR is the standard inclusion (Per.IsCrisDR). And Per (Formalism.Per, file .../InstIsAdicCompleteOCSpanSingletonSetCast.lean:33-72) defines Fl = tilt of O_C, Ainf = W(Fl) with Fontaine θ, AinfP = Ainf[1/p], BdRp = B_dR^+ (completion at ker θ), BdR, Acris0 = PD subalgebra (acrisZero), Acris, Bcris (Reformulated.AutPerF.Bcris l O_C). The B-action used in finiteCocycles is actC g = ρ(g) ⊗ (PP l).actQ(φ g). Also: (2) Schur: TypeI ⇒ ρ irreducible on G0 ⇒ every IsForm matrix is a scalar multiple of one P, and R(c)X = -P^-1 X^T P is scale-invariant; (3) FrobEig uses r(σ^-1) with σ arithmetic Frobenius and eigenvalues q^((n-1)/2)·Satake (Harris–Taylor normalization); a dual/twist mismatch would only replace π by π^∨, π^c or a twist and leave ad unchanged. Question: with this, is FAITHFUL-KNOWN (special case of Thorne arXiv:2207.04925 Thm 5.1: CM F, RACP (π,χ), r|G_{F(ζ_{p^∞})} irreducible ⇒ H^1_f(F+, ad r) = 0) now justified? Remaining blocking objections only, each with a fix; otherwise say 'AGREE'.

## Answer

AGREE
