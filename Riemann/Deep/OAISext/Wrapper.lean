-- oaisext (Riemann/Deep/OAISext/REPORT.md): OpenAI's Lemma 18.1, z = 0 clause (ZeroAt), at κ = 1,
-- without the hypothesis 51/100 ≤ β_*, from the betaEff fork of the port. UNBUILT (see REPORT.md).
import Riemann.Deep.OAISext.OAI.NumberTheory.DirichletL.Energy.CertifiedExistence

namespace OAI

noncomputable section
open scoped Classical BigOperators SchwartzMap ContDiff
open Filter

namespace SevenEighths.OAISextWrapper
open HeckeFamily CenteredMomentEnergyState CenteredMomentEnergyBands
open CenteredMomentEnergyCappedWidthInduction CenteredMomentEnergyCertifiedExistence
open CenteredMomentNaturalFixedRaySource CenteredMomentFiniteProfileExceptional QuadraticInitialBound
local notation "O"=>HeckeFamily.O

instance neZero_top_ideal : NeZero (⊤:Ideal O) := ⟨by
  intro h
  have h1 : (1:O) ∈ (⊤:Ideal O) := Submodule.mem_top
  rw [h, Ideal.zero_eq_bot, Ideal.mem_bot] at h1
  exact one_ne_zero h1⟩

/-- terminal_certificate at κ = 1 (slot type Empty, ambient modulus and subgroup ⊤, slot weight 0)
composed with certified_terminal, zero clause: no hypothesis on the zero supremum. -/
theorem zeroAt_unconditional (a b radial Bmask L Mcap ε:ℝ)
    (ha:0<a)(haUpper:a≤1/4)(hb:1≤b)(hrad:0<radial)(hmask:0≤Bmask)(hMcap:0<Mcap)(hε:0<ε):
    ∃degree:ℕ,∃S:Finset (ℕ×ℕ),∀η₀:Character,∀Q:Ideal O,
      internalQ Q η₀≠0→internalQ Q η₀≠⊤→internalQ Q η₀≤Ideal.span {(72:O)}→
      ∃C:ℝ,0<C ∧ ∀ᶠZ:ℝ in atTop,1<Z ∧
        ZeroAt (internalQ Q η₀) a b radial Bmask L Mcap ε Z degree S C:=by
  have hcert:=terminal_certificate (α:=Empty) (⊤:Ideal O) ⊤ le_top (fun _=>0)
    1 1 a b radial Bmask L 0 0 Mcap 1 ε one_pos (by intro x hx; simp at hx) contDiff_const
    zero_le_one ha haUpper hb hrad hmask hMcap (by norm_num) hε
    HeckeZeroSupremum.fiftyOne_le_betaEff (by linarith [HeckeZeroSupremum.betaEff_le_one])
  obtain ⟨degree,S,h⟩:=certified_terminal (α:=Empty) (⊤:Ideal O) ⊤ le_top (fun _=>0)
    1 a b radial Bmask L 0 0 Mcap 1 ε hmask hMcap.le (by norm_num) hε hcert
  refine ⟨degree,S,fun η₀ Q hQ0 hQt hQ72=>?_⟩
  obtain ⟨Czero,_,hCzero,_,hev⟩:=h η₀ Q le_top hQ0 hQt hQ72
  exact ⟨Czero,hCzero,hev.mono fun Z hZ=>⟨hZ.1,hZ.2.1⟩⟩

/-- The same with ZeroAt unfolded: the statement uses only unmodified port definitions. -/
theorem zero_clause_unconditional (a b radial Bmask L Mcap ε:ℝ)
    (ha:0<a)(haUpper:a≤1/4)(hb:1≤b)(hrad:0<radial)(hmask:0≤Bmask)(hMcap:0<Mcap)(hε:0<ε):
    ∃degree:ℕ,∃S:Finset (ℕ×ℕ),∀η₀:Character,∀Q:Ideal O,
      internalQ Q η₀≠0→internalQ Q η₀≠⊤→internalQ Q η₀≤Ideal.span {(72:O)}→
      ∃C:ℝ,0<C ∧ ∀ᶠZ:ℝ in atTop,1<Z ∧
        ∀s:NaturalState Z Bmask radial,s.fixedModulus=internalQ Q η₀→s.width≤Mcap→
        ∀p:Profiles a b,∀t X₁ X₂:ℝ,0<X₁→0<X₂→X₁≤Z^L→X₂≤Z^L→
          s.plainEnergy p t X₁ X₂≤C*diagonalControl s.radial.profile*(p.control S)^2*
            (1+‖t‖)^degree*Z^(s.width+ε):=
  zeroAt_unconditional a b radial Bmask L Mcap ε ha haUpper hb hrad hmask hMcap hε

end SevenEighths.OAISextWrapper

end

end OAI
