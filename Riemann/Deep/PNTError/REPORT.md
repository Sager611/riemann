DISCOVERY: first sorry-free Lean proofs of ψ(x) = x + O(x^{7/8+ε}) and M(x) = Σ_{n≤x} μ(n) = O(x^{7/8+ε}) (from OpenAI's zero-free half-plane Re s > 7/8); new as a formalization only, mathematically the classical consequence.

Status: all three rungs proved. No sorry/admit/native_decide/axiom in Riemann/Deep/PNTError/.
#print axioms (check file .scratch/pnt/ax.lean): [propext, Classical.choice, Quot.sound] for every main theorem.
Build: lake build Riemann.Deep.PNTError.Psi (green; Psi.lean ~2.5 min, other files 10 s – 3.5 min each).

Main theorems (namespace Fayed.PNTError, file Psi.lean):
- logDeriv_bound_seven_eighths: ∀ε>0 ∃C, Re s ≥ 7/8+ε, |t| ≥ 2 ⇒ ‖ζ'/ζ(s)‖ ≤ C log|t|        (rung 1)
- psi_sub_self_le: ∀ε>0 ∃C ∀x≥2, |ψ(x) − x| ≤ C x^{7/8+ε}; psi_sub_self_isBigO (Landau form)    (rung 2)
- mertens_le: ∀ε>0 ∃C ∀x≥2, |M(x)| ≤ C x^{7/8+ε}; also inv_zeta_bound_seven_eighths: 1/ζ ≪ |t|^η (rung 3)
ψ is Mathlib's Chebyshev.psi; M is mertens x := Σ_{n∈Ioc 0 ⌊x⌋} μ(n), defined in Mertens.lean.

General versions (any zero-free half-plane Re s > σ₀, 1/2 ≤ σ₀ < 1; ZeroFree σ₀ in Basic.lean):
logDeriv_bound(_two), inv_zeta_bound, psi_error_of_zeroFree, mertens_error_of_zeroFree. Psi.lean only feeds in
Fayed.QRH.zeta_ne_zero_of_seven_eighths_lt_re.

Files: Basic.lean (Perron for ψ, contour shift), LogDeriv.lean (rung 1), General.lean (ψ assembly),
InvZeta.lean (1/ζ bound), Mertens.lean + MertensMain.lean (M(x)), Psi.lean (7/8 instances).

Method:
- Perron: ψ(x) − (1/2π)∫_{−T}^{T} F(4+iy)dy, F = −ζ'/ζ·x^s/s, at x with dist(x,ℤ) ≥ 1/2, via the existing kernel
  lemma perron_kernel (Inputs/PerronKernel.lean); error ≤ 4x⁴(x+1)/(πT)·π²/6. General x: ψ, M constant on [N,N+1).
- Shift to Re s = σ₀+ε/2 with ResidueTheoremOnRectangleWithSimplePole' (Z23 port); residue x at s=1 for ψ
  (Mathlib log_deriv_riemannZeta_add_inv_sub_bounded), none for 1/ζ (Mathlib inv_riemannZeta_eq_sub_mul). T = x⁵.
- Rung 1: Borel–Carathéodory (BorelCaratheodoryDeriv, LogOfAnalyticFunction from the Z23 PNT+ port) on the disc
  around 2+it; inputs ‖ζ(w)‖ ≤ 15|Im w| (Convexity/Partial.lean) and ‖ζ(2+it)‖ ≥ 1/2 (from L(1,s)L(μ,s)=1).
- 1/ζ ≪ |t|^η: Hadamard three circles around 3+it, proved from Mathlib's three-lines theorem applied to J∘exp;
  gives ‖log ζ‖ ≪ (log t)^θ₀ with θ₀ = log(3−σ₀−δ)/log(3−σ₀−δ/2) < 1.

Imports / clash check: QRH Bridge clashes with Z23 FromPNTPlus.ZetaBounds (duplicate
instCoeForallRealForallComplex_riemann via Z23 Fourier.lean), so ZetaBounds/PerronZeta are NOT used. Checked
co-import of QRH.Bridge + Inputs.PerronKernel + Z23 StrongPNTPrefix + Convexity.Partial: OK.

External check (GitHub API trees + code search, 2026-10-20): openai/math (classical x·exp(−c√log x) PNT only, no
power saving), facebookresearch/atlas-lean (qualitative PNT, Perron kernel), anthropics/formal-math zeta23
(MediumPNT), PNT+ and Mathlib: no ψ or M(x) power-saving theorem found. Code search is default-branch only.

Gaps / caveats: constants are not explicit-optimal; exponent 7/8+ε (not better). Not imported by Riemann.lean
(separate target, like QRH). Disk after build: 8.9 GB free.
