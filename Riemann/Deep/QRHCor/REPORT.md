DISCOVERY: first Lean proofs, from OpenAI's ported 7/8 theorem, of Lindelöf right of 7/8 (ζ(σ+it) ≪ |t|^ε; L(s,χ) and 1/L(s,χ) ≪ (q(|t|+8))^ε for every χ ≠ 1) and of the effective bound L(1,χ) ≥ e^{-101}/log q. These are new as formalizations only; the mathematics is the classical consequence of a zero-free half-plane.

Status: all three rungs are proved, plus two extras (1b, and the real form of rung 3). No sorry, admit, native_decide or axiom in Riemann/Deep/QRHCor/.
#print axioms gives [propext, Classical.choice, Quot.sound] for all 13 main theorems (check file .scratch/qrhcor/ax.lean).
Build: lake build Riemann.Deep.QRHCor.{Dirichlet,Lindelof,LOne,LindelofL} is green. Each file takes 1–5 min and up to 3.6 GB.
LindelofL's strip lemma needs maxHeartbeats 10^6 and peaks at 6.5 GB. These are separate targets, not imported by Riemann.lean (same as QRH and PNTError).

Theorems (namespace Fayed.QRHCor):
Rung 1, Lindelof.lean
- lindelof_seven_eighths: σ > 7/8, ε > 0 ⇒ ∃C ∀t, 2 ≤ |t| → ‖ζ(σ+it)‖ ≤ C|t|^ε (the brief's statement).
- lindelof_seven_eighths_uniform: uniform on Re s ≥ 7/8+δ. zeta_lindelof: for any zero-free half-plane σ₀ ∈ [1/2,1).
Rung 1b, LindelofL.lean
- LFunction_lindelof: ∀δ,ε > 0 ∃C ∀q ∀χ ≠ 1 ∀s, Re s ≥ 7/8+δ ⇒ ‖L(s,χ)‖ ≤ C(q(|t|+8))^ε and ‖L(s,χ)⁻¹‖ ≤ C(q(|t|+8))^ε.
Rung 2, Dirichlet.lean
- LFunction_ne_zero_of_seven_eighths_lt_re: OpenAI's theorem, which was already in the port, restated.
- LFunction_real_ne_zero: no real zero β > 7/8, so no Siegel zeros. It assumes χ ≠ 1 or β ≠ 1, because Mathlib's value of L(1,χ₀) at the pole is an arbitrary placeholder.
- zeroFree_region_log: L(s,χ) ≠ 0 for Re s ≥ 1 − 1/(8 log q), for every t and every q ≥ 3.
- exists_absolute_real_zero_gap: OpenAI's SiegelZeros challenge statement. Here it is a 6-line corollary of 7/8 (c = log 3/8) and needs neither primitive nor real χ.
Rung 3, LOne.lean
- norm_LFunction_one_ge: χ ≠ 1 mod q ≥ 3 ⇒ ‖L(1,χ)‖ ≥ e^{-101}/log q.
- LFunction_one_re_ge and LFunction_one_lower_real_primitive_re: for real χ, L(1,χ) is real and ≥ e^{-101}/log q (the brief's statement).
- Tools: LFunction_eq_mul_Ichi (L(s,χ) = s∫₁^∞ A(t)t^{-s-1}dt on Re s > 0), norm_LFunction_le (‖L(s,χ)‖ ≤ q‖s‖/Re s), norm_LFunction_ge_of_one_lt_re.

Methods:
- Rungs 1 and 1b: PNTError's three-circles proof (InvZeta.lean), copied with a new last step: ‖ζ(s)‖ ≤ ‖ζ(3+it)‖·e^{‖J‖}.
  For L: Re J ≤ log(3q(|t|+8)) by norm_LFunction_le, and ‖L'/L‖ ≤ ΣΛ(n)/n² on Re s ≥ 2.
- Rung 3: Borel–Carathéodory for log L on the disc of radius 1/log q + 1/10 around 1 + 1/log q (zero-free, since Re ≥ 9/10).
  The real form uses the intermediate value theorem on [1, X].
- Continuation: Mathlib's LSeries_eq_mul_integral, with the integral written as a Mellin transform, then the identity theorem. Character partial sums are ≤ q.

External check (2026-10-20, GitHub API trees, contents and code search; code search covers the default branch only):
- openai/math: seven_eighths has no users outside its own files. SiegelZeros (a separate proof) gives the zero gap and c_ε·q^{-ε} ≤ L(1,χ).
  Elsewhere there is only an ineffective q^{-ε} (PrimeGaps) and 1/(64q) (TwoPoint). Nothing gives c/log q or Lindelöf right of 7/8.
- facebookresearch/atlas-lean has only definitions. PNT+ and Mathlib have nothing.

Caveats: the constant e^{-101} is crude; rungs 1b and 3 need χ ≠ 1; the 7/8 theorem is OpenAI's, trusted through the port.
Next: the prime number theorem in arithmetic progressions with error x^{7/8+ε}, uniform in q. It needs L'/L ≪ log(q|t|), which the same template gives.
