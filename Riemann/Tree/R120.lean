/-!
# Round 120 — triage of 8 recent Zenodo claims (Sep–Oct 2026; Zenodo API sorted by date; scouts downloaded and read them)

No discovery, no significant refutation. Every flaw below is visible on first reading.
* 23334982 Islam, "Formal Proof in Lean 4": no axiom or sorry, but it proves only the reformulation
  (LeastErasure (Rec ZetaZeros) → LeastErasure ZetaZeros) ↔ RiemannHypothesis. Rec projects every zero onto
  Re s = 1/2 (reg p = (0, p.2)). The needed bridge ("LookBack") is equivalent to RH by the file's own text.
  Lean-fidelity failure (R96 test 8).
* 23028070 Islam, "machine-checked, one bit": circular. The premise structure SelfGrounding already contains RH,
  and Faces carries the equivalences as fields.
* 23441946 Mehennaoui: conditional only (windowed determinants with real zeros converging to Ξ, then Rouché).
* 23161134 Caceres: RH-equivalent energy criteria; explicitly not a proof. A scout objected to
  ψ(x) − x ≪ x^Θ log² x with Θ = sup Re ρ, but that bound is standard (explicit formula with T = x; Ingham).
* 23120700 Saito, a claimed disproof: "π(x) − li(x) is unbounded both ways" contradicts nothing.
  Von Koch allows O(√x log x), and Littlewood's Ω± is consistent with RH.
* 23512195 Karaoğlan: the trace uses ζ(2s), so zeros map to Re s = 1/4 under RH, not 1/2.
* 23086975 Tsiokos: "conditional closure"; its premise is equivalent to RH (its own Proposition 6.2).
* 23102418 Bahtaji: a framework only; files restricted.
-/

namespace Fayed.R120
end Fayed.R120
