import Mathlib

/-!
# Round 39 — pinning the failing step of 8 recent claims

| claim | step | class |
|-------|------|-------|
| Velez, Zenodo 21173721 | "exact stabilization": positive cutoff forms ⇒ positive full form | finite → limit |
| Li-criterion, Zenodo 16900110 | phase alignment + O(n log² n) remainder ⇒ signs of λ_n | remainder swamps the zero terms |
| Yakaboylu, arXiv 2408.15135v11 | eqs. (17)–(18): biorthogonality ⇒ complete basis, resolution of identity | false general lemma (non-Riesz systems) |
| Bhattacharjee–Nandi–Frederick | Thm 3: "d > 0" read as contradicting "d = 0 ⇔ RH" | logic error (nb_no_contradiction) |
| Robin threshold, Zenodo 21869761 | §5 (18)–(19): partial Euler products ⇒ q₁ > ζ(3) | finite → limit (partial_lt_limit) |
| Q(A) form, Zenodo 18923112 | running-sum identity; GPU checks to 10¹³ | unclear / numerics |
| Abe, Cambridge Open Engage v25 | text not accessible | unclear |
| Gershon, Preprints 202604.1513 | Thm 1.1: strict log-concavity (TP2) vs Laguerre–Pólya | TP2 ⇏ real zeros (R038 core); Orca audit running |

Note: a scout claimed e^{-t⁴} has non-real transform zeros — false (Pólya 1927: all real). Our witness is e^{-t²}(1+t⁴/4).
-/

namespace Fayed.R39

/-- "d = 0 ⇔ RH" plus "RH fails" gives d ≠ 0: no contradiction can be extracted. -/
theorem nb_no_contradiction (d : ℝ) (RH : Prop) (h : d = 0 ↔ RH) (hn : ¬ RH) : d ≠ 0 :=
  fun hd => hn (h.mp hd)

/-- A strictly increasing sequence stays strictly below its limit: partial Euler products of
ζ(3) can never exceed ζ(3). -/
theorem partial_lt_limit {a : ℕ → ℝ} {L : ℝ} (hmono : StrictMono a)
    (hlim : Filter.Tendsto a Filter.atTop (nhds L)) (N : ℕ) : a N < L :=
  lt_of_lt_of_le (hmono (Nat.lt_succ_self N)) (hmono.monotone.ge_of_tendsto hlim (N + 1))

end Fayed.R39

#print axioms Fayed.R39.partial_lt_limit
