import Mathlib

/-!
# Round 45 — "finite base + uniform tail": the tail is always the RH-equivalent

| family | base verified | best unconditional tail | verdict |
|--------|---------------|-------------------------|---------|
| Li λ_n ≥ 0 | n ≤ 10⁵ numerically | none (no eventual sign) | CIRCULAR |
| Jensen J^{d,n}, D_r(n) | d ≤ 3 all n (GORZ); each d for n ≥ N(d) | no uniformity in d | CIRCULAR |
| Weil windows | finite truncations, L ≤ log 15 claimed | none | CIRCULAR |
| Báez-Duarte c_k | ~10⁷ | O(k^{-1/2} e^{-c√log k}) | CIRCULAR |
| Robin | n ≤ 10¹⁰ | (e^γ + o(1)) log log n | CIRCULAR |
| Riesz R(x) | via zeros to 3·10¹² | x^{1/2} e^{-c(log x)^{3/5}…} | DEAD |
| Mertens | |M(x)| < √x to 10¹⁶ | x e^{-c(log x)^{3/5}…} | CIRCULAR |
| de Bruijn–Newman | — | Λ ≤ 0.2 | CIRCULAR |

Both refuted 2026 claims (Desogus, Gershon) had this shape and failed exactly at the tail lemma.
Lean: once a base is verified, the full statement is equivalent to the tail — the base buys nothing.
-/

namespace Fayed.R45

/-- With the base verified below N, the statement for all n is equivalent to the tail. -/
theorem tail_iff (P : ℕ → Prop) (N : ℕ) (hb : ∀ n, n < N → P n) :
    (∀ n, P n) ↔ ∀ n, N ≤ n → P n := by
  constructor
  · intro h n _
    exact h n
  · intro h n
    by_cases hn : n < N
    · exact hb n hn
    · exact h n (not_lt.mp hn)

end Fayed.R45
