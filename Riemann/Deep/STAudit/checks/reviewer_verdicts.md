# Reviewer runs (2026-11-12, builtin reviewer, briefs inline, no tools)

## Run A — Theorem 1′ vs [R3] §3.5–§3.7 (brief: reviewer_A.md)
## Run B — Theorem MS given Theorem 1′ and Lemma L (brief: reviewer_B.md)
2/2 succeeded

=== Task 1: reviewer ===
VERDICT: CORRECT
ROUTINE
“Σ_{q≤E₀}q^{−1}σ_{−1}(q)𝓛(q,ψ)²”: the displayed §3.7 estimate requires precisely the twisted-sum bound and this weighted moment.
“ψ̃ := ψ·(·/q)”: these characters are real; conjugating the (B2) sum proves the same bound for b̄, with the same 𝓛_b.
“with cκ in place of κ”: legitimate; squarefree support forces (m,cκ)=1, and nonadmissible cκ contributes zero.
“V^{1/2}𝓛(q,ψ)”: multiplying the two (B2) bounds replaces V by 𝔙 without further assumptions.
“d(·)² by d(·)^{C_b}”: harmless using d(n)^A ≪_{A,ε} n^ε, including divisor factors coupled to 𝓛_b(q,ψ)².
“by partial summation over dyadic ranges”: the bound is ≪ C_bX′^{3δ}E^θlog²(2X′) ≪_{δ,C_b} X′^{4δ}E^θ, uniformly including θ=0.
“N ≤ X′^{1−4δ}”: choosing X₀ sufficiently large ensures 23N ≤ X′, so the required moment range is available.
“V ≤ X^{1/2}”: needed only for [R3]’s final simplification; Theorem 1′ instead retains √X′𝔙E^θ.
“(B2) is a hypothesis to be verified for each sequence”: this is a routine conditional extension, not a new cancellation theorem.

=== Task 2: reviewer ===
VERDICT: GAP — SIGNIFICANT
- “S₂ is treated identically” omits the changed outer coefficient: its bound is \(d_{2m+2}(s)s^{-1/2+\delta}\), not \(d_{2m+2}(s)s^{-1/2}\).
- “gives the extra factor N^{4δ}” accounts only for the inner \(k\)-sum; the outer squared coefficient also changes by \(s^{4\delta}\) relative to S₁.
- “and N = 1” exposes this omission for single-term blocks: \(N^{4\delta}=1\), while the outer coefficient still changes.
- “𝔖 := Σ_{s∈𝒮, s≤4P}d_{2m+2}(s)s^{−½}” supplies an elementary repair: retain these common Cauchy–Schwarz weights and insert \(s^{2\delta}\) inside the weighted sum.
- “σ = ½ − δ” then gives sieve factor \(N^{2\delta}\); since \(sN\le4P\le8P₀\), the total additional factor is at most \((8P₀)^{2\delta}\).
- “|ε_m|² ≤ C(q_m(|v|+3)^{2m+2})^{−2δ}” absorbs that repaired loss. Thus the claimed exponent survives; the theorem is not refuted.
- “k := the product of the p ∥ n with p ∉ 𝔅𝔓” gives the required unique coprime factorisation; the stated multiplicativity applies.
- “P₀ := X^{(m+1)/2+4δ}” exceeds the stated effective AFE length by a power \(X^\delta\); Lemma W(ii) correctly includes the single-term blocks.
- “X′ := p₀²X” ensures \(s\le X′\), \(N\le X′^{1-4\delta}\), and \(|v+\tau|\le X′^\delta\) for sufficiently large \(X\).
- “≤ 16X^{3(m+1)/4+50δ}” has correct exponent arithmetic; the remaining \(10\delta\) absorbs logarithms. The significant “Lindelöf on average” conclusion follows after the explicit S₂ repair.

## Run C — the S₂ repair (brief: the four-step repair of DERIVATION.md D4 with STLemP NOTE l.199–220 pasted)
VERDICT: GAP — ROUTINE
Step (3) overlooks \(N=1\) when \(2P<s\le4P\).
Use \(N\le4P/s\): then \(\sum_s w_2(s)N^{2\delta}\le(4P)^{2\delta}\mathfrak S\).
This recovers exactly the claimed \((4P)^{3\delta}\mathfrak S^2\) bound.
After \(|\varepsilon_m|^2\), the exponent is \(-\tfrac12(m+1)\delta+12\delta^2\).
It is nonpositive provided \(0<\delta\le(m+1)/24\); make that restriction explicit.
The quoted final \(50\delta\) comparison already needs still smaller \(\delta\).
Thus the repair works with mechanical clarifications; MS’s exponent is unchanged.
