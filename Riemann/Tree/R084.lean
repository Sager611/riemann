/-!
# Round 84 — decay laws for λ_min(L), the smallest eigenvalue of the Weil form on windows of length L

| item | provable without RH? | verdict |
|------|----------------------|---------|
| upper bound λ_min(L) ≤ exp(−c e^L) via explicit test functions | yes | quantitative (true facts about decay) |
| lower bound λ_min(L) ≥ exp(−C e^L) for all L | ≡ RH (stronger) | CIRCULAR |
| Landau–Widom-type asymptotics for the Weil operator | yes, if from operator asymptotics | decay only |
| near-null eigenvector vs low zeros | finite windows: yes | descriptive |
| which zeros dominate at window L | partly | descriptive |
| off-line zero at height γ ⇒ λ_min < 0 once L ≳ log γ | yes (explicit formula + localized tests) | gives "RH up to γ ⇐ positivity up to L ≈ log γ" |
| convexity / log-concavity of λ_min in L | open; not automatic | curiosity |
| trade-off L₀(T) from RH verified to T = 3·10¹² | yes | quantitative (our 1.75 certificate) |

All unconditional facts about λ_min concern how fast it decays; positivity for every L remains exactly RH.
-/
