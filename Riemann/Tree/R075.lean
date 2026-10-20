/-!
# Round 75 — pin the false statement of 8 claims (scout level)

| target | false / unproved statement | confidence |
|--------|----------------------------|------------|
| Ryazanov, arXiv 1808.10774 (rev. 2026) | eq. (4.1): φ*_n claimed n!-periodic "because each a_k is", and vanishing on every [n!k, n!k+1) | high → Orca audit |
| arXiv 2605.30767 "Musings on RH" | §2: neighbouring Re = 0 / Im = 0 contours cannot cross since a mixed Dirichlet/Neumann problem forces f ≡ 0 | high |
| Velez, Zenodo 21173721 | (not pinned in time) | — |
| Li-criterion Zenodo 16900110 | (not pinned in time) | — |
| Robin-threshold Zenodo 21869761 | eq. (19): q₁ = ∏_{p≤P}(1-p^{-3})^{-1} > ζ(3) — false for every P (Fayed.R39.partial_lt_limit) | medium |
| Q(A) Zenodo 18923112 | Thm 4.8, eq. (38): diagonal annihilation treated as controlling the intermediate band | high |
-/

/-! Correction (2026-10-20, Riemann/Deep/Pinned/Refute.lean): Robin-Zenodo eq. (19) refuted in Lean (not_robinEq19, P = 13). Q(A) eq. (38) is a correct identity; the error is in how it is used (diagonal cancellation does not bound the rest). Velez Lemma 7.1 follows if his Prop. 6.1 holds; the gap is elsewhere. Li-Zenodo 16900110: the deposited PDF has only definitions, no claim. -/
