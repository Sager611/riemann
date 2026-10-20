/-!
# Round 43 — certified numerics for ζ in Lean: feasibility

| topic | state of the art | plan | effort (lines) |
|-------|------------------|------|----------------|
| interval libraries | norm_num; LeanCert (2026, Lean 4.34) | port or reuse rational intervals | 5k–15k incl. infrastructure |
| ζ(1/2+it), t ≈ 14 | Euler–Maclaurin N = 16, m = 8: error ≲ 1e-16 | rational enclosures of √n, log n, sin, cos | 1k–2k |
| Hardy Z sign change | θ via Stirling remainder | Z(14.13), Z(14.14) opposite signs + IVT | 1.5k–3k |
| Turing's method to T₀ ≈ 100 | Trudgian 2011/2016 explicit S(t) bounds | Gram points + counting | 2k–5k |
| other provers | CoqInterval, Isabelle approximation, Flyspeck | reflection tactics, Taylor certificates | 2k–5k core |
| Arb certificates | Arb/FLINT balls | small verified checker in Lean | 2k–5k |
| rectangle counting | Rump, Arb | winding-number certificates | 5k–10k |
| zero-free box 0 < Im s < 14 | argument principle needed (Z alone misses off-line zeros) | ξ boundary enclosure | 2k–5k |

Shortcut via our Hardy files: Λ(1/2+it) is real (Riemann/Deep/Hardy/Repr.lean), so a certified sign change of
Λ(1/2+it) on [14.13, 14.14] plus the intermediate value theorem gives the first zero. Queued for an Orca deep dive.
-/
