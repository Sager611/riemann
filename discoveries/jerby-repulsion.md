# Jerby's Gram-point repulsion conjecture is false

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** Jerby conjectured (arXiv 2406.06548, conj:rep, TeX l. 313, and Conjecture 8.1, l. 777; 2310.14415 "Repulsion G-B-G", l. 266; restated as open in 2610.21952, l. 323): if g_n is a bad Gram point ((−1)^n Z(g_n) < 0, θ(g_n) = πn) whose neighbours g_{n±1} are good, then |Z′(g_n)| > 4|Z(g_n)|. This fails at n = 1,565,734,083 (n odd; g_n = 568,009,639.698978187581685…): Z(g_{n−1}) = 0.2763400589567841 (good), Z(g_n) = 18.39714778558964 (bad), Z(g_{n+1}) = 4.076959407739048 (good), Z′(g_n) ∈ [73.00385040031974 ± 4.4·10⁻¹⁵] < 4|Z(g_n)| ∈ [73.58859114235854 ± 3.0·10⁻¹⁵], i.e. μ = |Z′/Z|(g_n) ∈ [3.968215684906503 ± 2.3·10⁻¹⁶]. It fails again at n = 2,565,007,238 (n even; g_n = 906,089,969.7526938212…; Z(g_{n−1}) = −4.26923, Z(g_n) = −13.98639, Z(g_{n+1}) = −0.32154; μ ∈ [3.944872238693607 ± 4.9·10⁻¹⁶]) and at n = 3,315,629,581 (n odd; g_n = 1,155,454,823.1075850854…; Z(g_{n−1}) = 4.06689, Z(g_n) = 12.80566, Z(g_{n+1}) = 0.28173; |Z′(g_n)| = 50.23508562 < 4|Z(g_n)| = 51.22265784, μ = 3.922880048686967).

**Before.** Jerby supported the bound numerically to n ≈ 2.5·10^7 (2610.21952, l. 323). Our conjfals verified it for every isolated bad Gram point with n < 2.54·10^8 (smallest μ = 4.21553 at n = 69,784,845; RECORDS), with a normalised margin that shrank with height.

**Method.** jerbyfals: a Riemann–Siegel scanner with double-double phases relative to block centres (gram_win.c; against mpmath at 20 random Gram points per decade, relative μ error ≤ 3·10⁻⁷ up to 5·10^11). 64 log-uniform windows in t ∈ [10^8, 5·10^11] (29.0 million Gram points, 3.1 million isolated bad) gave no μ < 4 (minimum 4.2534). A scan targeted at the largest-|ζ| regions (proxy Σ_{p≤29} p^{−1/2} cos(t log p) + ½Σ_{p≤7} p^{−1} cos(2t log p)) of t ∈ [1.02·10^8, 1.0·10^9] (27.7 million Gram points, 2.5 million isolated bad) found all three counterexamples.

**Certificate.** Arb (python-flint, 256 bits; certify_cx.py): interval-Newton enclosures of g_{n−2}, …, g_{n+2} from θ(g_k) = kπ, with θ from the continuous log Γ branch; ζ(½ + ig_k) and its t-derivative as 2-term power series; Z = Re(e^{iθ}ζ). Radii are ≤ 4·10⁻²² for the g_k and ≤ 5·10⁻¹⁵ for Z and Z′. Arb also excludes Z = 0 on each block [g_{n−1}, g_{n+1}] (31,630 subintervals in the first case): every block is empty, a violation of Rosser's rule.

**Checks.** The orchestrator's own scripts. mpmath siegelz at 40 and 60 digits, with findroot for the Gram points, gives μ = 3.9682156849065032291 and the same neighbour values. Own Arb code (python-flint 0.9, 200 bits, acb_series ζ and log Γ) gives μ ∈ [3.9682156849065032290 ± 4.8·10⁻²⁰], for n = 2,565,007,238 μ ∈ [3.94487223869360749 ± 1.9·10⁻¹⁸] and for n = 3,315,629,581 μ ∈ [3.92288004868696710 ± 4.7·10⁻¹⁸], with the same good/bad patterns. The definitions (θ(g_n) = πn, Gram's law, the hypotheses of the conjecture) were read in the TeX of 2406.06548.

**Review.** Orchestrator's reviewer (statement, values and conventions pasted inline; it read certify_cx.py): CORRECT ("lgamma supplies the continuous branch normalized at θ(0) = 0; interval-Newton containment certifies the Gram roots, and the series correctly computes the t-derivative"), SIGNIFICANT ("refutes a published universal bound and establishes a barrier to arguments relying on it — not a disproof of RH").

**Novelty.**
- arXiv monthly listings Oct–Nov 2026 (lst.py on Gram, Jerby, repulsion, Hardy Z, viscosity, Rosser; 6–7 categories): only Jerby's 2610.21952, which restates the conjecture as open.
- arXiv API (Jerby; "bad Gram point"; Gram point and repulsion): Jerby's own papers only.
- OpenAlex: Jerby's Physica Scripta 2025 paper "A new discriminant for the Z-function and a Dyson-like repulsion for consecutive zeros" (doi 10.1088/1402-4896/adf3f1) is cited only by his own two later papers.
- Zenodo "Gram points": only Jerby's dataset 23390338 (250,000 Gram points).
- HAL: nothing on Gram-point repulsion.
- GXL evidence API (agent): only Jerby's own statements of the bound.
- Our files: only ConjFals (the bound holds below 2.54·10^8).
- No earlier counterexample found.

**Caveat.** All three counterexamples sit in zero-free Gram blocks (Rosser-rule violations), not in the "confined pairs" of Jerby's proposed mechanism. A refined bound, μ > 4 when the block holds two zeros, survives the data (minimum 4.81) but was not tested beyond it. RH is untouched. Lean: refuted in Lean 4 + Mathlib (JerbyLean.not_jerbyConj; named hypotheses: the Arb records and θ′ > 0 on [7, ∞); sorry-free, standard axioms; [LEAN.md](../Riemann/Deep/JerbyFals/LEAN.md)).

**RH relevance.** N8: Gram and Rosser counting; the extremes of μ are Rosser-rule violations. N13: Jerby's proposed repulsion mechanism fails as a universal bound.

**Files.** [JerbyFals report](../Riemann/Deep/JerbyFals/REPORT.md); in Riemann/Deep/JerbyFals/: certify_cx.py, certify_zerofree.py, gram_win.c, logs/certify_cx*.log, certify_zerofree*.log; Lean/ and [LEAN.md](../Riemann/Deep/JerbyFals/LEAN.md) (Lean refutation).
