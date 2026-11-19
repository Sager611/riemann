AUDIT PASSED: every item of the jerby-repulsion claim holds as written — Jerby's conjecture is stated as the page says, and the three counterexamples (μ = 3.96822, 3.94487, 3.92288, each with good neighbours and a bad centre) are reproduced to 19 digits by independent code; five cosmetic corrections (a conjecture number in the Claim, four outside it; §4).

# JerbyAudit — line-by-line referee of discoveries/jerby-repulsion.md

Orca agent jerbyaudit, 19 Nov 2026, 20:16–20:50 UTC. One process at a time, single-threaded, peak 210 MB (memcap); no Lean, no commits. Proof chain with source lines: [FULLPROOF.md](FULLPROOF.md).

## 1. The claim, item by item

| # | Item (page, Claim) | Verdict | Evidence |
|---|---|---|---|
| 1 | Jerby conjectured it in 2406.06548 conj:rep, TeX l. 313, and "Conjecture 8.1", l. 777 | HOLDS | l. 313–317 and l. 777–781 (label con:8.1); abstract l. 169 ("≥ 4", "for any n ∈ ℤ"). Cosmetic: "8.1" is the TeX label; by the document's counters (theorem counter per section, conj shares it, l. 777 lies in §4 after one theorem) it prints as Conjecture 4.2 |
| 2 | 2310.14415 "Repulsion G-B-G", l. 266 | HOLDS | l. 266–268: |μ(g_n)| > 4 with μ = Z′/Z (l. 258–262) |
| 3 | restated as open in 2610.21952, l. 323 | HOLDS, with nuance | l. 323–326 reports μ > 4 at isolated bad points as an experimental property up to n ∼ 2.5·10^7 (not labelled a conjecture). Its Remark l. 616–664 notes that isolated badness does not guarantee two zeros (Rosser's rule fails, Trudgian) and phrases the stronger law for confined pairs (l. 653–654). The page's Caveat already says the counterexamples do not touch that version |
| 4 | the statement: g_n bad ((−1)^n Z(g_n) < 0, θ(g_n) = πn), g_{n±1} good ⇒ |Z′(g_n)| > 4|Z(g_n)| | HOLDS | definitions 2406 l. 164, 184, 204–211, 731–737, 769; 2310 l. 196–207. Jerby's typos (θ with "log t", l. 188; Z = ζe^{−iθ} in the abstract, l. 164) admit only the standard reading. Index and parity are confirmed by reproducing his first bad points g_126 and g_134 and his block {9807960, …, 9807963}. No height range or extra hypothesis |
| 5–7 | n = 1,565,734,083 (odd), g_n = 568,009,639.698978187581685…, Z(g_{n−1}) = 0.2763400589567841, Z(g_n) = 18.39714778558964, Z(g_{n+1}) = 4.076959407739048, Z′ ∈ [73.00385040031974 ± 4.4·10⁻¹⁵] < 4|Z| ∈ [73.58859114235854 ± 3.0·10⁻¹⁵], μ ∈ [3.968215684906503 ± 2.3·10⁻¹⁶] | HOLDS | independent: μ = 3.9682156849065032291 (inside the ball), Z′ = 73.003850400319735626 three ways (FULLPROOF S3) |
| 8 | n = 2,565,007,238 (even), g_n = 906,089,969.7526938212…, Z = −4.26923, −13.98639, −0.32154, μ ∈ [3.944872238693607 ± 4.9·10⁻¹⁶] | HOLDS | μ = 3.9448722386936074881 |
| 9 | n = 3,315,629,581 (odd), g_n = 1,155,454,823.1075850854…, Z = 4.06689, 12.80566, 0.28173, |Z′| = 50.23508562 < 4|Z| = 51.22265784, μ = 3.922880048686967 | HOLDS | μ = 3.9228800486869671047; Z′ = −50.2350856201 (sign irrelevant) |

Margins: 4 − μ = 0.0318, 0.0551, 0.0771, equivalently 4|Z| − |Z′| = 0.585, 0.771, 0.988. Numerical errors are below 10⁻¹⁵, and the discrepancies of Jerby's own approximations are below 0.002 in μ.

## 2. Certificate audit (certify_cx.py, certify_zerofree.py): every step holds

- θ is Im acb_lgamma(1/4 + it/2) − (t/2) log π, the analytic log Γ: continuous with θ(0) = 0 and no branch slip. My own Stirling series and mpmath siegeltheta put all 15 Gram points within the printed balls (≤ 3.3·10⁻²¹). A 2πj slip would shift the index by 2j, which keeps the parity, so it could not matter in any case.
- θ′ = Re ψ(1/4 + it/2)/2 − (log π)/2 is correct. θ′ > 0 on [7, ∞) because Re ψ(1/4 + iy) increases in y and θ′(7) ∈ [0.05358983564 ± 2.2·10⁻¹²] (Arb, zf_mv.py). So the interval-Newton root is THE Gram point g_k.
- Interval Newton N(X) ⊂ X with X = [m ± 10⁻⁴⁵] gives existence, uniqueness and the enclosure. The printed radii, up to 4.33·10⁻²¹, are decimal rounding from str(30); the computed radius is below 10⁻⁴⁵.
- In the series [1/2 + it, i] → ζ, the coefficient z1 = iζ′ = dζ/dt. Z′ = Re(e^{iθ}(iθ′ζ + dζ/dt)), and the θ′ term contributes −θ′ Im(e^{iθ}ζ) = 0, so it cannot be mishandled. Evaluating on the Gram ball covers the true g_k. Good/bad and μ < 4 are certified Arb comparisons.
- Zero-free cover: it starts 10⁻¹⁵ below g_{n−1} and ends 10⁻¹⁵ above g_{n+1} (decimal rounding of the endpoints ≤ 5·10⁻¹⁶), bisects at exact binary points (no gaps), and accepts a piece only when the enclosure of Z on a superset excludes 0. It holds.
- Re-certified by a different algorithm: mean-value form Z(c) + [−r, r]·Z′([a, b]), own Gram enclosures, no mpmath (zf_mv.py: 414, 473 and 512 contiguous pieces, 232 s, 14 MB). Negative control: at n = 6708 the cover stops at the Lehmer zero 7005.0628661751.
- Sampled own-RS check (zerocheck.py): no sign change on any block; min |Z| = |Z(g_{n±1})| = 0.2763, 0.3215, 0.2817, rising inward. The two missing zeros sit in the adjacent Gram intervals: zero counts [1, 3, 0, 0, 1, 0], [1, 1, 0, 0, 2, 2] and [1, 1, 0, 0, 3, 1] on j = n−3 … n+2. These are length-2 Rosser violations, as the page says.
- Lean (read-only): JerbyConj is conj:rep word for word, over Mathlib ζ and Γ, with θ = ∫₀ᵗθ′. The hypotheses are exactly the certified records plus ThetaDerivPos (true, see above). The axioms are propext, Classical.choice and Quot.sound, with no sorry or axiom in the 696 lines (FULLPROOF S5).

## 3. Reviewer (most delicate step: Arb records ⇒ rigorous counterexample)

A reviewer subagent received the code lines, formulas and outputs inline. Verdict CORRECT: "Interval-Newton inclusion N ⊆ X, with derivative bounded away from zero, proves existence in X and enclosure in N … The analytic log-Γ branch, real on the positive axis, fixes the standard θ normalization; no unknown 2πj shift remains … The series coefficient is iζ′(s) = dζ/dt … Outward-rounded evaluation on N encloses the values at its true Gram point … |Z′(g_n)| < 73.004 < 73.588 < 4|Z(g_n)|."

## 4. Corrections (cosmetic, outside the claim; not edited: discoveries/ and JerbyFals/ are read-only for this agent)

1. Page, Certificate: "Radii are ≤ 4·10⁻²² for the g_k" is false as printed: certify_cx.log prints 4.27·10⁻²² and certify_cx3.log up to 4.33·10⁻²¹. Fix: "printed radii ≤ 4.4·10⁻²¹ (decimal rounding; the Newton enclosures have radius < 10⁻⁴⁵)".
2. Page, Method: the targeted scan "of t ∈ [1.02·10^8, 1.0·10^9] (27.7 million Gram points, 2.5 million isolated bad) found all three counterexamples". The third (t = 1.155·10^9) came from scan R3, t ∈ [1.0023·10^9, 1.2023·10^9] (REPORT §5). Fix: "t ∈ [1.02·10^8, 1.20·10^9] (33.95 million Gram points, 3.02 million isolated bad)".
3. Page, Claim: "Conjecture 8.1, l. 777" is the label con:8.1, which prints as Conjecture 4.2 (§4; the theorem counter is shared). Same in JerbyFals/LEAN.md. Orchestrator check (arXiv PDF of 2406.06548, 19 Nov): the printed numbers are Conjecture 1.1 (l. 313) and Conjecture 7.1 (l. 777, §7 of 9 sections), not 4.2; the page and LEAN.md now say 7.1.
4. Page, Before: "verified it for every isolated bad Gram point with n < 2.54·10^8" was double-precision screening, with only the minimum Arb-certified (RECORDS says so). Add "(double-precision screening)".
5. RECORDS.md l. 27: "false from n = 1,565,734,083" suggests that this is the first counterexample. Minimality is unproved: between n = 2.54·10^8 and 1.57·10^9 only 64 windows and proxy-selected bins were scanned. Fix: "false at n = 1,565,734,083 (first found)".

## 5. Side observations

- Jerby's printed viscosities are approximate. μ(g_730119) is 4.4577335 exactly; he prints 4.46023 (2310 l. 1270). μ(g_9807962) is 0.0425003; he prints 0.0750883 (2406 l. 757). Neither his Z_{[t/2]}(t; 1̄) sum nor the bare RS main sum reproduces his numbers (jerbynum.py). The pattern and isolation status agree with his. This is harmless for the refutation: the conjecture is stated for Hardy's Z, and under both of his approximations the first counterexample stays below 4 (μ = 3.96819 and 3.97012).
- Novelty re-check (19 Nov 2026): arXiv monthly listings 2026-10 and 2026-11 (math.NT, math.GM, math-ph, math.CV, math.CA; 4,025 titles), the arXiv API (title "Gram"; abstract "Jerby", "bad Gram", "Gram block") and Zenodo ("Jerby", "bad Gram", Gram-point repulsion/viscosity). The only Gram-point paper newer than the discovery (11 Nov) is Goulden 2611.23332 (16 Nov, "Gram's law for moving targets"); its abstract does not mention Jerby, viscosity or repulsion. No earlier or competing counterexample was found.
- The conjecture's second sentence ("corrupt Gram points must be non-isolated", 2406 l. 781; 2310 con:8.1, l. 753) is qualitative ("dramatically below C"). Our μ ≈ 3.92–3.97 refutes the quantitative bound, not that sentence.

## 6. Dependents (rg -n -i jerby discoveries/ TREE.md COROLLARIES.md RECORDS.md REFUTATIONS.md DISCOVERIES.md FORMALIZATIONS.md NOGO.md)

All stand:
- DISCOVERIES.md l. 12.
- discoveries/jerby-repulsion.md: corrections 1–4.
- RECORDS.md l. 27: correction 5. Its 4.21553 at n = 69,784,845 was rechecked: μ = 4.2155318, isolated bad.
- FORMALIZATIONS.md l. 13.
- TREE.md l. 12 (N8: "false: certified counterexamples … zero-free length-2 Gram blocks") and l. 477 (R441).
- NOGO.md l. 599 (conjfals).
- COROLLARIES.md has no entry.

REFUTATIONS.md l. 37–38 (Jerby's PCC proof; 'corrected Gram's law ⟺ RH') are independent results. They were not re-audited here.

## Files (Riemann/Deep/JerbyAudit/)

- indep.py, logs/indep_cases.log, logs/indep_fidelity.log: independent μ (mpmath) and convention checks; 18 MB peak, 29 s.
- jerbynum.py, logs/jerbynum_jerby.log, logs/jerbynum_ours.log: Jerby's approximations; 130 MB, 31 s.
- zerocheck.py, logs/zerocheck.log: sampled zero counts; 210 MB, 22 s.
- zf_mv.py, logs/zf_mv.log, logs/zf_mv_test.err (negative control): mean-value Arb cover; 14 MB, 232 s.
- FULLPROOF.md.
- The .err files hold the /usr/bin/time -l and memcap lines.

RH relevance: N8, N13
