# JerbyAudit — the refutation chain, step by step, with source lines

Sources: arXiv e-prints fetched 19 Nov 2026 (2406.06548: Nonlinearity-10-2024-YJ.tex, 916 lines, gz dated 26 Mar 2025; 2310.14415: Paper-RH-2.tex, 1500 lines; 2610.21952: Clouds.tex, 2359 lines).

## S1. Jerby's definitions (what is refuted)

| Item | 2406.06548 | 2310.14415 | 2610.21952 |
|---|---|---|---|
| Z(t) = e^{iθ(t)} ζ(1/2 + it) | l. 184 (abstract l. 164 prints ζ·e^{−iθ}: sign typo) | l. 180, 288 | — |
| θ(t) = arg Γ(1/4 + it/2) − (t/2) log t | l. 188 (typo: log π intended; "Riemann–Siegel θ", refs Edwards, Ivić) | l. 184 (same typo) | — |
| g_n: unique solution of θ(g_n) = πn | l. 164, 204–207 | l. 196–198 | l. 110 |
| Gram's law / good: (−1)^n Z(g_n) > 0; bad: fails | l. 164, 209–211 | l. 202–207 (bad = "fails to meet": (−1)^n Z ≤ 0) | l. 119 |
| first bad Gram point n = 126 (Hutchinson) | l. 214; g_1^bad = g_126, g_2^bad = g_134 (l. 743) | l. 211 | — |
| μ(g_n) = |Z′(g_n)/Z(g_n)| | l. 731–737 | l. 258–262 (signed Z′/Z; conjecture uses |μ|) | l. 316–322 |
| Gram block; isolated bad = middle of a block with N = 2 | l. 769 | l. 745 | l. 616–627 |

Conjectures (verbatim content):
- 2406.06548 l. 313–317 (conj:rep): g_n bad with good consecutive neighbours g_{n±1} ⇒ |Z′(g_n)| > 4|Z(g_n)|. l. 777–781 (Conjecture 8.1, label con:8.1): g_n isolated bad ⇒ same; "In particular, corrupt Gram points must be non-isolated". Abstract l. 169: "If g_n is an isolated bad Gram point then |Z′(g_n)| ≥ 4|Z(g_n)|. We conjecture that this ... holds for any n ∈ ℤ." No height range, no other hypothesis.
- 2310.14415 l. 266–268 ("Repulsion G-B-G conjecture"): g_n bad, g_{n−1}, g_{n+1} good ⇒ |μ(g_n)| > 4. l. 753–755 (con:8.1): corrupt points non-isolated (qualitative).
- 2610.21952 l. 323–326: "Initial experiments ... suggest that isolated bad Gram points satisfy μ(g_n) > 4, to height ... n ∼ 2.5·10^7" (an empirical statement, not labelled a conjecture). Remark l. 616–664: isolated badness does not guarantee a two-zero ("confined") configuration since Rosser's rule fails (Trudgian); no Rosser violation among the first 250,000 Gram points; l. 653–654: "the data suggest more than the lower bound μ(g_n) > 4 for isolated bad Gram points associated with a confined pair of zeros" (+ liminf law l. 656–662).

Reading: conj:rep (l. 313) ≡ Conjecture 8.1 (l. 777) ≡ G-B-G (2310 l. 266): the hypothesis is exactly "g_{n−1}, g_{n+1} good, g_n bad" (block of length 2), the conclusion strict, for every n. The standard Riemann–Siegel θ (log π) is the only reading under which Z is real and g_126 is the first bad Gram point; the log t in l. 188 is a typo. The index convention is the standard one (θ(g_0) = 0, g_0 ≈ 17.8456): it is fixed by l. 214/743 (first bad n = 126, second 134), rechecked numerically in S3.

## S2. Certificate audit (Riemann/Deep/JerbyFals/certify_cx.py, certify_zerofree.py; logs/certify_cx{,2,3}.log, certify_zerofree{,2,3}.log)

C1. θ. `theta(t) = Im acb(1/4, t/2).lgamma() − (t/2) log π` (certify_cx.py l. 11). Arb's acb_lgamma is the analytic log Γ, continuous off (−∞, 0] and real on (0, ∞) (Arb docs: "differs from log(Γ(z)) in general"), so Im lgamma(1/4 + it/2) is the continuous arg with value 0 at t = 0: the Riemann–Siegel θ (log π, not Jerby's typo log t). Independent check (S3): Gram points from my own Stirling series θ_S(t) = (t/2)log(t/2π) − t/2 − π/8 + 1/(48t) + 7/(5760t³) + … agree with the Arb log midpoints to ≤ 3.3·10⁻²¹ (inside the printed radii) at all 15 points; θ_S − mpmath siegeltheta ≤ 2·10⁻³¹ at g_n. A 2π branch slip would move g_k by two Gram spacings (≈ 0.69). No slip. HOLDS.
C2. θ′. `dtheta = Re ψ(1/4 + it/2)/2 − (log π)/2` = d/dt Im log Γ(1/4 + it/2) − (log π)/2. Correct. θ′(g_n) = 9.15989, 9.39339, 9.51494 (S3). HOLDS.
C3. Gram enclosures (gram_ball, l. 13–19): 8 midpoint Newton steps from float(mpmath.grampoint(k)), then the interval Newton operator N = m − (θ(m) − kπ)/θ′(X) on X = [mid ± 10⁻⁴⁵], asserting N ⊂ X and rad N < 10⁻⁴⁵. Standard theorem: N(X) ⊂ X ⇒ exactly one root of θ(t) = kπ in X, and it lies in N (θ′(X) ∌ 0 is implied, else the quotient is unbounded and containment fails). By ThetaDerivPos (θ′ > 0 on [7, ∞): Re ψ(1/4 + iy) is increasing in y > 0 since ∂_y Re ψ(x + iy) = Σ_k 2y(x + k)/((x + k)² + y²)² > 0, and θ′(7) = 0.0536 > 0, S3) this root is THE g_k. HOLDS.
C4. Printed radii. The logs print g with str(30) and Z, Z′ with str(16); Arb's printed ball contains the computed ball (arb_get_str guarantee), so the printed radii (≤ 4.3·10⁻²¹ for g_k: 4.33e-21 at k = 3315629579; ≤ 4.7·10⁻¹⁵ for Z, Z′) are decimal-rounding artefacts; the computed Newton radii are < 10⁻⁴⁵. Page text "Radii are ≤ 4·10⁻²² for the g_k" is FALSE AS WRITTEN (cases 1 and 3 print 4.27·10⁻²² and 4.33·10⁻²¹); harmless (any radius < 10⁻¹⁰ suffices). COSMETIC.
C5. ζ and its t-derivative (ZdZ, l. 20–23): acb_series [1/2 + it, i] = s(ε) = 1/2 + i(t + ε); its zeta() has coefficients z0 = ζ(1/2 + it), z1 = iζ′(1/2 + it) = dζ(1/2 + it)/dt. t is the Gram ball, so the enclosures hold for the true g_k. Z = Re(e^{iθ}z0), Z′ = Re(e^{iθ}(iθ′z0 + z1)) = d/dt Re(e^{iθ}ζ). Note: since e^{iθ}ζ(1/2 + it) = Z(t) is real, Re(e^{iθ}·iθ′ζ) = −θ′ Im(e^{iθ}ζ) = 0, so the θ′ term cannot be mishandled; the logs bound |Im e^{iθ}ζ| ≤ 1.8·10⁻⁶⁵. HOLDS (independently: own RS termwise derivative = mpmath siegelz(derivative=1) = Richardson difference quotient of siegelz to 20 digits, S3).
C6. Classification (l. 30): good = (Z > 0 if k even, Z < 0 if k odd) = (−1)^k Z > 0, bad = the strict opposite, both as Arb certified comparisons; iso = good(n−1) ∧ bad(n) ∧ good(n+1); "certified mu<4" = Arb comparison |Z′|/|Z| < 4 on balls. HOLDS. Log values: case 1 (n odd) Z = 0.27634 (n−1 even: good), 18.39715 (bad), 4.07696 (good); case 2 (n even) −4.26923 (good), −13.98639 (bad), −0.32154 (good); case 3 (n odd) 4.06689, 12.80566, 0.28173 (good, bad, good). All three match the page.
C7. Zero-free cover (certify_zerofree.py): start a0 = arb(nstr(g_{n−1}, 25)) − 10⁻¹⁵, b0 = arb(nstr(g_{n+1}, 25)) + 10⁻¹⁵; nstr rounding ≤ 5·10⁻¹⁷ (case 3: ≤ 5·10⁻¹⁶ for 10 integer digits, still < 10⁻¹⁵), so a0 < g_{n−1} and b0 > g_{n+1} (lesson 277: the cover starts below the left endpoint and ends above the right one). Bisection splits [a, b] at the exact binary m = mid((a+b)/2) into [a, m], [m, b] (no gaps); each piece is accepted only if the ball t = arb((a+b)/2, (b−a)/2) ⊇ [a, b] (arb(mid, rad) adds |rad| to the midpoint ball) gives a Z-enclosure excluding 0 (θ and ζ evaluated on the same ball separately: a superset of Z([a, b]), safe). A piece is rejected by SystemExit if width < 10⁻¹²; the logs end with the success line (31,630 / 40,539 / 49,181 pieces) and a common sign. HOLDS. The "min |Z| lower bound 3.5·10⁻⁵" is a ball artefact; the true minimum on each block is |Z(g_{n±1})| = 0.2763 / 0.3215 / 0.2817 (S4).

## S3. Independent recomputation (no Arb, no JerbyFals code): indep.py, logs/indep_cases.log, logs/indep_fidelity.log

Method: mpmath 1.4.1 at 40 digits. θ_S = own Stirling series (to t⁻⁷); g_k = findroot(θ_S − kπ) and, separately, findroot(siegeltheta − kπ) (agree to ≤ 2.4·10⁻³⁸); Z, Z′ by (a) own Riemann–Siegel sum 2Σ_{k≤a} k^{−1/2}cos(θ − t log k) + (−1)^{N−1}a^{−1/2}(C0 + C1/a + C2/a²), a = √(t/2π), termwise t-derivative; (b) mpmath siegelz, siegelz(derivative=1); (c) Richardson central difference of siegelz (h = 10⁻⁶, 5·10⁻⁷). Peak 18 MB, 29 s.

| n | parity | Z(g_{n−1}), Z(g_n), Z(g_{n+1}) (siegelz) | pattern | Z′(g_n) (a) = (b) = (c) | μ (mpmath) | page enclosure | 4 − μ |
|---|---|---|---|---|---|---|---|
| 1,565,734,083 | odd | 0.276340058956784105, 18.3971477855896357, 4.07695940773904781 | G B G | 73.003850400319735626 | 3.9682156849065032291 | 3.968215684906503 ± 2.30·10⁻¹⁶ (diff 2.29·10⁻¹⁶: inside) | 0.0317843 |
| 2,565,007,238 | even | −4.26922836496585988, −13.9863948395329112, −0.321539448947598733 | G B G | 55.174540721880914329 | 3.9448722386936074881 | 3.944872238693607 ± 4.89·10⁻¹⁶ (diff 4.88·10⁻¹⁶: inside) | 0.0551278 |
| 3,315,629,581 | odd | 4.06689291641759336, 12.805664459952126, 0.28173441376197048 | G B G | −50.235085620125960222 | 3.9228800486869671047 | 3.922880048686967 ± 1.05·10⁻¹⁶ (diff 1.047·10⁻¹⁶: inside) | 0.0771200 |

Own RS vs siegelz: |ΔZ| ≤ 2.2·10⁻¹⁸ at all 15 points (C3 truncation); μ(a) − μ(b) ≤ 5·10⁻¹⁹. 4|Z(g_n)| − |Z′(g_n)| = 0.58474, 0.77104, 0.98757. The page's g_n digits, neighbour values, |Z′(g_n)| = 50.23508562 and 4|Z(g_n)| = 51.22265784 (case 3) all agree. The orchestrator's 200-bit μ (page, Checks) agrees to its 19–20 printed digits.

Fidelity of conventions (indep_fidelity.log, 5.5 s): with the standard θ and θ(g_n) = nπ, the bad Gram points with 0 ≤ n ≤ 140 are exactly {126, 134} (Jerby 2406 l. 214, 743: g_1^bad = g_126, g_2^bad = g_134); g_0 = 17.8455995. n = 6708: G G B G G (Jerby 2310 l. 1265: isolated bad). n = 9,807,962: pattern g_{n−2..n+2} = G B B G G, i.e. the block {9807960, …, 9807963} of Jerby 2406 l. 771. n = 730,119: isolated bad (2310 l. 1270). n = 69,784,845: isolated bad, μ = 4.2155318 (= conjfals/RECORDS 4.21553). So Jerby's index, parity and good/bad conventions are the ones the page uses.

Jerby's own viscosity numbers are not exact (jerbynum.py, logs/jerbynum_*.log): μ(g_730119) = 4.4577335 exactly (he prints 4.46023), μ(g_9807962) = 0.0425003 (he prints 0.0750883); neither his Z_{[t/2]}(t; 1̄) sum (4.4582271, 0.0275339) nor the bare RS main sum (4.4884539, 0.0398511) reproduces them. His numerics thus carry errors ≈ 6·10⁻⁴ relative (and ≈ 7·10⁻³ absolute in Z′ at the corrupt point). The conjecture is stated for Hardy's Z (2406 l. 182–184), which is what the page refutes; and the counterexample is robust under his approximations: at n = 1,565,734,083 the RS main sum gives μ = 3.9701184 and Z_{[t/2]}(t; 1̄) (2.8·10⁸ terms) gives μ = 3.9681944, both < 4 (logs/jerbynum_ours.log, 130 MB peak).

## S4. Zero-free blocks, independent sampled check (zerocheck.py, logs/zerocheck.log; not a proof — the Arb cover C7 is)

Own RS (main sum + C0) on a 10⁻⁴ grid over [g_{n−3}, g_{n+3}] (≈ 20,000 points; 8 spot checks vs siegelz: ≤ 2.1·10⁻⁸). Sign changes per Gram interval [g_j, g_{j+1}), j = n−3, …, n+2: case 1 [1, 3, 0, 0, 1, 0]; case 2 [1, 1, 0, 0, 2, 2]; case 3 [1, 1, 0, 0, 3, 1]. Each block [g_{n−1}, g_{n+1}] has no sign change; |Z| is smallest at a block end (0.27634, 0.32154, 0.28173) and rises inward; the two "missing" zeros sit in the adjacent Gram intervals (3 zeros in [g_{n−2}, g_{n−1}) in case 1; 2 + 2 right of the block in case 2; 3 in [g_{n+1}, g_{n+2}) in case 3): a length-2 Gram block with 0 zeros, i.e. a Rosser-rule violation, as the page says. Since min |Z| on each block is ≈ 0.28–0.32 ≫ any numerical error, the zero-freeness is robust.

## S5. Lean statement fidelity (read-only; logs/lean_axioms.out #print output, LEAN.md; not rebuilt)

`JerbyConj := ∀ n g, IsGram n g → Bad n g → ((∃ g', IsGram (n−1) g') ∧ ∀ g', IsGram (n−1) g' → Good (n−1) g') → (same for n+1) → 4|hardyZ g| < |deriv hardyZ g|`; `IsGram n t := 7 ≤ t ∧ theta t = nπ`; `Good n t := 0 < (−1)^n hardyZ t`; `Bad := ¬Good`; `theta t := ∫₀ᵗ thetaDeriv`, `thetaDeriv t := Re ψ(1/4 + it/2)/2 − log π/2`; `hardyZ t := Re(e^{iθ(t)} ζ(1/2 + it))` (Mathlib riemannZeta; `hardyZ_eq`: imaginary part 0; `gamma_phase`, `theta_unique`: θ is the continuous phase with θ(0) = 0). This is conj:rep (2406 l. 313) and Conjecture 8.1 (l. 777) verbatim, with Jerby's typos resolved the only consistent way (S1). The strict conclusion matches l. 315/779; the counterexample theorems give |Z′| < 4|Z| strictly, so the abstract's "≥ 4" reading (l. 169) falls too, though only the ">" version is stated in Lean. Hypotheses: the Arb records (`GramCert k a b zl zu := ∃ t ∈ [a, b], θ t = kπ ∧ Z t ∈ [zl, zu]`, intervals = printed mid ∓ rad, e.g. certAM: t ∈ [568009639.698978187581685032968666, …032969334], Z ∈ [18.39714778558963574, …64426], Z′ ∈ [73.00385040031973562, …74438]; exactly what C3–C5 certify) and `ThetaDerivPos` (θ′ > 0 on [7, ∞): true, S2 C3). Axioms printed: [propext, Classical.choice, Quot.sound] for every main theorem; no axiom/sorry/native_decide/admit/implemented_by in Lean/*.lean (rg). HOLDS.

## S6. Falsification of the audit target (corners)

| Corner | Test | Result |
|---|---|---|
| index/parity convention | bad Gram points in 0..140 with θ(g_n) = nπ; Jerby's blocks | {126, 134} = Jerby's g_1^bad, g_2^bad; {9807960..9807963} reproduced: same convention |
| 2π branch slip in θ | Stirling θ vs Arb/mpmath Gram points | agree to 10⁻²¹; a slip would shift k by an even number (parity, hence good/bad, unchanged) — harmless anyway |
| π slip / log t typo | literal l. 188 θ with log t | e^{iθ}ζ not real; (−1)^n Z > 0 meaningless; standard θ forced |
| neighbour itself bad | sign margins | |Z(g_{n±1})| ≥ 0.2763 ≫ radii 5·10⁻¹⁸; parity of n ± 1 opposite to n |
| Z′ sign | case 3 has Z′ < 0 | conjecture uses |Z′| (2406) or |μ| (2310): irrelevant |
| θ′ term in Z′ | Re(e^{iθ}iθ′ζ) = −θ′ Im(e^{iθ}ζ) = 0 | cannot be mishandled; Z′ matches FD and siegelz(derivative=1) to 20 digits |
| strict vs non-strict | 4 − μ = 0.0318, 0.0551, 0.0771 | refutes > and ≥ |
| Jerby's numerics instead of Z | his Z_{[t/2]}(t;1̄) and bare RS main sum at case 1 | μ = 3.96819, 3.97012 < 4: robust |
| "isolated" = confined pair? | 2406 l. 769, 2310 l. 745 define isolated by signs at Gram points only | not a hypothesis of conj:rep; 2610.21952's confined-pair version (l. 653–654) is NOT refuted (blocks are zero-free), as the page's Caveat says |
| zero-free block | own mean-value Arb cover (zf_mv.py) + own sampled RS (zerocheck.py) | certified zero-free (414 / 473 / 512 pieces, cover starts 10⁻¹² below g_{n−1}); negative control n = 6708 stops at the Lehmer zero 7005.0628661751 |
| Jerby's quoted μ | μ(g_730119), μ(g_9807962) | exact 4.4577335, 0.0425003 vs printed 4.46023, 0.0750883 (his numerics are approximate; side note) |
