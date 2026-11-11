AUDIT PASSED + BARRIER (model-based): no flaw in the checked parts of Mori's unconditional e^{−4πμ} bound (inputs, constants, Tables 1 and 3 rebuilt); RH with zero counting and every known multiplicity or pair statistic cannot give lim (−log λ_min)/μ = 4π — in RH-consistent models each zero lost below 4πμ (a multiple zero, or a pair closer than δ* ≈ e^{−3.2μ}) multiplies λ_min by (y/(2+√(4−y²)))⁴, y = γ/2πμ (conjectured lost-sample law from Mori's own dual certificate; numerics μ = 6–20, certified at μ = 5, 10), so a proof needs o(μ) such zeros for ζ, which no known result gives.

# MoriRate — Mori's e^{−4πμ} law for windowed Weil forms: audit, rate, constant (morirate, Orca, 11 Nov 2026)

Sources: R. Mori, Zenodo 24116425 v4 (11 Nov 2026, "sharp-4π") and 24047156 v4 (8 Nov 2026, "paper III": endpoint-smoothed prolate vector,
Thm F/G), both unrefereed. Window [λ⁻¹, λ], a = log λ, μ = λ², s = log μ; our window [−(ln c)/2, (ln c)/2] has c = μ. Under RH,
W(u) = Σ_ρ m_ρ|û(γ_ρ)|² over both signs (Mori Lemma 2.1). Galerkin: explore19's wmat.py (CCM Fourier basis, closed forms, Arb
entries), eigenvalues by mpmath at 4πμ/ln 2 + 190–260 bits; every kept run single-threaded, peak ≤ 83 MB where measured (logs/peak_rss.txt; the unmeasured lostsample μ ≤ 14 and
pairtest runs are smaller instances, ps ≤ 45 MB). One
discarded run broke the 300 MB cap: prolate.py at c = 1000π reached 644 MB (scipy returned all 6000 eigenvectors) and was killed
after 3 min; fixed (3 eigenvectors, c/π ≤ 400 guard; rerun at c = 10π: 55 MB).

## Short answer
1. (i) Audit: inputs ([4]'s prolate defect bound, [5, L3.3], III L3.1–L3.2), Trudgian-derived constants, the conversion 237522697 and the
   final constants, Table 1 (rows G(a), G(b)) and the construction (Table 3 norms to 6–8 digits) all reproduce. The printed Table 3
   W-values are 9–25 % low (Galerkin under-convergence; §8 is not used in proofs). Not re-derived: paper III §§4–6. No FLAW.
2. (ii) BARRIER, model-based. λ_min is a monotone functional of the zero multiset; deleting one zero at height γ₀ costs e^{−4g₂(γ₀/2πμ)} within 0.4 nats,
   g₂(y) = log((2+√(4−y²))/|y|) (|y| < 2), 0 (|y| ≥ 2): measured − predicted = −0.38, −0.26, −0.08, −0.12 (even) at μ = 6, 10, 14, 20.
   A close pair costs this only if δ < δ*; a density of lost samples below 4πμ, or merged bands above it, change the rate. Known
   results allow ≍ μ log μ lost samples below 4πμ, so RH plus counting and pair statistics cannot give the limit 4π (the models are
   not ζ: this obstructs proofs, it does not show non-implication); (SEP_B), (EXT), (CENT) follow from nothing known.
3. (iii) −ln ε_c − 4πc + 5 ln c = −14.9 ± 0.4 for every certified c (5 … 47): λ_min ≈ e^{15±1}μ⁵e^{−4πμ}; Mori's O(log μ) is −5 log μ.

## (i) Audit of the unconditional upper bound
Audited: λ_min ≤ 1.362·10⁹μ²³e^{−4πμ} (μ ≥ 5; [16, Thm B], 24116425 (1.1)) and paper III Thm G: ≤ 1.437·10²³μ^{9/2}(log μ)⁴e^{−4πμ}
(50 ≤ μ ≤ 10⁴), ≤ 1.24·10²³μ⁵(log μ)⁵e^{−4πμ} (μ ≥ 50); vector k̃ = E(h̃)|_{[−a,a]}, h̃ = L(x/λ)(Ĩ₄h₀ − Ĩ₀h₄), h_n = ψ_n(x/λ)/√λ,
ψ_n prolate at c = 2πμ, L = (1 − e^{−c²(1−u²)/4})². In paper III, W(k̃) is a zero sum (radical lemma), so no prime-sum or
archimedean estimates enter; the zero side uses Trudgian's count, Mossinghoff–Trudgian and RH to 3·10¹².

| check | where | result |
|---|---|---|
| Trudgian constants | III v2 note, App. A | published 0.112, 0.278, 2.510: C_N = 7/8 + 2.510 + 0.2/e = 3.385 + 0.2/e ✓; 𝔫(x) = (1/4π + 0.224)x + 0.556 log x + 5.02 + 0.4/e = 2× Trudgian plus the main term of a length-½ window ✓ (read on the PDF page image: pdftotext drops the parentheses) |
| 1 − λ₄(c) ≤ 6621c^{9/2}e^{−2c}, c ≥ 10π | [4]; B ≤ U(c)H² | d₄/bound = 0.139, 0.151, 0.160, 0.178, 0.180 at c/π = 10, 14, 20, 100, 200 (d₄ = 2.60e−18 … 8.29e−531, 77–590 digits); d₄/Fuchs = 0.762 … 0.987 → 1, so the ratio → 1210/6621 = 0.183 ✓ |
| A_n² ≤ 4cd_n/λ_n | [5, L3.3] | ratio 0.248–0.250 (n = 0), 0.231–0.249 (n = 4) ✓ |
| I₀² ≥ πℓλ/c | III L3.2 | β₀²/(πℓ/2c) = 6.3, 7.5, 9.0, 20.2, 28.6 ✓ |
| \|ψ_n(1−x)\| ≤ A_n e^{2c√x} | III L3.1 (first bound) | max ratio 0.994, 0.991, 0.988, 0.940, 0.884 on 36 points of (0, 1] ✓ |
| ‖k‖² ≥ 0.109H² (μ ≥ 50) | [5, Thm D] → III Prop 3.3(c) | ‖k‖²/H² = 0.22269, 0.22144, 0.21949, 0.21937, 0.21931, 0.21927 at μ = 5, 7, 50, 100, 200, 500 (float64; μ = 5, 7 equal the 40-digit values): margin 2.0× ✓ |
| 6621(2π)^{9/2}/0.1089 ≤ 237522697 | III §7 | 237522696.91 ✓; ×6.048·10¹⁴ = 1.4365·10²³ < 1.437·10²³ ✓; ×5.2·10¹⁴ = 1.2351·10²³ < 1.24·10²³ ✓ |
| Table 1, rows G(a), G(b) | III p. 3 | 1.26, 6.22e−5; 30.1, 1.17e−2, 4.69e−11, 1.88e−28 = printed (1.26, 6.2e−5; 30.1, 1.2e−2, 4.7e−11, 1.9e−28) ✓; row G(c) needs P_sm, not recomputed |
| Table 3: vector rebuilt from §1.3, W from the zero side (zeros to 1000 at 40–49 digits, + tail) | III §8 | ‖k‖²/H² = ‖k̃‖²/‖h̃‖² = 0.22268534 (μ = 5), 0.22143867 (μ = 7) = printed 0.222685, 0.221439 ✓. W(k)/B = 1.76, 2.05 and W(k̃)/B̃ = 2.40, 3.07 at μ = 5, 7 (printed 1.6, 1.8 and 2.1, 2.3); the partial sums to height 1000 (1.63, 1.82; 2.34, 2.68) are already lower bounds, so the printed Galerkin values are 9–25 % low, beyond the 'few percent' stated. Harmless (Thm F allows 6·10¹⁴(log μ)⁴·B). Rayleigh quotients 1.51e−17, 2.06e−17 at μ = 5 against λ_min ≈ 9.16e−18 (ChekAudit, N = 1600): within 1.7–2.3× of optimal |
| bound values | | μ = 5: Thm B 8.38e−3; μ = 47: 1.23e−209 (certified ε₄₇ = 1.87e−242); μ = 100: Thm B 2.4e−491, G(a) 1.15e−511, G(b) 4.6e−510; truth ≈ 1.4e−529 (fit), Wiktorowicz core bound 4.58e−530 |

Scope: AUDIT PASSED for the inputs, constants, final arithmetic, Table 1 and the construction. Thm B's proof is in [16] (Zenodo, not in
the two records audited): only its values were checked. The §§4–6 chain (cancellation lemma, high-
frequency bounds, boundary corrector, zero sums above 3·10¹², strip integrals, C_L, C_M, C_T, C_B) was not re-derived; the bound sits
10^15 (μ = 5) to 10^18 (μ = 100) above the truth, so only a structural error could make it false.

## (ii) Rate: where the hypotheses enter, and why RH alone cannot give 4π
Where (24116425). (1) Lemma 2.1: λ_min ≥ s·A((a/π)Γ) for ordinates counted once: multiplicity is dropped. (2) Lemma 2.3: thinning to an
h-separated set loses N − P_h points. (3) Prop. 2.4 (Voss–Hall): A ≥ exp[−c − (D + 3)(2 log(1/q) + 4 log(1/ε) + 15.8)]: each unit of deficit
D costs 2 log(1/q) + …, q = separation. (4) Central range [−e^{B_E}μ, e^{B_E}μ]: D ∋ A_μ = P_{2π/s³}, cost c₂(A_μ + s³) log s in (1.4),
≈ 20 log s per pair (Remark 6.3) ⇒ (CENT). (5) Cμ < |t| ≤ H(μ): (EXT) makes the thinned density ≥ d + 1 on intervals of length ≥ s².
(6) Above: RH + |S(t)| ≤ (1/4 + o(1)) log t/log log t ⇒ clusters ≤ (1/2 + ε) log t/log log t ⇒ no deficit above exp(μ^ϖ), ϖ > 1/2. Thm 1:
(SEP_B) bounds X_μ(h) at h = exp(−μ^{ϖB}); cost s³μ^{ϖB} = o(μ) iff B < 2. Thm 3 (RH alone): deficit ≍ s·e^{μ^ϖ} ⇒ exp(−exp(μ^ϖ)).

Real or artefact? Tested on models that keep positivity (RH): the true zeros with zeros deleted, merged or moved, W_model = W − 2|û(γ_k)|²
(+ 2|û(γ′)|²), a sum of squares.

T1. Lost-sample law (lostsample.py, every zero to y = 3, both sectors). Deleting ρ₀ = ½ + iγ₀ (with ρ̄₀) multiplies λ_min by e^{−4g₂(y₀)+o(1)}:
| μ | N | mean(measured − 4g₂), 0.3 ≤ y ≤ 1.8, even / odd | range, even | first zero: measured / 4g₂ | mean for y ≥ 2.2 |
|---|---|---|---|---|---|
| 6 | 52 | −0.38 / −0.29 | [−0.41, −0.25] | 9.18 / 9.43 (y = 0.375) | 0.04 |
| 10 | 71 | −0.26 / −0.34 | [−0.28, −0.19] | 11.37 / 11.50 (y = 0.225) | 0.04 |
| 14 | 108 | −0.08 / −0.12 | [−0.10, −0.01] | 12.89 / 12.85 (y = 0.161) | 0.03 |
| 20 | 170 | −0.12 / (even only) | [−0.13, −0.07] | 14.29 / 14.28 (y = 0.112) | 0.03 |
N×1.5 at μ = 6 moves the log-ratios by ≤ 0.03. Mechanism (DERIVATION.md): first-order sensitivity of Mori's extremal problem (Thm 4).
The deficit is μ(−log|y|)₊dy, −log λ_min ≈ 2μP with P = 2π, and the dual certificate is the harmonic measure ω₂ of ℂ∖((−∞,−2]∪[2,∞))
at 0, with Green function g₂ (∫g₂(−log|y|)dy = 2π to 15 digits). A lost sample adds unit masses at ±y₀, so P grows by 2g₂(y₀)/μ and
−log λ_min by 4g₂(y₀); zeros above 4πμ cost nothing singly (supp ω₂). Example: at μ = 47 the first zero is worth (0.0479/4.00)⁴ ≈ 2·10⁻⁸
of λ_min. This is the quantitative form of explore19's "visible region γ ≲ 4.4πc" (ExploreAntenna).
Certified instance (certdel.py, Arb; the first zero enclosed by a certified sign change of Z): an explicit even f has
(W − 2|û(γ₁)|²)(f) ≤ 2.2375·10⁻²¹‖f‖² at μ = 5 and ≤ 2.2668·10⁻⁴⁸‖f‖² at μ = 10, while W ≥ ε₅ = 7.447·10⁻¹⁸ and ≥ ε₁₀ = 9.134·10⁻⁴⁴
on even f (certified; ChekAudit): removing the first zero lowers the bottom by ≥ 3.3·10³ and ≥ 4.0·10⁴ (law: e^{8.69} = 5.9·10³,
e^{11.50} = 9.9·10⁴; converged λ_min/X = 4.1·10³ at μ = 5, Galerkin 8.7·10⁴ at μ = 10). Rigour rests on wmat.py's closed forms (checked against quadrature in explore19).

T2. Close pair (pairtest.py, μ = 8, γ₆ moved to γ₅ + δ, y = 0.75): λ(δ) = λ_merged + Kδ² for δ ≤ 10⁻¹² (K = 2.30·10⁻¹³, quadratic to 6
digits over 8 decades), crossover near 10⁻¹⁰, and λ(δ) = 1.00–1.78 λ_orig for 10⁻⁸ ≤ δ ≤ gap: a resolved close pair helps.
λ_merged = λ_deleted = 1.98·10⁻³λ_orig (log-ratio 6.23; 4g₂ = 6.56). Crossover scale (kstar.py, first-order K, y ≈ 0.75):
ln δ* = −15.9, −22.9, −28.9, −35.4, −41.8 at μ = 6, 8, 10, 12, 14, so δ* ≈ e^{−3.24μ+3.6}. Exact dependence:
λ_min(δ) ≈ min{λ_res·e^{−4g₂(y)} + 2F₀′(γ)²δ², λ_res}: a single pair never drives λ_min to 0 (λ_min(Γ) ≥ λ_min(Γ minus a label)), it costs
at most e^{−4g₂(y)}.

T3. Several lost samples (μ = 10): super-additive, as convexity of the extremal value predicts: {1, 2}: 21.45 vs Σ4g₂ = 21.40; {1..5}: 48.12 vs
47.03; {1..10}: 86.12 vs 81.37; {2, 4, …, 20}: 68.92 vs 64.04.

T4. Rate change (exttest.py, even, extra −log λ_min). Every 4th zero in (0, 4πμ] deleted: 26.5, 42.7, 60.7, 79.1, 98.3 at μ = 6 … 14 (slope
≈ 9, rising). Every 8th zero merged into its lower neighbour (double zeros; 25 % of labels multiple, Σm² = 1.25N): 15.8, 22.8, 32.4, 40.2, 50.6
(slope ≈ 4.4). Every 16th merged: 11.4, 18.9, 28.6 at μ = 6, 10, 14 (2, 3, 5 merges; slope ≈ 2.2). The extra rate is ∝ p, as Σ4g₂ ≈ 4πpμ log μ predicts
(slope 4πp(log μ + 1) = 10.4, 5.2, 2.6 at μ = 10 for p = 1/4, 1/8, 1/16).

T5. Exterior bands (exttest.py, even, N×1.7–1.8). Every 2nd zero with y ∈ [2.05, 4] deleted at μ = 8 (25 zeros, singles sum 1.32):
λ_min falls by e^{34.3}. Band y ∈ [2.05, 3.45]: extra 12.5, 24.9, 36.5, 48.2 at μ = 6, 8, 10, 12 (singles sum ≤ 1.3). Merged pairs above 4πμ leave
(1/4π)log(t/2π) < (1/2π)log μ distinct points (sub-Nyquist up to t = 2πμ²), and the collapse is collective. So exterior non-clustering
information is needed by the statement, not only by Mori's method ([19, Prop. 8.1] is about the method).

Verdict (ii): BARRIER, model-based (numerics + first-order potential theory; the lost-sample law is a conjecture). Positivity (RH)
enters only through the zero multiset, and the rate is not a function of RH and counting: (a) ≥ εμ lost samples at y ≤ 2 − η raise the rate by ≥ c(η)ε; (b) a proportion p of double
zeros below 4πμ gives −log λ_min ≳ 4πμ + 4πpμ log μ; (c) merged bands above 4πμ collapse λ_min. Sharpest form: merging a proportion ε/log T
of the zeros at height T changes every asymptotic statistic (N(T), S-bounds, Montgomery's F(α), simple and distinct proportions) by o(1),
yet ≈ 2εμ labels below 4πμ are lost, Σ4g₂ ≈ 4πεμ (the mean of g₂ on (0, 2) is π/2), and the rate becomes ≈ 4π(1 + ε). A proof of the
limit from RH must show that ζ has o(μ) (g₂-weighted) multiple zeros and pairs closer than δ* below 4πμ, and no sub-Nyquist clustering
above. RH-alone heuristic for Thm 3: clusters of the size the CCM S-bound allows ((1/2)log T/log log T) leave (1/π)log log T distinct
points per unit, below Nyquist (1/2π)log μ exactly for T < e^{√μ}; so ϖ > 1/2 is the natural limit of RH + S-bounds.
Necessary vs Mori: the harmful scale is δ* ≈ e^{−πμ}, not 2π/s³ or exp(−(log m)^B) (pairs at Mori's scales cost nothing; they help); his
count o(μ/log log μ) is of the right order (o(μ) is necessary).

Known about close pairs, and whether Mori's hypotheses follow: unconditionally > 67.25 % of zeros are simple and on the line and ≥ 83.62 %
distinct (Lamzouri, arXiv 2609.02882 v2, after a Claude-produced proof); under RH ≥ 19/27 simple (Bui–Heath-Brown, arXiv 1302.5018;
Montgomery 2/3; Cheer–Goldston 0.6727); RH ⇒ multiplicity ≤ (1/2 + o(1))log T/log log T; small-gap theorems only show that close pairs
exist; no lower bound for |γ − γ′| is known; Montgomery's pair correlation conjecture gives o(N(T)) = o(T log T) multiple zeros, not
o(T). All zeros below 3·10¹² are simple (sign changes, Platt–Trudgian), so for μ ≤ 2·10¹¹ the central range has no lost samples: the
barrier is asymptotic. (SEP_B) needs simplicity: open. (CENT) needs o(T/log log T) multiple zeros up to T: open even under RH + PCC.
(EXT) needs uniform local pair bounds up to H(μ): open. None follows from the known results.

## (iii) The constant
Certified ε_c (Chekmachev; c = 5, 10, 16 re-certified in ChekAudit, c = 47 in ExploreFMM): −ln ε_c − 4πc = −23.39, −26.56, −28.78,
−30.63, −31.10, −33.73, −34.02 at c = 5, 10, 16, 23, 25, 45, 47; plus 5 ln c: −15.35, −15.05, −14.92, −14.95, −15.00, −14.70, −14.77
(flat). Galerkin λ_min (WeilSimple U(c), c = 5–20): ln λ_min + 4πc = 5.00 ln c + 15.86 (residual ≤ 0.24). Hence −log λ_min = 4πμ − 5 log μ
− 15.9 ± 0.6 on 5 ≤ μ ≤ 47: Mori's O(log μ) is −5 log μ. Unconditional bounds: −log λ_min − 4πμ ≥ −23 log μ − 21.0 (Thm B; 35–75 nats
loose at c = 5–47), ≥ −5 log μ − 5 log log μ − 53.2 (G(b): same power μ⁵, ≈ 10^16(log μ)⁵ loose); G(a)'s μ^{9/2} is compatible only on
its finite range (it would cross e^{15.9}μ⁵ near μ ≈ 10^40). Under Mori's hypotheses the lower side is −log λ_min ≤ 4πe^{2/s}μ +
c₂(A_μ + s³)log s, a second-order term ≈ +8πμ/s instead of −5 log μ.

## Review
reviewer (subagent, 11 Nov 2026): INCOMPLETE. (1) GAP: the dual certificate gives only a lower bound for the perturbed extremal value
(rigorous: p(σ) ≥ 2π + ∫g₂ dΔ for every admissible σ, so in Mori's multiplier problem lost samples cost at least Σ2g₂/μ), no matching
upper bound, and no proof that the multiplier value fixes the spectral shift to o(1): the law is a conjecture. (2) GAP: model multisets are
not zeta zero sets, so this is an obstruction for proofs using positivity, counting and pair statistics, not non-implication from RH;
many deletions need uniform nonlinear control (T3–T5 are evidence only). Significance: NEW candidate quantitative conjecture. Weakest
point: the transfer from multiplier cost to the spectral minimum. Labels above follow this verdict.

## Novelty
Lost-sample law, crossover δ* ≈ e^{−πμ}, exterior collapse: not in 24116425/24047156 (Thm 4 gives the extremal problem and certificate,
no perturbation law; Remark 6.3 charges ≈ 20 log s per central pair; [19, Prop. 8.1] concerns the method), not in our files (rg on
TREE.md, NOGO.md, Riemann/Deep/*/REPORT.md), no hits in arXiv listings 2026-10/11 (NT, CA, CV, SP, PR, FA, CO; .scratch/lst.py) or
Zenodo ("Weil quadratic form"), nor in arXiv searches (Weil form, zeros, sampling, slit-plane Green function). Candidate for
DISCOVERIES only after a proof of the first-order law (two-sided, uniform spectral comparison); for now a conjecture with evidence.

## Caveats
Galerkin values are Rayleigh–Ritz upper bounds; log-ratios move ≤ 0.03 under N×1.5. Models are not zeta functions; they keep positivity
and show what RH plus counting can and cannot force. The asymptotic law is extrapolated from μ ≤ 20 (offset ≤ 0.13 nats for μ ≥ 14; whether it tends to 0 is open). The first-order law is
derived from Mori's extremal problem by LP sensitivity, not proved.

RH relevance: N1 (rate of Weil positivity; quantifies why window records are finite zero verification; says what a proof of the 4π law needs).

## Files
lostsample.py (T1), certdel.py (certified instance), pairtest.py (T2, T3), kstar.py (δ*), exttest.py (T4, T5), prolate.py (inputs),
testvec.py (Table 3), normk.py ([5, Thm D]), mr_numbers.py (Table 1, constants, ∫g₂), stats_ls.py, table_ls.py, DERIVATION.md; logs/.
Rerun: uv run --no-project --with mpmath --with python-flint python3 lostsample.py 10 3 (≈ 1 min); prolate.py/testvec.py/normk.py need
--with numpy --with scipy.
