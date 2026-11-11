AUDIT PASSED: Baluyot's proof of the mollified mean value Thm 1.2.2 (θ < 4/7) and its deduction of Thm 1.2.1, Conrey's slope-8/7 density estimate N(σ,T) ≪ T^{1+θ(1−2σ)} log T, is correct line by line. Every quoted input matches its source. The one misprint is harmless.

# Baluyot: referee report on the only proof of Conrey's slope-8/7 zero-density estimate

Source: S. A. C. Baluyot, *On the zeros of Riemann's zeta-function*, PhD thesis, University of Rochester 2017, hdl 1802/32855, chapter 3 (printed pp. 93–144). We read all 52 pages from page images. Equation numbers are the thesis's.

## 1. Map

Full map: [MAP.md](MAP.md). In short, the mollifier M(s) = Σ_{n≤y^{1+ε}} μ(n)n^{−s}W(y/n)/W(y), y = T^θ, is a smoothed Selberg mollifier with W(ξ) = 1 + Aξ^{1−2σ} + O(ξ^{−3/4}), A = Γ(2−2σ)ζ(2−2σ)/ζ(2σ). After Gaussian smoothing, the mean square is g(u) = 𝔐_1 + 𝔐_2 + 𝔇_1 + 𝔇_2 + 𝔈. The diagonal 𝔐_1 = 1 + O(y^{1−2σ}) is obtained by residues. 𝔐_2 ≪ y^{1−2σ}, and 𝔇_j are exponentially small. 𝔈 is bounded by Vaughan's identity plus two Kloosterman inputs: Deshouillers–Iwaniec's bilinear bound gives T^{−1/2}y^{7/8} (type II), and Weil's bound gives T^{−1}y^{7/4} (type I). The second is the square of the first, so both are small iff θ < 4/7. Littlewood's lemma, Selberg's σ-averaging, a dyadic sum and Ingham's bound for σ ≥ ½ + ε_0 then give Thm 1.2.1.

## 2. Numerical tests (falsification first)

The theorem's σ-range [½ + 1/log y, ½ + ε_0] is empty at every computable height, because (3.6.5) forces ε_0 < (1 − 7θ/4)/(4θ). At θ = 0.55 that means ε_0 < 0.017, so the range needs log y > 59. We therefore test the formula itself on σ ∈ (½, 1).

| # | Test (script → log) | Result |
|---|---|---|
| N1 | W(ξ;σ) three ways: Mellin series Σ_b b^{1−2σ}/(e^{b/ξ}−1), residue expansion, contour quadrature (selftest_mollifier.py) | agree to 12 digits (σ = 0.52, 0.6, 0.75; ξ = 1 … 10) |
| N2 | 𝔐_1 by divisor sums vs the direct double sum Σ_{h,k} a_h a_k/[h,k]^{2σ} (selftest_mollifier.py) | agree to 15 digits |
| N3 | 𝔐_1 exact for y = 10²…10⁵, σ = ½ + c/log y (c = 1, 2, 4), σ ∈ {0.6, 0.75} (main_scan.py → main_scan.log) | (𝔐_1 − 1)/y^{1−2σ} ∈ [0.54, 0.88] for σ ≤ ¾, with no 1/(2σ−1) loss. At c = 1, 𝔐_1 = 1.0822, 1.1000, 1.1110, 1.1184, rising to 1/(1−e^{−2}) = 1.1565. The thesis's 𝔐_11 = (1 + Ay^{1−2σ})/W(y)² differs from 𝔐_1 by −0.44, −0.31, −0.23, −0.19 × y^{1−2σ}, i.e. 𝔐_12 is bounded as claimed |
| N4 | ζ(σ+it) by vectorised Riemann–Siegel vs Arb (check_rs.py) | max error 6·10^{−4}, 9·10^{−5}, 1.4·10^{−5} at T = 10³, 10⁴, 10⁵ |
| N5 | full I = (1/T)∫_T^{2T} \|ζM\|² for T = 10³, 10⁴, 10⁵, θ ∈ {0.2, 0.4, 0.55}, ε = ½, σ ∈ {0.55, 0.6, 0.7, 0.8, ½ + 1/log y, ½ + 2/log y}, 63 cases, trapezoid step ⅛ (meanvalue.py → meanvalue_T1e{3,4,5}.log) | \|I − 𝔐_1 − 𝔐_2\| ≤ 2.9·10^{−3}, 1.2·10^{−4}, 9.5·10^{−5}, while the thesis's bound for this term is 1.65, 1.55, 1.46 at θ = 0.55. (I − 1)/y^{1−2σ} ∈ [0.37, 1.69]. Example at T = 10⁵, θ = 0.55, σ = ½ + 1/log y: I = 1.092683, 𝔐_1 + 𝔐_2 = 1.092674. Step ¼ changes I by ≤ 1.5·10^{−6} (meanvalue_T1e4_h025.log) |
| N6 | beyond the theorem: θ = 0.7, 0.9, 1.0 at T = 10⁴ (meanvalue_T1e4_beyond.log) | \|I − 𝔐_1 − 𝔐_2\| ≤ 9·10^{−4}. The 4/7 limit is invisible at this height, as expected for a technical barrier |
| N7 | identities the proof quotes (check_lemmas.py → check_lemmas.log, check_lemmas_34.log) | Estermann's functional equation with shifts at 3 points: 9·10^{−25}. Lemma 3.1.3 (= Conrey's Lemma 5) vs the direct series: 6·10^{−16}. Lemma 3.1.7 (= Conrey's Lemma 3): 1·10^{−16}. Lemma 3.1.2, \|D(0,2σ−1,0,H/K)\| ≤ 2K² for K ≤ 12: max ratio 0.37 |
| N8 | exponents in exact arithmetic (exponents.py → exponents.log, sympy) | DI expansion: 7 monomials, (1,2,3,4) not (1,2,3,3). The maximum is T^{−2}y^{7/2} iff y ≤ T^{8/13}. Lemma 3.1.4 contributes T^{−1+7δ/2}. 𝔈_11, 𝔈_32 → 0 iff θ < 4/7, and 𝔈_2 iff θ < 4/5. The type I maxima are y^{7/4} and y^{3/2}. §3.6: (3.6.5) ⇒ (3.6.6), (3.6.7). Ingham comparison: concave, slope 2θ − 4/3 at σ = ½ |

The main term passes every test. The off-diagonal is 10^{−5}–10^{−3}, far below its proven bound, so no finite computation can probe the Kloosterman estimates themselves.

## 3. Referee notes

| Step | Check | Verdict |
|---|---|---|
| Lemma 3.1.1 (W) | residues at s = 0 and 1−2σ, shift to Re s = −¾; (3.1.1) and y^{1−2σ} ≤ e^{−2} give \|Ay^{1−2σ}\| < ½; (iii) via Γ(s+3)/((s+1)(s+2)) with r = 1/ξ − 3 | correct |
| Lemma 3.1.2 | Hurwitz grouping (3.1.7); \|ζ(2σ−1, a/K)\| ≤ 2 + 2K^{2σ−1}, \|ζ(0, b/K)\| ≤ ½ | correct (N7) |
| Lemma 3.1.3 | Mellin of e^{−z}; residues at 0, 1−α, 1−β (Σ_a e(abH/K) = K·1_{K\|b}); Estermann's equation reduces to the classical one at α = β = 0 | correct (N7) |
| Lemma 3.1.4 | factor bounds: \|z^{iu}\| = e^{−λu}, Gaussian, \|Γ(s−2σ+1)\| ≪ (1+\|t\|)^{c−2σ+1/2}e^{−π\|t\|/2}, \|(−iϱ(z−1))^{s−1}\| = a e^{−bt}e^{πt/2}, ∫(1+t)^κ e^{−ϑt} ≪ ϑ^{−κ−1}; min(b, π−b) ≫ λΔ/log Δ near x = 1 | correct; the true size λ^{−5/2}Δ^{−7/2} is smaller by logs |
| Lemma 3.1.5 | e(nb̄/v) = e(nq b̄'/(vq)) with b̄' the inverse mod vq; (nq, vq) = q(n, v) | correct |
| Lemma 3.1.6 | the first display is DI's Lemma 1 (Acta Arith. 43 (1984), p. 310, read from the scan) with ϱ = 1 and (c, d, u, v) → (v, b, ℓ, a); DI's conditions (c, ϱd) = 1, (v, c) = 1 become (b, v) = 1, (a, v) = 1 | correct. Misprint: E lists (1,2,3,3) where the expansion gives (1,2,3,4). (UNV)²(N+A)·A^{−1}VU(N+A)(V+A²) contains U²N²V²·A·A^{−1}VU·A·V = AN²U³V⁴ (exponents.py). Harmless: (1,2,3,4) contributes T^{−2}y^{7/2}, exactly the maximum |
| Lemma 3.1.7 | χ(1−s) = (2π)^{−s}Γ(s)(e^{iπs/2} + e^{−iπs/2}); Mellin–Parseval with f(z) = (Δ/2√π)e^{−Δ²(log z)²/4}z^{−s_0}, g(z) = e(∓z) | correct (N7) |
| Lemmas 3.1.8–3.1.9 | Gaussian integral; 𝒥 = Γ(2σ−1)[(σ−½−iu)^{1−2σ} − (σ−½+iu)^{1−2σ}] + O((2σ−1)Γ(2σ−1)u^{1−2σ} + T^{−δ}); integrations by parts with H_1, H_2 | correct |
| §3.2, Lemma 3.2.1 | Littlewood's lemma on [σ, b] × [T, 2T]; ζM = 1 + O(2^{−σ}); arg bound by Jensen's theorem with M ≪ y^{1+ε}; Jensen's inequality | correct |
| §3.3, Lemma 3.3.1 | shift past s = 1 only; functional equation; L7 with β = 2σ−1, ω = mnh/k; rays ±λ; L3 with w = H(z−1)/K. We re-derived the prefactors: 𝔐_1 gets q^{−2σ}(HK)^{−2σ} = (h,k)^{2σ}(hk)^{−2σ}, 𝔐_2 gets (h,k)^{2−2σ}/(hk), and 𝓘 = 1 is the residue at z = 1 of the sector contour | correct. Immaterial slips: (3.3.15) drops the minus sign of L3's second D, and the residue in (3.3.6) drops a 2πi. Both terms are only bounded in absolute value |
| §3.4, Vaughan | 1/ζ = (1/ζ)(1−ζM_0)² + 2M_0 − ζM_0²; the coefficients of 1 − ζM_0 vanish on 1 < n ≤ y^{1/4} | correct |
| Lemma 3.4.1 (type II) | WLOG α_0 ≤ β_0 by symmetry; W(y/abq) split by its Mellin integral (Re z = 1/log y, (y/q)^z ≪ 1, ζ(1+z) ≪ log y); C(a, ℓ) ≪ 1; then L6. Bookkeeping recomputed: (AB)^{u+n−4}V^{v+n−4} ≤ y^{2n+u+v−8} (all exponents ≥ 0); with y^{1/4} ≤ A ≪ y^{1/2} the per-term bounds are T^{−2}y^{7/2} (three terms), T^{−4}y^{27/4}, T^{−4}y^{13/2}, T^{−4}y^6 and T^{−3}y^{9/2}. Their maximum is T^{−2}y^{7/2} iff y ≤ T^{8/13}, decided by (−1,4,3,4): T^{−2}y^{13/4} ≤ 1. A direct check without the TN trick gives the same T^{1/2}y^{7/8} at N ≈ ABV/T | correct |
| 𝔈_11, 𝔈_12 | sums over N converge since η = 2ε; 𝔈_12's shift to c' = 2σ+η−1 meets only the pole s = 1 of 1/sin πs, because each L there is a finite Dirichlet polynomial; L4 gives T^{−1+7δ/2} and T^{1−2σ+5δ/2} | correct |
| S_2, Lemmas 3.4.2–3.4.3 | trivial bound (y^{1/4}V)^{1+η}N^{−η}; type I: partial summation + L5 (with (ℓā, K) for (ℓā, K)^{1/2}, a weaker valid bound); Σ_{K,ℓ}(ℓ, K) ≪ (NV)^{1+ε}; qAV^{3/2} ≤ y^{7/4}, ABV^{1/2} ≤ y^{3/2} | correct |
| §3.5 | residues of the double Mellin integral; F_1(0,0) = F_1(1−2σ,0) = 1 (local factors (1−p^{−2σ})²/(1−p^{−2σ})² and 1); the 9/16 conditions hold on the shifted lines; no poles in r < Re s_1 ≤ ⅛; uniform in σ because ζ(2−2σ)/ζ(2σ) and 1/W(y) are bounded | correct (N3) |
| §3.6 | 𝔇_j ≪ e^{−½log²T}; the parameter choice in exact arithmetic: (3.6.6) is (3.6.5)/2 and (3.6.7) = (3.6.5) − 7δ/2 − 2ε_0θ; 20 000 random θ < 4/7 give no violation; Δ log T ≪ T^{1−δ/2} ≤ T^{1+θ(1−2σ)}. (3.6.12) is extended to ½ ≤ σ ≤ ½ + 1/log y by N(T) ≪ T log T (p. 142), so Selberg's averaging over [σ − 1/log T, σ] stays in range (a reviewer's GAP, closed by that paragraph). Ingham: 3(1−σ)/(2−σ) − 1 − θ(1−2σ) vanishes at ½, has slope 2θ − 4/3 < 0 there, and is concave | correct |
| §3.7 | lower bound with 0 ≤ v ≤ 1 | correct (not needed for Thm 1.2.1) |

Uniformity: the implied constants depend on θ, δ, ε and ε_0, never on σ. Of the inputs, DI's Lemma 1 (refereed, 1984, from DI Invent. Math. 1982 Thm 12) and Weil's bound are used as published, not re-proved; we checked the statements against the sources. Conrey's Lemmas 3 and 5 and Estermann's equation we re-derived and checked numerically.

Errata for the author: (a) Lemma 3.1.6, E: (1,2,3,3) → (1,2,3,4); (b) (3.3.15): the sign of the second D term; (c) (3.3.6): a factor 2πi in the residue; (d) Lemma 3.4.3: (ℓā, K) can be (ℓā, K)^{1/2}. None affects any bound.

## 4. Literature (11 Nov 2026)

- Conrey, Bull. AMS 20 (1989), Thm 2, is the announcement. OpenAlex W2116122835 lists 20 citing works, among them Simonič (JMAA 2020, 1910.08274: "announced"), "An explicit density estimate for Dirichlet L-series" (Math. Comp., doi 10.1090/mcom/2991), "An explicit log-free zero density estimate" (JNT 2024) and "Toward optimal exponent pairs" (Math. Comp. 2024). None gives a proof or reports a doubt.
- The thesis is not indexed in OpenAlex. Baluyot's 14 arXiv records (2017–2026) contain no publication of chapter 3; chapter 4 appeared as 2508.10857. Lamzouri 2611.14486 cites the thesis for the proof. LamzAudit found no other proof.
- Published results by the same technique back the 4/7 limit. Conrey, Crelle 399 (1989), Thm 2, is the mollified Levinson mean value for θ < 4/7. Deshouillers–Iwaniec, Acta Arith. 43 (1984), Thm 4, give ∫_0^T |ζM(½+it)|² ≪ T^ε(T + T^{1/2}M^{7/8} + M^{5/3}) for a_m = μ(m), hence M ≤ T^{4/7}. That is the same T^{1/2}y^{7/8} term that binds here.

## 5. Dependency map

rg '8/7|Conrey|Baluyot|0.7851|0.7892|smallgap' over RECORDS, TREE, COROLLARIES, DISCOVERIES, NOGO, discoveries/ and Riemann/Deep/:

| Entry | Rests on Thm 1.2.1? | After this audit | With Jutila's slope only (θ < ½) |
|---|---|---|---|
| RECORDS N8 small gaps λ ≥ 0.7851 ([SmallGap](../SmallGap/REPORT.md)); TREE N8 "Gaps" line | yes, via Lamzouri's Lemma ZeroDensity (θ < 4/7) | stands; the caveat "proved only in a thesis" can read "proved in Baluyot's thesis, audited (Baluyot)" | λ > 0.823899 (Selberg minorant, closed form; [LamzAudit](../LamzAudit/REPORT.md) fallback_theta.py) |
| TREE R421, LamzAudit (Lamzouri's 0.7892) | yes | stands; LamzAudit's sole caveat is resolved | same |
| Harvest6 l. 14 | context only (θ as the lever) | – | – |
| IsoZero, discoveries/isolated-zeros, linear-combinations-* | no: they use Conrey's published Crelle Thm 2 (via Bui–Conrey–Young 1002.4127) | unaffected | – |
| ZetaZeros, Harvest2, R018, R182 | no: "Baluyot" there is BGST (Acta Arith. 2024) | unaffected | – |
| NBLp, discoveries/nyman-beurling-lp, COROLLARIES l. 10, HPAudit | no: their 8/7 = 1/(7/8) comes from the zero-free half-plane | unaffected | – |

DISCOVERIES.md has no entry that rests on it.

## 6. Verdict (plain words)

The small-gap records need Conrey's 1989 estimate for zeros near the critical line, and its only written proof is Baluyot's 2017 thesis. We read that proof in full. The proof is sound. Its core is a mean-square estimate for ζ times a short Dirichlet polynomial. The main term is computed exactly and matches computation: at heights up to 10⁵ the true mean square agrees with the predicted main terms to within 10⁻⁴. The hard part is the error term, which rests on two published Kloosterman-sum bounds (Deshouillers–Iwaniec and Weil). They are quoted correctly, and the exponent arithmetic is right. It gives exactly the limit 4/7, the same limit as in the published Deshouillers–Iwaniec and Conrey work. One misprint in a list of exponents changes nothing. So Lamzouri's 0.7892 and our 0.7851 keep their footing. The remaining caveat is only that the proof never appeared in a journal.

Review: two verdict-only `reviewer`s on pasted evidence. Reviewer 1 returned PASS on the DI transcription, the misprint, the maximum T^{−2}y^{7/2} for y ≤ T^{8/13}, and the condition 4/7 < 8/13; its "GAP" asked only to say "the displayed upper bounds are o(1) iff θ < 4/7", adopted above. Reviewer 2 returned PASS on the residues and F_1 values, on uniformity, and on the numerics. Its GAP on the σ-averaging range is closed by the thesis's own extension of (3.6.12) to ½ ≤ σ ≤ ½ + 1/log y (p. 142).

RH relevance: N8 (load-bearing density input of the small-gap records)

## Files

Riemann/Deep/Baluyot/: REPORT.md, MAP.md; code/: zeta_rs.py (vectorised Riemann–Siegel), check_rs.py, mollifier.py (W, coefficients, 𝔐_1, 𝔐_2), selftest_mollifier.py, main_scan.py, meanvalue.py, check_lemmas.py, exponents.py, bench_zeta.py; logs/: main_scan.log/.json, meanvalue_T1e3.log, meanvalue_T1e4.log, meanvalue_T1e5.log, meanvalue_T*_h*.json, meanvalue_T1e4_h025.log, meanvalue_T1e4_beyond.log (run stopped at θ = 1 for memory), check_lemmas.log (sections 1–3; section 4 crashed on a gcd import, rerun in check_lemmas_34.log with a wider quadrature range for 3b), exponents.log.
Rerun: `cd code && uv run --no-project --with mpmath --with numpy --with scipy python3 -W ignore meanvalue.py 100000 0.2,0.4,0.55 0.125 0.5` (6 min, 120 MB); `uv run --no-project --with sympy python3 exponents.py` (4 s); `uv run --no-project --with mpmath python3 check_lemmas.py` (10 min).
