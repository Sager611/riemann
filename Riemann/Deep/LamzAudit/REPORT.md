AUDIT PASSED: Thm Main2 and Prop SmallSetZeros of Lamzouri 2611.14486v1 are correct line by line, and all inputs match their published statements. Sole caveat: the slope-8/7 density estimate is proved only in Baluyot's 2017 PhD thesis (Thm 1.2.1), and we did not re-derive its Kloosterman-sum mean-value core.

# LamzAudit: referee report on Y. Lamzouri, "The distribution of small gaps between zeros of the Riemann zeta function", arXiv 2611.14486v1 (10 Nov 2026)

Claim (Thm Main, l.195-197): D(λ) > 0 for λ ≥ 0.7892, unconditionally, where D(λ) = liminf (1/N(T))·#{0 < γ_n ≤ T : γ_{n+1} − γ_n ≤ 2πλ/log T}.
Line numbers refer to SmallGapsPaper.tex (e-print 2611.14486v1). J_θ(f) := f(0) + 2∫_0^θ αf(α)dα.

## 1. Proof check

| Step (TeX lines) | Check | Verdict |
|---|---|---|
| Lemma ZeroDensity (227-239): Σ_{β≥1/2} x^{β−1/2} ≪ N(T) for x ≤ T^θ, θ < 8/7 | x^{β−1/2} = 1 + log x ∫_{1/2}^β x^{σ−1/2}dσ; with ε = (8/7−θ)/2 the integrand is ≪ N(T)T^{−ε(σ−1/2)}; ∫ ≤ 1/(ε log T). Needs the log T factor and uniformity down to σ = 1/2 in (Eq:Conrey); both hold (§2) | correct |
| Lemma Fujii (242-249) | follows from moments of S(t+h) − S(t) with h log T ≤ 100, since N(t+h) − N(t) = O(1) + S(t+h) − S(t) | correct (input §2) |
| Lemma DoubleSumZeros (251-286) | blocks I_j of length 2π/log T; 1/(1+((γ−γ')log T)²) ≤ 2/(1+(j−k)²); weighted Hölder; ab^q ≤ a^{q+1} + b^{q+1}; I_{j+1} ⊂ (t, t+4π/log T] for t ∈ I_j (also j = J−1); Cauchy-Schwarz + Lemma Fujii with k = q+1 | correct |
| Integral identity (294-305) | BGST Lemma 3 uses only re-indexing ρ' → 1 − ρ̄' (BGST l.236-257), so it holds for every invariant sub-multiset; w(u) = 4/(4−u²) is the analytic continuation of Montgomery's weight. Numerically exact (N1) | correct |
| Kernel bound (307-320): f^(z) ≪ e^{2πθ\|y\|}/(1+x²) | two integrations by parts, supp f ⊂ [−θ, θ] | correct |
| Lemma BGST (322-329) | = BGST Lemma 5 (§2) | correct |
| Prop SmallSetZeros, reduction (340-360) | S_{Z∖E} = S_Z − S_E; cross term ≤ 2√(G(Z)G(E)) by Cauchy-Schwarz (needs f ≥ 0); G(Z) = (J_1 + O(1/√log T))N(T) | correct |
| \|w\| ≤ 4/3 (362-365) | \|4 − u²\| ≥ Re(4 − u²) = 4 + (γ−γ')² − (β−β')² ≥ 3 | correct |
| Grouping (370-375) | fold E onto E_≥ = {β ≥ 1/2}: each off-line class has 2 members, on-line 1; the four terms dominate 2T^{θ\|β−β'\|} + 2T^{θ(β+β'−1)}; each ≤ T^{θ(β+β'−1)} when β, β' ≥ 1/2 | correct (factor 8) |
| Separation (376-380) | T^{θ(β−1/2)}T^{θ(β'−1/2)} ≤ (T^{θ(2β−1)} + T^{θ(2β'−1)})/2, symmetric kernel | correct |
| Hölder (380-389), ε = (1−7θ/4)/2, q = ⌊1/ε⌋+1, p = q/(q−1) | p < 1/(1−ε); θ/(1−ε)² < 4/7 ⟺ (1−u)² > 0, u = 7θ/4; so 2θp² < 8/7 and Lemma ZeroDensity applies with x = T^{2θp²}; G(E) ≪ δ^{1/(pq)}N(T); 1/(pq) = m/(m+1)² ≥ ε/2 needs m = ⌊1/ε⌋ ≥ 3, true because ε < 1/16 for θ > 1/2 (not stated). Constants depend on θ, f only, not on E | correct |
| Thm Main2: E_λ(T) (404-411) | E_λ ∋ every multiple and off-line zero (zero gap to a copy or to 1 − ρ̄, same ordinate), so E_λ is invariant; count of short gaps ≥ #E_λ/2 | correct; invariance also uses that ρ_1 (excluded by n ≥ 2) is simple and on the line, true but unstated |
| Truncation g = ψ_ε f (413-431) | θ + ε = (θ + 4/7)/2 ∈ (1/2, 4/7); g = f on [−θ, θ], g ≥ 0; exponent 1/8 − 7(θ+ε)/32 = 1/16 − 7θ/64 (exact) | correct; the O(1/√log T) term of the Prop is dropped in (Eq:ReductionGoodSet) (harmless) |
| Positivity f ≥ g (430-441) | the integral identity holds for f with support [−1/λ, 1/λ] (each α, finite sums); f ≥ g ≥ 0 | correct |
| Upper bound (443-447) | Z_0 simple, on the line, gaps > 2πλ/log T; diagonal f^(0) = 1, off-diagonal f^ ≤ 0, w > 0 | correct |
| Contradiction (448) | choose δ, then T: N(T) ≥ (J_θ(f) − Cδ^c − C/√log T)N(T), J_θ(f) > 1 | correct |
| Kernel (451-470) | f = λR^(λ·), f^ = R(·/λ); R ≤ 0 for \|x\| > 1, R^ ≥ 0, R^ ∈ C_c²; closed form J = 1.000100147643 (Gauss-Legendre agrees to 12 digits); θ < 1/λ, so the closed form applies on [0, θ] | correct |

No gap found. Minor points for the author:
1. (Eq:ReductionGoodSet) l.425-431: restore the O(1/√log T) term.
2. l.404: say that ρ_1 is simple and on the line (needed for invariance of E_λ and for ρ_1 ∈ Z_0), or put γ_0 := 0.
3. l.389: note ⌊1/ε⌋ ≥ 3, needed for 1/(pq) ≥ ε/2.
4. l.322: BGST normalise by (T/2π)log T; replacing it by N(T) costs 1 + O(1/log T) (harmless).
5. Abstract: log γ_n; body: log T. The body's version is the stronger one.
6. l.246: also cite Tsang 1986, Thm 4, which gives Lemma Fujii directly.
7. Effectivity: at θ = 4/7 − 10^{−5}, q = 228572 and the δ-exponent is 1/16 − 7θ/64 = 1.094·10^{−6} (exact, N8). The proportion is ineffective, and the main-term margin J − 1 = 1.0·10^{−4} beats the lower-order term ≈ −1.64/log T (predicted and observed, N3) only once log T > 1.6·10^4. A remark would help readers.

## 2. Inputs

| Input | As used | As published | Status |
|---|---|---|---|
| Conrey density (l.221-226) | N(σ,T) ≪_ε T^{1−(8/7−ε)(σ−1/2)} log T uniformly for σ ≥ 1/2; N(σ,T) = #{β > σ, 0 < γ ≤ T} | Announced as Thm 2 in Conrey, Bull. AMS 20 (1989) 79-81 (AMS and Project Euclid block scripted access; the same form is quoted by Simonič 1910.08274 and arXiv 2208.11442). Proof: Baluyot, PhD thesis, Rochester 2017 (hdl 1802/32855), Thm 1.2.1: "If 0 < θ < 4/7, then N(σ,T) ≪_θ T^{1+θ(1−2σ)} log T uniformly for σ ≥ 1/2 and large T", same N(σ,T) | Matches (2θ = 8/7 − ε). Proof exists only in the thesis: Conrey's J. reine angew. Math. 399 (1989) 1-26 supplies the Kloosterman-sum technique but not this estimate (Baluyot p.17: "a proof has never been published"); no later journal or arXiv proof found (OpenAlex author list; arXiv search). We checked the thesis's §3.6 deduction (mean value ≤ T + O(T^{1+θ(1−2σ)}) for 1/2 + 1/log y ≤ σ ≤ 1/2 + ε_0 → Littlewood → σ-averaging → dyadic sum; Ingham for σ ≥ 1/2 + ε_0, valid since 3(1−σ)/(2−σ) < 1 − θ(2σ−1) there). Not re-derived: the mollified mean value Thm 1.2.2 (§3.1-3.5, §3.7; Kloosterman sums; the binding errors T^{−1/2}y^{7/8} and T^{−1}y^{7/4} in (3.6.2) both force θ < 4/7) |
| BGST pair correlation (l.294-297, 322-329) | Lemma 3 identity; Lemma 5: Σ f^(i(ρ−ρ')log T/2π)w(ρ−ρ') = (f(0) + 2∫_0^1 αf + O_f(1/√log T))N(T), f ∈ L¹ even real, supp ⊂ [−1, 1], Lipschitz at 0 | Baluyot-Goldston-Suriajaya-Turnage-Butterbaugh, Acta Arith. 214 (2024) 357-376, doi 10.4064/aa230612-20-3; arXiv 2306.04799v1: Lemma 3 (l.236-257) and Lemma 5 (l.349-372) as quoted, with (T/2π)log T for N(T); Thm 1: F(α) = T^{−2α}(log T + O(1)) + α + O(1/√log T), uniformly for 0 ≤ α ≤ 1 | Matches. Refereed, 17 citing works, no erratum found. We checked the off-line modifications in Lemma 4 (factor x^{δ+δ'} ≤ x^{2Θ−1}, Korobov-Vinogradov) and Lemma 5. Only \|α\| < 4/7 is used |
| Fujii moments (l.242-249) | ∫_0^T (N(t+h) − N(t))^{2k} dt ≪_k T, uniformly for 0 ≤ h ≤ 100/log T | Fujii, Trans. AMS 196 (1974), Corollary; h-range corrected in Trans. AMS 267 (1981), per footnote l.246 (originals not accessed: AMS blocks scripted access). Independent: Tsang, Acta Arith. 46 (1986), Thm 4, as quoted in STTB, Trans. AMS 375 (2022) (arXiv 2010.10675 l.151-159): ∫_T^{T+H} \|S(t+h) − S(t)\|^{2k} = c_k H log^k(2 + h log T) + O(H(ck)^k(k^k + log^{k−1/2}(2 + h log T))), for 0 < h < 1 and T^a < H ≤ T | Implied, via a dyadic sum and N(t+h) − N(t) = (ϑ(t+h) − ϑ(t))/π + S(t+h) − S(t). Classical, refereed |

No input is misquoted. Without the thesis, Jutila's published slope 1 − ε (1983; statement as quoted by Baluyot and Simonič) gives θ < 1/2. By our routine deduction the same proof then goes through (take ε with θ/(1−ε)² < c/2), giving D(λ) > 0 for λ > 0.823899 with the Selberg minorant (closed form), or ≈ 0.8209 with a float SDP (d = 31, uncertified) (fallback_theta.py, logs/fallback_theta.log). Both are still far below Zhao's 1 − 2.36·10^{−6}.

## 3. Falsification tests

All computed zeros are simple and on the line. For real data, Lemma ZeroDensity, the grouping step and both Hölder steps are therefore trivial, so we test them on synthetic off-line configurations. The only finite-T content left in real data is the BGST/Montgomery main term.

| # | Test (script, log) | Result |
|---|---|---|
| N1 | identity (l.302-305) for a multiset with 3 off-line pairs (β − 1/2 up to 0.45), a double zero and 2 simple on-line zeros, T = 1000, α ∈ {0.37, −0.81, 1.3} (synth_checks.py) | LHS = RHS to 1.7·10^{−11}, LHS real. Dropping one partner makes LHS complex and the identity fails; Montgomery's weight 4/(4+(γ−γ')²) also fails |
| N2 | closed forms of f^, g0^ (g0 = f·1_{[−θ,θ]}) and J vs Gauss-Legendre (num_pairs.py selftest) | agree to 2·10^{−14}; J = 1.000100147643 |
| N3 | Lemma BGST with g0 on Odlyzko's first 2,001,052 zeros (T = 1.13·10^6), pairs to distance 100, tail 2·10^{−7} (num_pairs.py, logs/num_pairs_zeros{1,6}.log, offset_fit.py) | S_g0(Z)/N = 0.8121, 0.8390, 0.8628, 0.8724, 0.8817, 0.8867 at T = 1.4·10³ … 1.13·10^6. The deficit J − S/N times log T = 1.36 … 1.58 rises toward the predicted 1.641 = (1 + log 2π)(2f(0) − J), which comes from the local density log(t/2π)/2π vs the scaling log T/2π. A GUE + local-density model is off by −0.0136, i.e. ≈ −2.65/log²T, and the gap shrinks. Consistent with J; a constant offset below about 10^{−2} cannot be resolved |
| N4 | Thm Main2 chain at T = 1.13·10^6 with the real E_λ (λ = 0.7892) | δ = #E_λ/N = 0.374, #Z_0/N = 0.626, S_f(Z_0)/N = 0.604 ≤ #Z_0/N, S_g0(Z_0)/N = 0.470 ≤ S_f(Z_0)/N, Cauchy-Schwarz ratio 0.69 ≤ 1. Since δ is large, no contradiction arises, as the theorem predicts. Short-gap fraction 0.194 |
| N5 | Prop SmallSetZeros: removing small E. Scan of E_λ' at T = 7.5·10^4 (num_pairs.py --scan); adversarial E (the δN zeros with the largest local kernel mass) and random E at T = 2.0·10^5 (adversarial_E.py) | removal cost per δ: −1.47 (E_λ', δ ≤ 0.045); adversarial −2.61, −1.84, −1.39 at δ = 10^{−3}, 10^{−2}, 10^{−1}; random −0.95. The cost is O(δ), well inside the claimed O(δ^c) |
| N6 | Lemma DoubleSumZeros, A_q/N for q = 1…6 | 1.085, 1.183, 1.295, 1.428, 1.589, 1.789 at T = 1.13·10^6; bounded, slowly approaching a limit |
| N7 | inequalities l.362-389 on 3 random invariant off-line configurations × 3 values θ ∈ {0.52, 0.55, 0.56}, T = 10^6 | \|w\| ≤ 4/3, grouping (×8), separation, both Hölder steps: all hold |
| N8 | exponent bookkeeping in exact rationals at θ = 4/7 − 10^{−5} | q = 228572, 2θ'p² < 8/7, 1/(pq) ≥ ε/2, 1/8 − 7θ'/32 = 1/16 − 7θ/64 = 1.094·10^{−6} |

Finite-T reading: at every computable height the main-term inequality is invisible (S_g0(Z)/N = 0.887 < 1 < J). The margin of 10^{−4} first beats the −1.64/log T term at T ≈ exp(1.6·10^4). This is expected for an asymptotic theorem with an O(1/√log T) error and refutes nothing, but no finite computation can test Thm Main.

## 4. Literature (11 Nov 2026)

- arXiv: 2611.14486 has only v1 (10 Nov 2026 18:46 UTC), 11 pages, no journal ref, statement "no AI generated text or ideas" (l.214). The math.NT listing for 2026-11 (282 entries, to 2611.14517) has no reply, comment or follow-up. OpenAlex has not indexed it yet; Zenodo and HAL show nothing relevant.
- BGST: published, 17 citing works on OpenAlex (e.g. 2610.23467, 2610.22927, 2609.33043), none reporting a problem. Their follow-up 2501.14545 cites an unpublished 2024 Conrey note ("Zeros near the critical line": > 99.9% of zeros within 8/log T of the line), consistent with the density input.
- Lamzouri's earlier BGST-based paper 2609.02882v2 (Sep 2026) acknowledges that Conrey, Goldston, Gonek, Granville, Milinovich and Tenenbaum read a draft. It uses neither Fujii nor Conrey's density, so Prop SmallSetZeros (removing a small set of zeros via Fujii plus density) first appears in 2611.14486. No independent check of it exists yet.
- Baluyot's thesis is cited for the proof (l.221 footnote). Simonič 1910.08274 calls Conrey's 8/7 "announced"; Baluyot (thesis p.17) says no proof had been published before.

## 5. Referee verdict (plain words)

The argument is sound. If almost all zeros were simple, on the line and well spaced, the unconditional pair-correlation formula of Baluyot-Goldston-Suriajaya-Turnage-Butterbaugh, applied with a Selberg-type kernel, would yield a contradiction. The new step is to show that a tiny set of bad zeros (off the line, multiple, or closely spaced) cannot change the pair-correlation sum much. Lamzouri proves this with Fujii's moment bounds and Conrey's zero-density estimate. This step is correct, and the exponents work exactly up to θ < 4/7, which is half of Conrey's density slope. Every quoted input matches its source. The one soft spot is that Conrey's slope-8/7 density estimate has never appeared in a journal; its proof is in Baluyot's 2017 thesis. If that proof failed, the same method with Jutila's published estimate would still give about 0.82 instead of 0.7892. The result is purely asymptotic and gives no explicit proportion. Recommendation: accept after minor revision (points 1-7 of §1).

For us: Riemann/Deep/SmallGap (λ ≥ 0.7851) and TREE N8 may keep resting on Thm Main2, with the stated caveat that the density input is proved only in a thesis.

Review: an independent `reviewer` read the TeX lines and returned PASS on (R1) the Hölder exponents, (R2) grouping and separation, (R3) the Thm Main2 logic and exponent identity, (R4) Lemma BGST vs BGST Lemma 5, and (R5) the Conrey density vs Baluyot Thm 1.2.1: "no mathematical error or misquotation".

## Files

Riemann/Deep/LamzAudit/: REPORT.md; num_pairs.py (pair sums on Odlyzko zeros, kernels, model); synth_checks.py (N1, N7, N8); adversarial_E.py (N5); fallback_theta.py (Jutila fallback); offset_fit.py; logs/: num_pairs_zeros1.log, num_pairs_zeros6.log, synth_checks.log, adversarial_E.log, fallback_theta.log, offset_fit.log.
Rerun: `uv run --no-project --with numpy python3 num_pairs.py <zeros.gz> 300000,1000000,2001052` (zeros6.gz from Odlyzko's tables; 3 min, peak RSS 110 MB); `uv run --no-project --with numpy python3 synth_checks.py` (5 s).
