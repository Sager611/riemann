NO IDEA: Erdős–Turán with Landau–Gonek inputs only re-derives the Hlawka–Fujii–Ford–Zaharescu law of {αγ}. Given the primes, its ordinate sums differ from Landau's value only by −2Σ_{δ>0} x^{iγ}(cosh(δ log x) − 1), which Selberg's density bound makes ≪ T log²x/log T. So the law holds for every zero set obeying Selberg's bound and certifies nothing beyond zero density, and it never reaches N8's local quantities. The RH-sensitive variant, the on-line share of Landau's sum rule (1 for ζ under RH; 0.38 at x = 6 and −1.86 at x = 2 for Davenport–Heilbronn at T = 2500), is not fixed by Landau plus Selberg, and no unconditional method for it is known (node N8).

# ExplorePotential: potential theory → N8 (explore35, 12 Nov 2026, 03:09–04:15 UTC)

Derivations: [DERIVATION.md](DERIVATION.md). Code: [code/](code/). Logs: [logs/](logs/).
Notation: zeros ρ = β + iγ, δ = β − 1/2, 0 < γ ≤ T, counted with multiplicity. S(x, T) = Σ x^{iγ} sums over ordinates only; L(x, T) = Σ x^ρ.

## 1. Choosing the bridge (each measured before choosing)

| bridge | measurement | verdict |
|---|---|---|
| N1: CPD of a log kernel | In x-space Weil's archimedean kernel is −κ(\|x\|) = −e^{\|x\|/2}/(2 sinh\|x\|) ~ −1/(2\|x\|). That is a Riesz kernel with s = d and symbol 2θ′(t) ~ log(\|t\|/2π), not a logarithmic kernel. On mass-zero even f, the conditional (CPD) sector, Arch ≥ 0 holds up to c = 3.316 (a = 0.5995; float Galerkin, N = 50 and 90 agree to 4·10⁻⁴; logs/cpd_arch*.log). The margin there is 0.27 at c = 2.5 and 0.088 at c = 3, below the n = 2 atom ln 2/√2 = 0.49. Arch + pole already fails at c = 2.27 (explore32). The bound uses only the gamma factor, so DH satisfies it too. | rejected |
| N1: energy duality / λ_min | Mori's bounds are already potential theory: the harmonic measure ω₂ and the Green function g₂ of ℂ∖((−∞,−2] ∪ [2,∞)), audited in MoriRate and LostLaw. Chekmachev's lattice identity is exact up to node selection, so it has no relaxation whose dual could cap it; only cost limits it (WeilDefl). | already done |
| N10: integer Chebyshev dual | The dual of the resultant method is a measure τ with ∫ log\|Q\| dτ ≥ −ℓ deg Q for every integer Q. Zero distributions of integer polynomials are feasible but cap the method only at the known upper bound 0.4226847; equilibrium measures of capacity-1 sets give only t_Z ≤ 1. No RH content. | rejected |
| N8: Erdős–Turán + Landau–Gonek | §2–§5 | chosen |

## 2. Transferred statement (proofs in DERIVATION §1–5)
- Erdős–Turán: D ≤ N/(K+1) + 3 Σ_{k≤K} \|S(x_k, T)\|/k with x_k = e^{2πkα}.
- Lemma 1 (exact identity, no hypothesis): S(x, T) = x^{−1/2}L(x, T) − 2E_off with E_off = Σ_{δ>0} x^{iγ}(cosh(δ log x) − 1).
  The on-line version is S_on = x^{−1/2}L − Σ_{δ>0} 2x^{iγ} cosh(δ log x).
- Landau–Gonek: L = −(T/2π)Λ(x) + O(x log(2xT) log log 3x + …). Lemma 2 (from Selberg's density bound): E_off ≪ T (log x)²/log T for x ≤ T^{1/8}.
- Corollary: N(T)^{−1} S(x_k, T) = −Λ(x_k)/(√x_k log(T/2π)) + O(1/log² T), with no hypothesis. Under RH the discrepancy is
  D/N ≤ (2πα + 3Σ_k Λ(x_k)/(k√x_k) + o(1))/log T. This is Hlawka's 1975 bound (RH); Fujii has log log T/log T unconditionally, and
  Ford–Zaharescu's 1/log T-order limit measure (unconditional for their class of h) rests on the same two lemmas (TeX math/0405459, Lemmas 1 and 2).

## 3. Tests
Smallest case: α = log 2/(2π), k = 1, x = 2. Predicted S/N = −(T/2πN) log 2/√2.

| zeros | N | S(2)/N measured | predicted | max over x ≤ 40 of \|S − P\| (absolute) |
|---|---|---|---|---|
| ζ, 0 < γ < 2·10⁴ (ExploreHelio) | 22,491 | −0.06926 | −0.06936 | 6.5 |
| ζ, t ≈ 10⁸ (gram_win -zeros, 330 windows, 21 s, 11 MB) | 20,458 | −0.02970 | −0.02956 | 7.7 |
| ζ, t ≈ 9.8·10¹⁰ (explore34's 34 windows, H/2π = 89.8) | 2,107 | −0.01750 | −0.02088 | 17.8 |

Every x ≤ 40 matches, prime powers and composites alike: the absolute error stays O(log T) (≤ 7.7 against \|P\| up to 2,300 at 2·10⁴ and 10⁸; ≤ 17.8 at 10¹¹, where log T = 25 and \|P\| ≤ 66; logs/zeta_ft_*.log). Discrepancies:

| α | D/N at 2·10⁴ | D/N at 10⁸ | Erdős–Turán bound / D | Landau first-order D / measured D |
|---|---|---|---|---|
| log 2/2π | 0.0517 | 0.0221 | 7.9, 8.7 | 0.95, 0.94 |
| log 6/2π (Λ(6^k) = 0) | 0.0034 | 0.0045 | 18, 15 | — |
| √2/10 | 0.0039 | 0.0041 | 15, 15 | — |

The ratio 0.0517/0.0221 = 2.34 matches 1/log T scaling (predicted 2.35). Generic α gives D below √N (87 against 150): the zeros are rigid.

## 4. Where it breaks (no new bound on a node quantity)
- Global only. ET at level K resolves {αγ} at 1/K after summing over all heights. A count in [t, t + h] needs α ≲ 1/h. Then
  either x_k → 1, where S is governed by N(T) and S(T) and gives D_local ≲ √h log T against a count of h log T/2π, or
  x ≥ T, where Gonek's √x log(xT) exceeds the main term TΛ(x)/√x. Gaps, multiplicities and on-line proportions (the N8 line) are out of reach.
- Zero density only (Lemma 1, DERIVATION §5). Given Landau's value of L, the law constrains off-line zeros only through E_off. Per zero, E_off grows like
  (δ log x)²/2 while δ log x is small and like x^δ/2 once δ log x ≳ 1. A moved pair leaves S unchanged and changes L by O_x(1). A family at
  fixed distance with ≍ T zeros would shift S by ≍ T, but Selberg's bound already excludes it for ζ. Every ζ zero set that obeys
  Selberg's bound obeys the law, so the law certifies nothing beyond zero density. This is a statement about this law, not about every
  conceivable ordinate statistic: L itself depends on β.

## 5. Adversarial: Davenport–Heilbronn and moved pairs
DH zeros up to 2500 (dh_ft.py, dh_offsearch.py): 2,406 on-line (sign changes); 109 off-line right members (218 zeros), found from
the count deficit (217.4); δ quartiles 0.016/0.118/0.198/0.259/0.398, none with β > 1. Its Landau term uses c(x), the
coefficients of −f′/f; these are nonzero at composites too (c(6) = 1.936, while Λ(6) = 0).
- With all zeros the sum rule closes: the residual is ≤ 2.1 at every T ≤ 2500 (dh_share.log; T = 420: 9.04 → 1.48 when the 9 known pairs are added).
- Share of the sum rule carried by on-line ordinates (ζ under RH: 1; one moved ζ pair at T = 2·10⁴ changes it by ≤ 2/1560 = 0.13%):

| T | DH on-line / off-line zeros | x = 2 | 3 | 4 | 6 | 7 | 9 | 11 | 13 | 14 |
|---|---|---|---|---|---|---|---|---|---|---|
| 420 | 304 / 18 | −0.86 | −0.34 | 0.73 | 0.66 | −0.12 | 0.77 | 0.85 | 0.06 | 0.65 |
| 1000 | 848 / 56 | −1.30 | −0.80 | 0.70 | 0.57 | −0.29 | 0.71 | 0.90 | 0.07 | 0.61 |
| 2500 | 2406 / 218 | −1.86 | −1.20 | 0.59 | 0.38 | −0.48 | 0.61 | 0.87 | −0.05 | 0.48 |

- The off-line zeros (8% of all zeros) are phase-locked. For γ log 2 mod 2π, 82 of the 109 fall in [3π/4, 5π/4] (uniform: 27); for γ log 3, 83 fall near 0.
  Those are the phases 2^{−iγ} ≈ −1 and 3^{−iγ} ≈ +1, where the leading terms of f = 1 + κ2^{−s} − κ3^{−s} − 4^{−s} + … all oppose the 1.
  The on-line zeros of ζ carry its Landau bias through the same mechanism (Rademacher's predominance).
- The all-ordinate law for DH misses its own Landau value by −2E_off: 0.2–19% of P, up to 36 in absolute terms at x = 39 (dh_alllaw.log; Lemma 1 checked numerically).
- ζ model (zeta_family.py): take the 3,470 ordinates γ < 2·10⁴ with cos(γ log 2) < −0.9 (15.4%), a phase-locked family, as split double
  zeros 1/2 ± δ + iγ. The ordinate law does not change. By Lemma 1, x^{−1/2}L moves by \|2E_off\| = 0.16, 16.1 and 146 at x = 2 for
  δ = 0.01, 0.1 and 0.3, against Landau's error scale of 6.5. So the primes see the family only for δ ≳ 0.064 (δ log T ≳ 0.6). The
  on-line ordinates lose 2Σ_F x^{iγ}, which is −4.30 times the Landau term at every δ.
Verdict: the all-ordinate law separates ζ from DH only through E_off (a zero-density effect) and through c(x) (an Euler-product effect).
The on-line share separates them strongly. For ζ it is 1 under RH, but Landau plus Selberg do not fix it: near-line off-line zeros allowed by
Selberg's bound could carry phase-coherent mass ≫ T. Nothing shows that such mass exists, and no unconditional method for the share is known.
Moved pairs are invisible to both.

Review (reviewer subagent, verdict only, claims inline): Lemma 1, Lemma 2 and the corollary are CORRECT (the constant may depend on x; Λ(x) = 0
off the prime powers). The universal reading "ordinate statistics see only zero density" was INCORRECT as stated and is repaired above to a
statement about the Landau law. "Quadratic" holds only for δ log x small. The non-determination claim is qualified as above.

## 6. Novelty
- §2: known. Hlawka 1975 (ET plus uniform Landau, RH, 1/log T); Fujii, Comment. Math. Univ. St. Pauli 42 (1993) (unconditional, log log T/log T); Ford–Zaharescu, math/0405459, whose
  unconditional measure rests on Gonek–Landau and Selberg's density; Ford–Soundararajan–Zaharescu, 0805.2745 (pair correlation, the
  Selberg class); 2103.01956 (a GL₂ density estimate proved for exactly this use). The TeX was read for the first two.
- §5: an arXiv search (DH off-line zeros, Landau formula for linear combinations, Bombieri–Hejhal) found no on-line-share or
  phase-locking statement for DH. It is a numerical observation, not a theorem; gate not met, no DISCOVERIES entry.

## Suggested NOGO line
- explore35 (Orca, exploration: potential theory — Erdős–Turán discrepancy, Frostman/equilibrium duality, CPD of the log kernel, Fekete/integer Chebyshev → N8, N1, N10) NO IDEA: Erdős–Turán with Landau–Gonek inputs re-derives Hlawka–Fujii–Ford–Zaharescu's {αγ} law (ζ zeros at 2·10⁴ and 10⁸ match Landau's −(H/2π)Λ(x)/√x within 7.7 for all x ≤ 40; ET bound 8–18× the true discrepancy). Since Σ_{all γ} x^{iγ} = x^{−1/2}Σx^ρ − 2Σ_{δ>0} x^{iγ}(cosh(δ log x) − 1) and Selberg makes the last sum ≪ T log²x/log T, the law holds for every zero set obeying Selberg's bound and certifies only zero density. It is blind to moved pairs and to near-line families (a phase-locked 15% family moved to δ = 0.01 shifts Landau's sum by 0.16 against an error of 6.5), and it is global, so no N8 gap or multiplicity bound follows. Davenport–Heilbronn's on-line zeros carry only 0.38–0.61 of its Landau sum at x = 4, 6, 9, and −1.86 at x = 2 (T = 2500). Its 109 off-line pairs close the sum rule (residual ≤ 2.1) and are phase-locked (82/109 with 2^{iγ} ≈ −1). For ζ that share is 1 under RH, but Landau plus Selberg do not fix it. N1: Weil's archimedean kernel is −κ(\|x\|) ~ −1/(2\|x\|), not logarithmic; its mass-zero sector is positive only to c = 3.32 with margin below the n = 2 atom; λ_min's potential theory is Mori's (g₂); N10's dual is a polynomial construction [ExplorePotential](Riemann/Deep/ExplorePotential/REPORT.md)

RH relevance: N8 (ordinate statistics); also N7 (the law is a zero-density statement), N1 (archimedean CPD measured).

## Files
code/: zeta_ft.py (Fourier coefficients of {αγ}, Erdős–Turán vs measured D), parse_gw.py (gram_win -zeros → ordinates; zeros at 10⁸:
`gram_win -zeros 31 16 $(cat logs/nlist_1e8.txt)`, 21 s, 11 MB; zeros file not kept), dh_ft.py (DH on-line zeros, coefficients c(n)),
dh_off.py (T = 420 check with ExploreHelio's 9 pairs), dh_offsearch.py (off-line zeros from the deficit), dh_share.py, dh_alllaw.py,
zeta_family.py (phase-locked ζ family model),
cpd_arch.py with weil_matrix.py (copied unchanged from ../WeilSimple). logs/: every output and /usr/bin/time block. Data:
logs/dh_ft_2500.json (2,406 on-line DH zeros), logs/dh_off_2500.json (109 off-line right members β + iγ below 2500).
Peak RSS ≤ 81 MB for every process (cpd_arch 81, dh_ft 63, dh_offsearch 57, zeta_ft 39, gram_win 11); one process at a time.
