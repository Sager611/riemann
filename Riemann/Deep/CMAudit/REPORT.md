AUDIT PASSED: CMQuartic Lemmas 1–4, Theorem A (for 0.6 ≤ σ ≤ 1), Proposition C, Theorem B, CM2 Theorems 1–3 and CM3 Lemma P_F and Theorems 1–3 are correct as stated; the defects are cosmetic (the discovery file drops Theorem A's σ-range), and no counterexample appears in 22 test fields with d_K up to 10⁸ or in the family sieve over all 1429 CM fields over Q(√5) with 500 < N𝔡_{K/F} ≤ 16000.

# CMAudit: referee of the CM least-prime chain (CMQuartic, CM2, CM3), 2026-11-12

Brief: `.scratch/orca/cmaudit-brief.md`. Long derivations: [DERIVATION.md](DERIVATION.md) (D1–D10). Scripts and outputs: [checks/](checks/). Line numbers refer to [CMQuartic NOTE](../CMQuartic/NOTE.md) ([CMQ]), [CM2 NOTE](../CM2/NOTE.md), [CM3 NOTE](../CM3/NOTE.md). Refereed Q-analogues used as black boxes: [R2] Prop 1, [R3] Lemmas 3–9, [R5] Lemma P ([RC5Audit](../RC5Audit/REPORT.md), [REPORT2](../RC5Audit/REPORT2.md)), [C] Lemmas 2–5 ([RCAudit](../RCAudit/REPORT.md)).

## Verdict table

| Step | NOTE lines | Inputs | Verdict |
|---|---|---|---|
| [CMQ] Facts (a)–(c): Fact (a) on 𝒵_K, gamma factor Γ_C(s)ⁿ, h_K ≤ C(n)√d_K log^{2n−1} | CMQ 35–41 | class number formula, Louboutin 2000, Friedman 1989 | correct, uniform in K (D3, D4) |
| (F4), (F5), (I6) | CMQ 43–47 | IK 5.7, 5.3–5.4, 5.11 | standard, used within hypotheses |
| (HT) over F | CMQ 48–54 | HT TeX cor:MVT l.841–850, bound_int l.949–952 | matches: V = Q^{m}T^{m²[F:Q]/2+1}\|𝒮\| at θ_m = 0 (GRC holds for π_χ), any base field (D5) |
| (TZ) | CMQ 55–58 | TZ TeX Thm 3.2, l.1137–1146 | matches; used once, for the fixed ζ_{H_F⁺} (D7) |
| Lemma 1: Δ = ββ′ᶜ − βᶜβ′ ∈ N_{K/F}(𝔠)𝔇, \|NΔ\|^{1/2} ≥ N𝔠√d_K/d_F | CMQ 64–75 | local monogenicity, f′(y) = y − yᶜ | correct (D1); 3300/3300 pairs in 550 random ideals 𝔠, 455 non-principal |
| Lemma 2: O_F-primitive ideals per class | CMQ 79–99 | Lemma 1, Dirichlet units | correct; strengthened form (1 + 4π(N/√D)^{1/n})ⁿ, no units or regulator (CMQ 299–308), correct (D2); max_Cπ_C/bound ≤ 0.042 |
| Remark 2.2: Z-primitive weights fail by d_K^{1/4−o(1)} | CMQ 105–110 | — | confirmed: unweighted and Z-primitive constants exceed the weighted one by up to 3.8× and 2.5× (d_K > 10⁷, N = √d_K/8), growing with d_K |
| Remark 2.3: Nβ ≥ d_K/(4ⁿd_F²) for β ∈ O_K∖F | CMQ 111–112 | Lemma 1 | correct; equals Miller–Peckner's key lemma (TeX l.618–661, written there only for O_K = O_{K⁺}[√β]) |
| Lemma 3: content-weighted large sieve, constant C_Fε⁻¹log^{2n−1}(3d_K)(N + √d_K), 0 < ε ≤ ½ | CMQ 118–128 | Lemmas 1–2, Fact (c) | correct (D3); exact sharp constant = top eigenvalue of the dual Gram = h_K·max_CΣ_{𝔞∈C}N𝔤(𝔞)^{−1−2ε}, measured ≤ 0.88κ_K(N + √d_K) |
| Lemma 3′: content-weighted divisor sums | CMQ 130–135 | d_k(𝔤O_K) ≤ d_k(𝔤)² | correct |
| Theorem A(a): Σ_χN(σ;t,χ) ≪ (d_KTⁿ)^{2(1−σ)/(3−2σ)+ε} | CMQ 139–168 | [R2] Prop 1 (refereed), Lemma 3, (2.1′) | correct for 0.6 ≤ σ ≤ 1; e(σ), the Tⁿ and the dual term ξ(𝔞ᶜ) re-derived (D4) |
| Theorem A(b) | CMQ 141 | (F4), Fact (c) | correct |
| Proposition C: X^{10(1−σ)/(3−2σ)} over the δ₀-rough family | CMQ 178–189 | HT with m = 2, (2.1′) | correct: C(π_χ) = 12ⁿd_F²D, multiplicity ≤ 6, exponent and T^{(3n+2)/5} recomputed (D5) |
| Lemma 4: genus characters = ψ∘N, ψ ∈ 𝒢(K), kernel {1, η} | CMQ 193–209 | class field theory | correct (D6); \|Cl_K/(1−c)Cl_K\| = \|𝒢(K)\|/2 in 26/26 fields |
| Theorem B (A > 3/2) and (5.2) | CMQ 211–263 | Thm A, Prop C, Lemma 4, HT (m = 1), TZ, ZFR | correct; parameters at A = 1.55 recomputed: φ(σ₀) = 0.96412, ε′ = 0.00897, s₁ = 0.00558, E_b margin 0.0034 at κ = 10⁻³ (D7) |
| Corollary B1 | CMQ 265–271 | Thm B, Remark 2.3 | correct |
| [CM2] Lemma 1 (local orders, criterion C_f, dyadic primes) | CM2 41–57 | DVR lattices | correct |
| [CM2] Lemmas 2–3 (ring class groups, Hecke relations with j_𝔣(𝔭)) | CM2 59–81 | — | correct (split, inert, ramified, 𝔭 \| 𝔣 re-derived) |
| [CM2] Lemma 4: dual identity with ω and R(m) | CM2 83–120 | Lemmas 1–3 | correct: [O_𝔣^× : U] = q\|μ(O_𝔣)\|/2, R(m) well defined mod U² (D8) |
| [CM2] Lemmas 5–6 (detector over F) | CM2 133–156 | [R3] Lemmas 3–4, CMQ Lemma 3 | correct; M_𝔱 run over squarefree relative norms, hence primitive ideals |
| [CM2] Lemma 7 (boxes in ideal lattices: only Π_vZ_v matters) | CM2 180–187 | Π\|k_v\| ≥ N𝔞 | correct; 800 random boxes, skew 10^{±6}: max count/bound 0.167 |
| [CM2] Lemma 8, §4.3 (large 𝔴, source of N³X⁻¹) | CM2 189–201 | Lemma 7, [R3] §3.3 | correct; no longer load-bearing after CM3 |
| [CM2] Lemma 9 (Poisson over the ideal lattice, hyperbolic region) | CM2 212–229 | Lemma 7, Lemma 1 | correct: main-term constant, dual lattice 𝔮^{−1}𝔇_F^{−1}, Jacobsthal factor at 𝔭 \| 𝔠 re-derived (D8) |
| [CM2] Lemma 10 (local densities incl. dyadic) and (S) | CM2 231–252 | Hensel | correct (one vacuous sentence in (b)) |
| [CM2] §4.5 Mellin, Ǩ_n = Φ̌(1−z)B(z,½)ⁿ | CM2 254–272 | ζ_F | correct |
| [CM2] Lemma 11 (Σ_{u∈U/U²} turns [ν_𝔫/𝔮] into ray class characters) | CM2 282–297 | [U:U²] = 2ⁿ | correct |
| [CM2] Lemma 12, §4.7 (classical large sieve suffices) | CM2 299–323 | convexity at Re s = −ε | correct: Q² ≤ L/(NκNκ′) |
| [CM2] Lemma 13 (pole killed by (Z)), Lemma 14 (Jutila over F) | CM2 327–371 | quadratic Hecke characters | correct |
| [CM2] Theorem 1 (family sieve, X^{3/2} + N³X⁻¹) | CM2 160–165, 373–382 | above | correct, for the detector sequences only |
| [CM2] Theorems 2–3 ((13−8σ)/(7−2σ); A > (3+√69)/10) | CM2 388–434 | Thm 1, Prop C | correct: σ_c = 0.846688, A₀ = 1.130662 (D10) |
| [CM3] Lemma P_F (relative Lemma P over all orders O_F + 𝔣O_K) | CM3 20–67 | CM2 Lemmas 2–3 | correct (D9); lattice division at 𝔭 \| 𝔢 needs v_𝔭(𝔠𝔯) = 0, which l.40 and l.43 ensure |
| [CM3] Corollary P_F | CM3 69–79 | Lemma P_F, Φ ≥ 0 | correct: the identity is an interchange of finite sums, positivity gives the inequality |
| [CM3] Theorem 1 (family sieve X^{3/2+C₁δ}N^{1−2σ}, N ≤ X^{1−4δ}) | CM3 83–102 | CM2 §4.4–4.8 at 𝔴 = O_F | correct: N𝔠√M₀ ≤ 4N, E₀ ≥ 2v^{−1/2}, E₀² ≤ L |
| [CM3] Theorem 2 (Σ_KN_K(σ) ≪ X^{2−σ+ε}) | CM3 108–117 | Thm 1, CMQ Lemma 3 | correct |
| [CM3] Theorem 3 (every A > 1) and sharpness remark | CM3 121–156 | Thm 2, Prop C, CMQ Remark 2.3 | correct: φ = 1 − (A−1)(1−σ) ≤ 1 − 8ε′ |
| [CM3] §5 (every place where F ≠ Q changes) | CM3 158–173 | — | correct |

## Findings (none breaks a proof)

1. **σ-range (discovery file).** Claim 2 of discoveries/cm-least-primes.md has no σ-range; CMQ Theorem A(a) is for 0.6 ≤ σ ≤ 1 (l.140). For σ ≤ ½ the claim follows from Theorem A(b), because 2(1−σ)/(3−2σ) ≥ ½ there; on ½ < σ < 0.6 it is not proved (the [C] Lemma 5 detector needs σ ≥ 0.6) but routine with CM2 Lemma 6 (σ ≥ ½ + 2δ). Fix: add "0.6 ≤ σ ≤ 1".
2. **Wording (discovery file).** "Rough" means odd δ₀-rough (CMQ §4 family); a prime relative discriminant is automatically odd. The κ of the d_K-form is half the NOTE's κ (l.219–220).
3. **Regulator.** "No regulator factor of F" holds for the strengthened Lemma 2 (l.299–308). The NOTE's Lemma 2 and 3 statements (l.79–84, 118–121) still carry e^{nr_F}: true, weaker. For Q(√5) the strengthened bound is 48× the true count at large N, as predicted by (4π)²d_F/((2π)²R_K/(w_Kζ_F(2))).
4. **Lemma 3, "both terms attained up to F-constants".** The N-term is attained up to κ_K = d_K^{o(1)}, not up to a constant depending only on F.
5. **Tests in the NOTE.** CMQ §6 tested Lemma 1 only for 𝔠 = O_K. Here: 550 random ideals 𝔠, 455 of them non-principal; all pass.
6. **Effective but astronomical X₀** (as in RC5Audit): the Poisson step needs X^{3δ} ≥ 4C₃ with δ ≤ (A−1)/800. No numerical least-prime bound follows.
7. **Not re-refereed:** [R] §2 bookkeeping (margins only, as in RC5Audit), and the proofs of HT and TZ (statements checked in the TeX).

## Numerical falsification tests (all passed)

**Per field.** Own PARI/GP 2.19 code ([checks/cmls.gp](checks/cmls.gp)): K = F(√−α), class groups by bnfinit, conjugation c from nfgaloisconj, 𝔇_{K/F} = 𝔇_K/f′(y), primes grouped into c-orbits, all ideals of norm ≤ N_max enumerated with class, O_F-content and Z-content. 22 fields: Q(√5) 13 (9 with prime D, 2 with D = p₁p₂, 2 with even 𝔡), d_K = 1.0·10⁴ … 1.0·10⁸, h_K = 3 … 239; Q(√2) 6, d_K = 2.8·10⁴ … 9.9·10⁷; Q(√10) 3 (h_F = 2, so 𝒵_K ≠ Ĉl_K), d_K = 1.7·10⁶ … 9.9·10⁷, h_K ≤ 180. N/√d_K ∈ {1/8, …, 16}; N ≤ 1.6·10⁵, up to 2.1·10⁵ ideals per field. Peak memory 51 MB, single-threaded ([out_q5.txt](checks/out_q5.txt), [out_q2.txt](checks/out_q2.txt)).

| Quantity | Result |
|---|---|
| Lemma 1, Δ ∈ 𝔠𝔠ᶜ𝔇_{K/F} | 3300/3300 pairs; \|NΔ\|/((N𝔠)²D) an integer ≥ 16 in the samples; product identity error ≤ 6·10⁻³⁹ |
| Lemma 2, max_Cπ_C(N)/(1 + 4π(N/√D)^{1/2})² | 0.006–0.042 (all fields and N) |
| Lemma 3 sharp constant K(N)/(N + √d_K), ε = ½ | 0.026–1.66 |
| K(N)/((N + √d_K)κ_K), κ_K = res ζ_K ∈ [0.29, 2.77] | 0.077–0.877; at N = 16√d_K it is an F-constant, ζ_F(4)/ζ_F(2)·(max/mean): 0.82–0.88 (Q(√5)), 0.70–0.74 (Q(√2)), 0.60–0.62 (Q(√10)) |
| slope of log K/(N + √d_K) in log d_K at fixed N/√d_K | −0.013 … +0.029 (no growth) |
| growth in N, N = √d_K/8 → 16√d_K | factor 3.6–9.1, saturating: K ≈ c₁h_K + c₂κ_KN, both terms present |
| ε = 0.1 instead of ½ | within 2.4% (Q(√5)), 7% (Q(√2)), 12% (Q(√10)) |
| 𝒵_K constant (Q(√10), \|j(Cl_F)\| = 2) | ≤ the Ĉl_K constant |
| unweighted / weighted constant | 1.0–3.83 (largest at d_K > 10⁷, small N); Z-primitive / weighted ≤ 2.52 |
| exact Gram diagonalisation | h = 11 (Q(√5)) and 36 (Q(√10)): top eigenvalue = h·max_CΣw⁻¹ to 10⁻³⁵ |
| CM2 Lemma 7, 800 boxes in ideals of Q(√5), Q(√2), Q(√10) | max count/bound 0.167 |
| CMQ Lemma 4, 26 fields over Q(√5), Q(√2), Q(√3), Q(√10) | \|Cl_K/(1−c)Cl_K\| = \|𝒢(K)\|/2 in 26/26 ([out_genus.txt](checks/out_genus.txt)) |

Reading: the content-weighted constant is ≍ κ_K(N + √d_K) with no power of d_K beyond κ_K, as Lemma 3 states (its log^{2n−1}(3d_K) is generous, as in RCAudit's Lemma 2 test). The weight is needed: without it the constant grows with d_K at N < √d_K (Remark 2.2).

**Family sieve (CM3 Theorem 1)** ([checks/famcm.gp](checks/famcm.gp), [out_fam.txt](checks/out_fam.txt)). F = Q(√5); all CM K/F with X/2 < N𝔡_{K/F} ≤ X (one per squarefree ideal (α), α ≫ 0; the X = 2000 count, 91 fields with Σh_K = 455, matches CM2 §7). b = the detector sequences of CM2 §3: admissible 𝔫 (squarefree, prime factors of norm ≥ 19), a_V with V = √X/2, σ = 5/6, Y = X, pole killed by the 4^{1−σ}-term (or not). LHS = Σ_Kh_KΣ_C|S_C|² = Σ_KΣ_χ|Σ_𝔫b(𝔫)λ_χ(𝔫)|²; R1 := LHS/(Σ_Kh_K·Σ|b|²), the family-size normalisation; R_X := LHS/(X^{3/2}Σ|b|²), Theorem 1's.

| X | fields | Σh_K (= 0.0051X^{3/2}) | pole killed, N ≥ 4√X: R1 (R_X) | χ = 1 share | pole present: R1 at N ≈ X/2 |
|---|---|---|---|---|---|
| 1000 | 45 | 157 | 0.85–1.16 (0.0042–0.0057) | 20–65 % | 38.6 |
| 2000 | 91 | 455 | 0.88–1.38 (0.0045–0.0070) | 25–52 % | 47.2 |
| 4000 | 186 | 1295 | 0.94–1.64 (0.0048–0.0084) | 13–45 % | 55.4 |
| 8000 | 367 | 3680 | 0.94–1.51 (0.0048–0.0078) | 27–34 % | 80.6 |
| 16000 | 740 | 10318 | 0.96–1.57 (0.0049–0.0080) | 17–39 % | 100.2 |

- Pole killed: slope in N on [4√X, X/2] between −0.33 and 0.00; slope in X at fixed N/X = 1/16, 1/8, 1/4, 1/2: −0.07, +0.07, −0.03, +0.15 (the last is N = X/2, outside Theorem 1's range N ≤ X^{1−4δ}). LHS stays at the family size, as Theorem 1 says.
- Pole present: slope in N 0.66 → 0.92, R1 up to 100: the trivial character, about N/√X times the family size (CM2 §4.9 Remark (1)).
- N < 4√X (where a_V removes most of b): pole-killed R1 = 0.9–5.0, flat in X, with χ = 1 carrying 19–90 %, as for F = Q in RC5Audit.
- Peak memory 127 MB; X = 16000 took 7.7 min.

## Literature check

- `lst.py 'CM field|class group|least prime|large sieve|Hecke|ideal class|zero density|totally real|Chebotarev' --zenodo '"CM field" "least prime"'`: 45 arXiv title hits (Oct–Nov 2026, 7 categories), 17 on class groups or Hecke L-functions (e.g. 2611.00924, 2610.04045, 2610.15989, 2610.18491, 2610.26835), none on least primes in ideal classes or class-group large sieves or densities over families of fields. Zenodo: 117 loose hits, none relevant.
- Pascadi–Thorner 2508.14888 (TeX v4): Thm 1.3, (D_F^{−n(n+1)}Q^{2n+1}T^{[F:Q]n(n+1)+4}\|S\|⁴)^{(1−σ)/(3−2σ)+o(1)}; Thm 1.1, (NQ)^{o(1)}(N + Q^n\|S\|), Q^n replaceable by max C(π×π̃′)^{1/2}. For Ĉl_K: d_K^{3(1−σ)/(3−2σ)} and N + d_K. Weaker than (1)–(2); no class-group, CM or least-prime application.
- Cho–Lemke Oliver–Zaman 2512.24963 ("The least prime with a given cycle type"): fixed group G, Np ≪ \|Disc K\|^{α(G,C)+ε}; nothing on ideal classes, CM fields or families.
- Thorner–Zaman 1604.01750: Thm 3.2 and Thm 1.1 (D_K^{694}) as quoted. Lemke Oliver–Thorner–Zaman 2012.14422: almost-all Chebotarev for fixed G only; Hilbert class fields (Galois group of order 2h_K → ∞) not covered.
- Humphries–Thorner 2103.05634: base field arbitrary; θ_m = 0 under GRC, as CMQ uses.
- Miller–Peckner 1003.5336: 1-level density of {L(s,χ)} over CM K with K⁺ = K₀ fixed (h(K₀) = 1), under GRH, supp φ̂ ⊂ (−1,1); their key lemma is the sharpness bound (known, as CM3 NOVELTY2 says).
- Verdict: none of (1)–(3) is stated, implied or contradicted; GRH implies (2) and (3), so no unconditional result can contradict them. Consistent with CMQuartic NOVELTY and CM3 NOVELTY2.

## Reviewer

Two `reviewer` runs ([checks/reviewer1.md](checks/reviewer1.md), [reviewer2.md](checks/reviewer2.md)).
- Run 1 (condensed excerpt, five items): Lemma 1 and Lemma 2 "CORRECT / ROUTINE".
  - G1: Lemma 3 needs ε ≤ 1.
  - G2: a counterexample to Lemma P_F step (b) with 𝔠 = 𝔣 = (3).
  - G3: the 𝔴 = 1 identity is not shown.
- Run 2 (the NOTE lines verbatim): "G1 CLOSED — Lemma 3 explicitly assumes 0 < ε ≤ ½"; "G2 CLOSED — The proposed counterexample has zero weight. For surviving terms, v_𝔭(𝔠𝔯) = 0"; "G3 CLOSED — … finite-sum interchange gives the exact identity … Lemma P_F supplies positivity"; "Overall CORRECT".

RH relevance: N11. This audit removes the last unrefereed link of the CM least-prime line (almost all CM fields over a fixed F, every A > 1, sharp). Records only; no RH node moves.
