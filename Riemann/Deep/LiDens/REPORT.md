RECORD: λ_n > 0 for every n ≤ 110.4757·T², T = 3000175332800 (was 106.2557 T²)

# lidens (node N3, records only): Kadiri–Lumley–Ng's zero density inserted into Li200k's certificate, 2026-11-20

## Short answer

- **Certified** (ball arithmetic, python-flint 0.9.0): with the Lean-checked half-plane Re s > θ₁₄₁ = 0.8733458561067185768…, RH to
  T = 3,000,175,332,800 (Platt–Trudgian), HSW's N(t) bound and KLN's explicit N(σ, T′) at H0 = 3·10¹² (rows 0.85–0.88),
  **λ_n > 0 for all 1 ≤ n ≤ 994,397,608,026,097,376,944,979,968 = 110.475709 T²** (fine grid; the coarse grid gives 110.460213 T²).
  +4.0% over Li200k's 106.255746 T² (reproduced: chain 1 of every run).
- **The brief's literal bound** (2·KLN simplified form as the count, slabs from 0.845): 109.947461 T². Two rigorous refinements add
  +0.5%: the count is min(2·KLN, HSW total count), and KLN's own log-form (Lemma 4.14, (bnd-N@H0)) replaces its simplified form.
- **Second row, B = 436666397/500000000** (QRT PR 144, pending): λ_n > 0 for n ≤ 994,414,124,818,031,229,790,060,544 = 110.477544 T²
  (fine grid; coarse 110.462046); B moves by 1.3·10⁻⁵, worth +0.002 T².
- **Robustness**: KLN Thm 1.1's printed "δ ≥ 1" enforced (δ = 1 in every row): 108.918572 T² (still a record).
- **Falsification passed**: Li200k's per-zero inequality re-derived and ball-checked at the real scale (400 zeros, γ ∈ [T, 3T], δ ≤ B − ½,
  n ≤ 120 T²; min relative slack 3.4·10⁻⁶: the bound is nearly tight, as expected when cos nφ ≈ 1); q_n and every slab difference f_j
  decrease in t (proof §1, 200 ball checks); KLN counts exactly what the sum needs (β > σ strict, γ < T strict, multiplicity; slab
  convention (a_j, a_{j+1}] matches). Reviewer: CORRECT given a_0 ≥ 0 (here a_0 = 0.30). NOGO l.546 (R400: "explicit zero density cannot
  raise it") is now refuted by a certificate (CeilAudit3 B2 had refuted it in floats).
- **What binds now**: the density near σ = B. At N the loss is 99.8% the slab part (top slab, σ_{J−1} within 10⁻⁴ of B: 2·KLN(σ_{J−1}, T + 10¹⁰)
  = 1.7·10¹⁰ zeros allowed in (T, T + 10¹⁰]), 0.2% the all-zero part at a_0 = 0.30. At N, d log Q/dc = 0.364 and d log M/dc = 0.009, so a factor κ
  less density loss is worth at most ln κ/0.355 T² (halving: +1.95 T²).

## 1. The inequality (every input with its source line)

Notation: ρ = β + iγ zeros with multiplicity, T = T0 = 3000175332800, δ = |β − ½|, g(x) = 2 − 2cos x, θ(t) = 2 arctan(1/(2t)),
ψ_n(t, a) = 4 sinh²(na/(2t²)) + 4na²/t³ (Li200k's q_n with a(t) ≡ a: REPORT l.45, verify.py l.135–138).

1. **Per zero** (Li200k REPORT l.42–48, re-derived in NOTES.md Step 1): for ρ = ½ + δ + iγ and ρ* = 1 − ρ̄ (same height and multiplicity),
   z = 1 − 1/ρ: |z_ρ| = r, |z_ρ*| = 1/r, ℓ = −log r ≤ δ/(γ² + (½ − δ)²) ≤ δ/γ², arg z_ρ = arg z_ρ* = φ, |φ − θ(γ)| ≤ 2δ²/γ³; per upper zero
   2 − 2cosh(nℓ)cos nφ ≥ g(nθ(γ)) − ψ_n(γ, δ). So λ_n ≥ M_n − L_n, M_n = Σ_{γ>0} g(nθ(γ)), L_n = Σ_{γ>T, upper zeros} ψ_n(γ, δ).
2. **Monotonicity** (analytic): ψ_n ↑ in a ≥ 0 and n, ↓ in t; for a′ > a ≥ 0, f = ψ_n(·, a′) − ψ_n(·, a) ≥ 0 is ↓ in t and ↑ in n
   (s = n/(2t²): d/ds[sinh²(a′s) − sinh²(as)] = a′ sinh 2a′s − a sinh 2as > 0; 4n(a′² − a²)/t³ likewise).
3. **Slabs**: a_0 < … < a_J dyadic, σ_j = ½ + a_j, σ_J = B⁺ ≥ θ₁₄₁ (no zero has δ > θ₁₄₁ − ½: half-plane + reflection). Cap a(ρ) = a_0 if
   δ ≤ a_0, a_{j+1} if a_j < δ ≤ a_{j+1}. For γ ≤ t_K: ψ(γ, δ) ≤ ψ(γ, a_0) + Σ_j f_j(γ)·1[δ > a_j], f_j = ψ(·, a_{j+1}) − ψ(·, a_j);
   for γ > t_K: ψ(γ, δ) ≤ ψ(γ, a_0) + ψ(γ, a_J). Hence **L_n ≤ Qall_n(a_0) + Σ_j S_j(n) + R_n** with:
4. **Qall_n(a_0)** = Li200k's Q_hi with a ≡ a_0 (verify.py l.140–172: 2ψ(T)E(T) + upper sum of ∫_T^{100T} ψ w + tail_bound l.125–133),
   valid by the Stieltjes step Σ_{γ>t1} q(γ) ≤ 2q(t1)E(t1) + ∫_{t1}^∞ q(Ñ′ + E′) (Li200k REPORT l.49–50), E = HSW Cor 1.2 with 9.4926 + 7/8
   (verify.py l.25–36; LiRangeAudit), w = (1/2π)log(t/2π) + 0.13/t ≥ Ñ′ + E′.
5. **S_j(n)** = f_j(t_K)U_j(t_K) + Σ_k U_j(t_{k+1})(f_j(t_k) − f_j(t_{k+1})) ≥ Σ_{T<γ≤t_K, δ>a_j} f_j(γ), because
   C_j(t) := #{T < γ ≤ t : δ > a_j} = 2N(σ_j, t⁺) (β > σ_j and the reflections β < 1 − σ_j; none with γ ≤ T), C_j(T) = 0,
   Σ = f_j(t_K)C_j(t_K) + ∫_T^{t_K} C_j(−f_j′) dt, C_j ↑, −f_j′ ≥ 0, and C_j(t) ≤ U_j(t) := min(2·KLN(σ_j, t), N⁺(t)),
   N⁺(t) = Ñ(t) − Ñ(T) + E(t) + E(T) ≥ N(t) − N(T⁻) (HSW).
6. **KLN(σ, t)** = (t − H) log t/(2πd)·log(1 + 2πd·A(σ)(log kt)^{2σ}(log t)^{4(1−σ)}t^{8(1−σ)/3}/(t − H)) + B(σ)(log t)², H = H0 − 1,
   H0 = 3·10¹², k = 1 (KLN arXiv 2101.12263 Lemma 4.14 (bnd-N@H0), TeX l.1548–1566; Thm 1.1 l.193–217; N(σ,T) defined l.121:
   σ < β < 1, 0 < γ < T; multiplicity via Littlewood's lemma l.259–261). A(σ) = exp(c0 + c1σ), B(σ) = q0 − q1σ from Ingham4/cert.py
   class KP l.120–138 and build_kln_H0('3e12') l.274–292 (admissibility: R2AAudit, R500 (1)); valid for σ > ½ + d/log H0 = 0.5117
   and every t ≥ H0; per slab the row of 0.85–0.88 smallest at t = T. Rows reproduce the brief (logs/kln_rows.log).
7. **R_n** = 2ψ(t_K, a_J)E(t_K) + tail_bound(n, a_J, t_K) ≥ Σ_{γ>t_K} ψ(γ, a_J) (same Stieltjes step from t_K = 2T).
8. **Blocks** (Li200k REPORT l.58–59, verify.py l.248–260): on [n1, n2], λ_m ≥ M_lo(n1, n2) − (bound at n2) (all ψ terms ↑ in n);
   chain 1 = Li200k's chain with a = a_J from 89; chain 2 continues with the new bound; n ≤ 88 by Li200k part 1 (REPORT l.56–57; rerun, --small).

## 2. Falsification checks (brief step 1)

| check | result | where |
|---|---|---|
| per-zero inequality, real scale: 400 balls ρ = ½ + δ + iγ, γ ∈ [T, 3T] (¼ at γ = T), δ ∈ (0, θ₁₄₁ − ½] (⅓ at the cap), n ≤ 120 T², prec 480 | holds, min relative slack 3.42·10⁻⁶ | cert.py falsify(), logs/cert_B141*.log |
| ψ_n(·, a) and f = ψ_n(·, a′) − ψ_n(·, a) decreasing in t, n ≤ 120 T², t ∈ [T, 2T] | 200 ball checks pass; analytic proof §1.2 | same |
| q_n(t) decreasing on [T, ∞) for every n (worst-case placement at T) | true for all n ≥ 1, t > 0 (both terms ↓ in t) | §1.2 |
| KLN's count vs the sum's need | N(σ,T): σ < β < 1, 0 < γ < T, multiplicity (TeX l.121, 259–261); sum needs δ > a_j, T < γ ≤ t: = 2N(σ_j, t⁺), bound continuous in t; zeros with δ = a_j exactly sit in slab (a_{j−1}, a_j], charged at a_j | §1.5 |
| KLN Table 1 constants hold only at σ = σ0 (TeX l.1690: "for σ = σ0"), not for σ ≥ σ0 as the brief says | we evaluate each row's A(σ), B(σ) at the slab's own σ (Lemma 4.14 holds for σ > 0.5117 with σ-dependent constants); A(σ) rises with σ (c1 ≈ 5.0–5.1), so the lemma does not give a fixed A(σ0) on [σ0, B] (R500's float model used it) | §1.6 |
| KLN Lemma 4.14 at H0 = 3·10¹² (printed with H0 = 3.06·10¹⁰) | input as audited (Ingham4 reading, Cully-Hugill–Johnston precedent; R2AAudit) | caveat |
| KLN Thm 1.1 prints δ ≥ 1 (rows use δ ≈ 0.303) | R2AAudit: no step uses δ ≥ 1; enforced anyway: 108.918572 T² | logs/cert_B141_delta1.log |
| reproduction of the old record | chain 1 of every run: 106.255746 T² (B = θ₁₄₁), 106.259560 T² (PR 144) | logs |
| Li200k part 1 (n ≤ 88) with the new variants | λ_n ≥ 0.00335 − 2.1·10⁻²³ | logs/cert_B141_best.log |

## 3. Reviewer (brief step 4)

Verdict-only reviewer subagent, steps 1.3–1.8 pasted inline with exact formulas (slab caps, factor 2, Stieltjes with C_j(T) = 0, min of
the two counts, blocks): "GAP: (b)–(e) require a_0 ≥ 0. Otherwise (d)'s factor-two identity fails, and f_j need not be nonnegative or
decreasing, invalidating (e). With a_0 ≥ 0, (b)–(g) are valid under the supplied counting bounds; endpoints, reflections, and
multiplicities are handled correctly." All runs have a_0 = s* − ½ ∈ {0.30, 0.33, 0.345} > 0 (asserted in cert.py): **CORRECT** for the
certificate as run.


## 4. Certificate (brief step 3; code/cert.py, code/runs.sh; one thread, ball arithmetic at 160 bits, real-scale checks at 480 bits)

| run (log) | B | count | slabs s* … B | t-grid (to 2T) | N (λ_n > 0 for all n ≤ N) | N/T² | time, peak |
|---|---|---|---|---|---|---|---|
| **B141_best** | θ₁₄₁⁺ | min(2·KLN_log, N⁺) | 0.80 (h 0.002) / 0.84 (h 10⁻⁴), J = 354 | K = 388 | **994,397,608,026,097,376,944,979,968** | **110.475709** | 175 s, 66 MB RSS |
| B144_best | 436666397/5·10⁸ | same | same | same | 994,414,124,818,031,229,790,060,544 | 110.477544 | 196 s, 66 MB |
| B141 | θ₁₄₁⁺ | min(2·KLN_log, N⁺) | 0.83, J = 176 | K = 172 | 994,258,122,450,992,374,040,494,080 | 110.460213 | 70 s, 40 MB |
| B144 | 436666397/5·10⁸ | same | same | same | 994,274,622,515,645,990,570,033,152 | 110.462046 | 77 s, 39 MB |
| B141_plain (brief's literal step 2) | θ₁₄₁⁺ | 2·KLN_simplified | 0.845, J = 189 | K = 172 | 989,642,813,205,584,442,678,050,816 | 109.947461 | 80 s, 40 MB |
| B141_delta1 (δ = 1 in every row) | θ₁₄₁⁺ | min(2·KLN_log, N⁺) | 0.83, J = 176 | K = 172 | 980,381,729,997,521,889,848,721,408 | 108.918572 | 60 s, 39 MB |
| Li200k baseline (chain 1) | θ₁₄₁⁺ / PR 144 | all zeros at a_B | — | — | 956,413,498,941,653,124,095,410,176 / 956,447,831,274,852,956,086,730,752 | 106.255746 / 106.259560 | 17 s |

θ₁₄₁⁺ = ⌈0.873345856106718576801184797·2⁶⁴⌉/2⁶⁴ (QRH decimal 0.873345856106718576801184796…). Each run: chain 1 (Li200k's own blocks
from 89 with a = a_B), then chain 2 (blocks with the new bound; 9–13 blocks, 132–196 bound evaluations), at the final N:
M_lo = 2.9783957·10²⁸ > Q = 2.9783857·10²⁸ (Qall(a_0) 5.05·10²⁵, slabs 2.9733·10²⁸, tail 2.7·10¹⁷). Float model (code/model.py,
code/scan.py, logs/scan.log; not a certificate) predicted 110.49 (refined, s* 0.80–0.83) and 109.95 (literal, s* 0.845); the grids lose
≈ 0.015 T². KLN rows: logs/kln_rows.log. Rerun: bash Riemann/Deep/LiDens/code/runs.sh all (≈ 11 min; peak 66 MB).
The final code was rerun end to end (runs.sh all, 07:10–07:21 UTC): all six logs identical up to run times.

## 5. Limits

- Only the constant moves: the loss is ≈ (effective count near T)·e^{c(B − ½)}; Li200k's effective count is ≈ 1.6·10¹¹ zeros of
  [T, 1.0126T] (CeilAudit3 B2), the slab sum acts like 1/4.8 of it (e^{(110.476 − 106.256)(B − ½)} = 4.8), so c gains ln 4.8/(B − ½).
- Not used: KLN rows re-optimised at σ ≈ B or at H0 = T (A(σ) of the rows used is ≈ 9.6 at B); a short-interval density near T
  (CeilAudit3 l.78: 2610.15989, unread); any lower half-plane than θ₁₄₁.
- Inputs as before: Platt–Trudgian (RH to T, refereed), HSW Cor 1.2 with 9.4926 (refereed + misprint fix), θ₁₄₁ (Lean-checked,
  unrefereed), KLN Lemma 4.14 at H0 = 3·10¹² (refereed paper; our evaluation at the higher RH height, audited in R2AAudit).
- Novelty (records-level, not rechecked): R500 (4) found no paper combining an explicit zero density with Li-positivity ranges.

## Suggested RECORDS.md line

| N3 Li coefficients | λ_n > 0 for all n ≤ 9.9439·10²⁶ = 110.4757 T² (T = 3000175332800) | θ₁₄₁ half-plane + KLN zero density (H0 = 3·10¹², log-form) + HSW total count in Li200k's loss sum; δ = 1 rows: 108.92 T²; PR 144 B: 110.4775 T² | lidens, 20 Nov 2026 (was 106.2557 T², Li200k) |

RH relevance: N3 (records only; ranges carry no information on zeros above T).

## Files

- Riemann/Deep/LiDens/REPORT.md — this report; NOTES.md — derivation notes (per-zero inequality, KLN count, slab bound).
- code/cert.py — certificate (imports Li200k/verify.py and Ingham4/cert.py read-only); code/runs.sh — all runs; code/kln_rows.py — rows;
  code/model.py, code/scan.py — float model and s* scan.
- logs/cert_{B141_best,B144_best,B141,B144,B141_plain,B141_delta1}.log, logs/kln_rows.log, logs/scan.log, logs/model_085.log, logs/test1.log.

