AUDIT PASSED: Theorem 1 (hard-core augmentation of Lamzouri's Thm Main2) holds exactly as written, given Lamzouri's audited Prop SmallSetZeros, and the λ = 0.768 certificate is valid (margin 1.46·10⁻³; 7.0·10⁻⁴ by an independent second route); no counterexample on 22,490 zeros or on adversarial hard-core sequences.

# HCAudit (Orca audit, N8, falsification first) — 12 Nov 2026, 17:20–17:55 UTC (agent hcaudit)
Target: discoveries/small-gap-hardcore.md (DISCOVERIES N8), Riemann/Deep/ExplorePointProc (explore40). Proofs: [FULLPROOF.md](FULLPROOF.md).
Code: hc_indep.py, hc_rig.py, hc_zeros.py, hc_synth.py, hc_mf.py. Logs: logs/. Source: Lamzouri arXiv 2611.14486v1 TeX (l.n = its lines).

## Short answer
- Theorem 1 is correct as stated; FULLPROOF §1–§4 is a complete proof whose only external inputs are Lamzouri's
  (Eq:ReductionGoodSet2), i.e. Prop SmallSetZeros (l.331–339; lamzaudit PASSED; Conrey density via Baluyot, BGST Lemma 5), and
  Riemann–von Mangoldt.
- cert768 is valid; D(λ) > 0 for every λ ≥ 0.768 stands. No step fails, so no weaker λ needs recomputing.
- Cosmetic errata only (§4). The proportion is ineffective: κ = 1/16 − 7θ0/64 = 1.1·10⁻⁶, so δ is astronomically small (as for Lamzouri).

## 1. Proofs (FULLPROOF.md)
| step | verdict | load-bearing facts checked |
|---|---|---|
| (i) reduction to Z_0, lower bound (§1) | complete | D(λ), E_λ, Z(T), N(T) all use zeros in (0, T] and FIXED L = log T (l.178, 183, 298, 405); an off-line zero's mirror has the same ordinate, so off-line and multiple zeros lie in E_λ; E_λ is invariant under ρ → 1 − ρ̄; Prop SmallSetZeros' hypotheses hold for g = ψ_εF and ψ_εf (even, ≥ 0, C_c², supp ⊂ [−θ', θ'], θ' = (θ + 4/7)/2); its proof (l.340–390) uses nothing else and its error C_fδ^κ is uniform in δ; the positivity identity holds for every even f ∈ L¹ (Cauchy convolution 2π/(4 + d²)) |
| (ii) step (iii′) (§2) | complete | ordered pairs: diagonal once, k = 1, 2 twice (f^ even), matched to the N_0 − 2 applications of (H) (×2 orders) + one leftover; k ≥ 3: distance > 3λ ≥ U strictly (gaps > λ in the same unfolding as f^), dropped; weights cost ≤ 4π²M_fN_0/L²; A telescopes to A(g_{N_0−1}) − A(g_1); β-term 2β(x_{N_0−1} − x_1 − N_0) ≤ 2β(#E_λ + 4N/L): needs exactly β ≥ 0; N_0(f^(0) + 2c) handled for either sign; contradiction for every large T; #small gaps ≥ #E_λ/2 gives D(λ) ≥ δ/2 |
| (iii) Lemma T′ (§3) | correct, redundant | F^ − f^ ≤ (ε_R − c_rκλ⁶)ξ⁻⁶ < 0 on [λ, ∞): (H) passes to F exactly (no o(1)); G^ ≤ −κ(λ/ξ)⁶ re-proved; R^ = FT R re-verified; J_θ(F) ≥ J_θ(f), F^(0) ≤ f^(0) + 2c_r. Not needed: Lamzouri's l.413–441 accept Schwartz f directly |
| certificate scheme (§4) | correct | Q1 = [λ, U]² bilinear (A-part exact on node-aligned cells; error (δ²/8)(2M2 + M2)); Q2 = y ≥ U, Q3 = x ≥ U by sign of f^ and the constant tail; Q1 ∪ Q2 ∪ Q3 = [λ, ∞)² |
The brief's [T, 2T] / varying-log γ concern does not arise: Lamzouri never uses the local density; only the total range TL/2π = N(T)(1 + O(1/L)) enters.

## 2. Falsification (FULLPROOF §6; none succeeded)
| test | data | result |
|---|---|---|
| (iii′) chain on actual zeros | 22,490 zeros ≤ T = 19998.76, Lamzouri's E_λ (#E = 26.2%), N_0 = 16604 | min gap 0.768032 > λ; k ≥ 3: all f^ ≤ 0; max Φ(g_i, g_{i+1}) = −0.347999 ≤ c_cert; telescoping exact to 9·10⁻¹²; S_w = 9284 ≤ bound 24419; S_w/N(T) = 0.413 ≥ J = 0.307 (the lower bound's direction; zeros escape through #E_λ, β-term 19,345) |
| exact finite inequality | 8 adversarial families (Tonks, all gaps λ, λ-runs + long gap, alternating, gaps at A's nodes, sup pair repeated, uniform, end clusters) | max(S − bound) < 0 everywhere (closest −0.158) |
| ground state of (H) | cyclic sequences n ≤ 30, coordinate ascent + restarts | max E_β = −0.3479921 at g ≡ λ, 1.3·10⁻⁴ below c_cert |
| sup Φ | 8× finer Arb lattice + local refinement | −0.3478845702 at (λ, 1.11709); c_cert above by 2.3·10⁻⁵ (best possible margin for this f, A, β: 0.001508) |

## 3. Certificate and code audit (FULLPROOF §5)
| check | result |
|---|---|
| reruns on /tmp copies (cert768, FULLARB=1, xcheck; cert770, cert775; shipped peaks 41 MB / 23 MB read first; neither script writes files) | identical to the shipped logs: margins 0.001462588425, 0.005142294469, 0.01374683343; all PASS |
| f > 0, f^ < 0 on [U, ∞), independent of Sturm | Arb root isolation: p's only real root −0.000218; q's largest real root 40.357846 < w0 = 40.358051 < 2πs²U² (f^ changes sign 6·10⁻⁶ below U, harmless: k ≥ 3 distances > U strictly) |
| M2 (checker: Cramér + Hermite ladder, M3 = 3536.5) | math re-derived; numerical max|f^'''| = 30.2: valid, loose |
| second route (no Cramér): M2 ≤ 4326.5 by Arb ball cover, corners on a dg/40 lattice in Arb | c = −0.347481068, margin 7.014·10⁻⁴ > 0 PASS |
| J | rigorous acb integration [0.306839281978626 ± 1.8·10⁻¹⁶] = closed form |
| cert_aug.py line by line | Hermite recurrence, p/q extraction, J and f^'' formulas, Sturm (gcd-terminated; endpoints not roots), F2 indexing k = i + j, node-aligned A interpolation, constant tail in R2/R3, M2 grid covering [λ, 2U]: all correct |

## 4. Errata (cosmetic; none affects validity)
- discoveries/small-gap-hardcore.md: "c_cert ≈ −0.3478618" mis-rounds −0.3478616532 (7 decimals: −0.3478617); the displayed number is
  1.5·10⁻⁷ below what the checker certifies (still a true bound on sup (H) ≈ −0.3478846, but not the certified one).
- cert_aug.py l.94 comment "+ 1e-15" is not implemented (float rounding of A, ~10⁻¹⁶, is covered by the final +10⁻¹² and FULLARB).
- Statement: "λ < U" is never used, and "A = β = c = 0, U = λ recovers Thm Main2" needs any U ∈ (λ, 3λ] instead of U = λ.
- Hypothesis can be weakened: f even, ≥ 0, L¹, C² near [−θ', θ'], sup u²|f^(u)| < ∞ (covers C_c² and Schwartz).
- Watchdog: xcheck_aug.py peaks at 266 MiB RSS (its last line builds 43 Hermite arrays of 400,001 floats); memcap's ~3 s sampling
  saw 23 MB and reports 0 MB for runs under one sample. Under the 300 MB cap, but short bursts escape the watchdog.

## 5. What Theorem 1 rests on
Lamzouri 2611.14486v1 Prop SmallSetZeros and l.413–441 (unrefereed; lamzaudit PASSED line by line; Conrey's slope-8/7 density only
in Baluyot's thesis), Riemann–von Mangoldt (only TL/2π − N(T) = O(N/L)), and the certificate's exact/Arb facts (re-verified by two
routes). Step (iii′) and the certificate scheme are elementary and now fully written out and checked.

## 6. Reviewer (verdict-only, step (iii′) with Lamzouri's exact definitions and the certificate values inline; logs/reviewer_iii.log)
CORRECT for multiplicities, the k ≥ 3 drop at U = 3λ, weights, telescoping, β-sign, fixed-L normalisation, step (5) and the charging
step; one GAP in this audit's own write-up: the explicit cutoff "≤ 4N/L for L ≥ 10" in §2(4) is unsupported by the O(log T)
remainder → replaced by "for all large T" (ratio → 1 + log 2π < 4). It recomputed the margin from the 10-digit display: 0.00146258838
(vs full-precision 0.001462588425); positivity unchanged. A second reviewer (proofread of all of FULLPROOF.md) timed out at
115 s without output; the proofread was done in-agent (fixed: multiplicity wording, a section reference, M_f = 0.776 in §6(a)).

## 7. Process
One process at a time, single-threaded, memcap 290 MB. Peaks (ru_maxrss / memcap sample): cert768 43/20 MB, FULLARB 46/22, xcheck
266 MiB/23 MB, cert770 41/0, cert775 41/19, hc_indep 47/0, hc_zeros 42/0, hc_synth 54/32, hc_rig 35/14, hc_mf 77/0. Each run ≤ 42 s.
Downloads (/tmp/hcaudit) deleted. No edits outside Riemann/Deep/HCAudit/ and .scratch/hcaudit/.

RH relevance: N8
