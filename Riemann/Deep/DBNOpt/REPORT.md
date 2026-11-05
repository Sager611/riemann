RECORD: Lambda <= 0.1383621 (was 0.13963), given the Lean-checked half-plane and HB's external inputs

# Re-optimised Hakimoglu-Brown barrier with the half-plane Re s > 0.87338196661320752 (node N2, 2026-11-05)

**Result.** Λ ≤ 1383621/10⁷ = 0.1383621 (< 0.13837). Rows: `table_theta.csv`. Check: `verify.py` (exact rationals plus Arb), all checks pass in 124 s with 37 MB peak memory on one core.
**Side result (needs only 7/8, our sorry-free QRH port):** Λ ≤ 695163/(5·10⁶) = 0.1390326 (`table_78.csv`, all pass). HB's own 24-row table left 0.0006 unused.
**Premises.**
1. ζ(s) ≠ 0 for Re s > θ = 5458637291332547/6250000000000000. This is Bebel's Lean-kernel-checked statement, tracker submission `2026-10-31-positive-gap-contact-return`, statement checked in R310. Its axiom list was not checked. The table needs only θ ≤ 0.873381966624.
2. HB's analysis in arXiv 2610.12803v1: Lemma head, Prop two-density, Lemma kernel-derivative and Theorem derivative-step. We audited it as sound in `../HBDBN`.
3. HB's external inputs, used as published and not replayed: the Polymath §8.4 claim (a) barrier, Platt's RH height 3.06·10¹⁰, Polymath Thm 1.5(i), and Props 3.1(ii), 6.1 and 6.2 with (37)–(41).

## 1. How θ enters (from the TeX, lines 98–147, 852–856 and 1187–1240)
- **One place only.** Re ρ ≤ θ plus ξ(s) = ξ(1−s) give |Im z| ≤ 2θ − 1 for the zeros of H₀, so M(0) ≤ q₀ = (2θ−1)². The paper says this itself (lines 132–137): the dynamics, the boundary input, the approximation band, the kernel inequality and all certificates are independent of θ.
- **In the table, q₀ touches only row 0:**
  - the Bernstein domain 0 ≤ r ≤ A₂ = q₀;
  - the own-pair terms 2λ/(B−q₀) and 2μ(B+q₀)/(B−q₀)² inside k₀;
  - the first gate E₂₀(k₀d₀)(2+k₀q₁) ≥ 2+k₀q₀.
- **The other table quantities do not depend on θ.** C_i, G_i, the nine scalar bounds (Z, ψ, ψ₂, δ₀–δ₂, ρ₀–ρ₂) and ε₀ depend only on the time window [t_i, t_{i+1}]. The constants p = 61/20, r₀ = 1/20, N* = 69098, N₁ = 65536, L/4 > 5.5716, X and b = 1/25 are fixed.
- **Why the fixed table gains nothing (R293).** In HB's fixed Table 2 a smaller q₀ only adds slack to gate 0. To gain, the whole chain must be recomputed: q₁…q_n, λ_i and μ_i (the Bernstein domains move with q_i), and the final step.
- **External inputs.**
  - Polymath's wall: H_t ≠ 0 on [X, X+1]×[1/5, 1]×[0, 1/5], with X = 6·10¹⁰ + 83951.5. It feeds Lemma head, which keeps zeros with |x| ≤ X below height 1/5.
  - Platt's RH verification: covers |x| ≤ X+1 at t = 0.
  - The approximation on 3 ≤ y ≤ 6: valid for 0 ≤ t ≤ 1/5 and x ≥ X.
  - None of these depends on θ or on the target Λ. They only impose t_n ≤ 1/5, final squared height b ≥ 1/25, and p − r₀ ≥ 3. We stay inside HB's range: p = 61/20, b = 1/25, t ≤ 0.1184, x ≥ X. The audit's caveat about the very end of Polymath's last step (t near 0.2) is not touched.

## 2. Re-optimisation (same theorem, new table)
| Change | HB 2610.12803 | here |
|---|---|---|
| q₀ | 9/16 | 5576563720/10¹⁰ ≥ (2θ−1)² = 0.55765637197 |
| time grid | 23 × 1/200, plus 463/10⁵ | 473 × 1/4000, plus 1121/10⁷ |
| C_i, G_i | Prop two-density, rounded to 10⁻⁴ | same Prop (t⁻ = t_i, t⁺ = t_{i+1}), from Arb upper bounds with denominator 10¹²; rounded to 10⁻⁸ |
| λ_i, μ_i | 3 decimals; every Bernstein coefficient > 13/10⁴; 8 U-pieces | maximise k_i by a linear program over the Bernstein coefficients; 16 U-pieces; coefficients > 0; 6 decimals |
| q_{i+1} | 6 decimals | least multiple of 10⁻¹⁰ that passes the exact E₂₀ gate |
| final time t_n | 11963/10⁵ | 1183621/10⁷ |

Gain decomposition (all with HB's analytic constants):
| Setting | Λ ≤ | Status |
|---|---|---|
| HB: q₀ = 9/16, 24 rows | 0.13963 | HB, audited |
| q₀ = 9/16, 477 rows (`table_78.csv`) | 0.1390326 | verified |
| new q₀, HB's grid of 24 rows, rows re-fitted (`table_theta_h200.csv`) | 0.1387728 | verified |
| new q₀, 474 rows (`table_theta.csv`) | **0.1383621** | verified |
| new q₀, continuum limit of this method (float model, steps → 0) | ≈ 0.13834 | model |

- θ accounts for 0.00067 of the gain: dΛ/dθ = 0.414 (cf. R293's 0.405). Re-gridding and re-fitting λ and μ account for 0.0006.
- R306's model 0.13897 covered θ alone. It is superseded.
- **The method is saturated in the kernel inequality** (float model `gen/sharp.py`). Replacing Lemma 10's bound K_h/h ≥ 2/(D+A₂+3r) by the exact minimum over admissible heights changes k by less than 10⁻⁵. The active constraints sit at r = A₂ (a second zero at the contact height), with D ≈ 1.6B and, in later rows, D → ∞. There the bound is exact.
- **Remaining levers.** All lie outside the paper: the wall height (b = 1/25 costs b/2 = 0.02), the band y ≥ 3 that fixes p, the analytic constants of Prop two-density, and θ itself.

## 3. Certificates (`verify.py`, independent of the generator)
- **Tools and self-tests.** It uses python-flint fmpq_mpoly for the Bernstein algebra; the generator uses its own dict polynomials and Fractions. Self-tests reproduce HB's row-0 scalar entries (strict upper bounds, within 10⁻¹¹) and HB's least Bernstein coefficient (row 21) exactly.
- **Table-independent scalars.** E̅₁₆₀(13929/625) < 113X/1420 and 69098² < 113X/1420, so L/4 > 5.5716 and N ≥ 69098. Also J(2) < 1, log N₁ < 12, and ∫_{N₁}^∞ log y · y^{−81/40} dy < 10⁻³.
- **Scalar bounds per row.** All 9 × 474 printed U-values are strict Arb upper bounds (128 bits). From them, max ε₀ = 0.017852 < 1/2, min(C^cert − C_i) = 2.8·10⁻¹¹ and min(G_i − G^cert) = 5.2·10⁻¹³.
- **Kernel positivity.** All 189,600 Bernstein coefficients are > 0; the least is 2.2·10⁻⁹. Hence Φ_i ≥ 0 on D ≥ 0, 0 ≤ r ≤ q_i.
- **Gates and chain.**
  - 8.73885 ≤ k_i ≤ 10.39438, all exact and positive.
  - Every gate holds exactly; the least margin is 5.8·10⁻¹³ (tiny, but these are exact rational inequalities).
  - The chain conditions hold: q₀ ≥ (2θ−1)², q₀ < B, q strictly decreasing, q_n = 1/25, t_n ≤ 1/5, and all parameters ≥ 0.
- **Conclusion.** Λ ≤ t_n + q_n/2 = 1383621/10⁷.
- **Commands.** Run from this directory:
  - `uv run --no-project --with python-flint==0.9.0 python verify.py table_theta.csv 5458637291332547/6250000000000000`
  - `… verify.py table_78.csv 7/8`
- **Regeneration.** `gen/build.py 4000 new|78 out.jsonl` then `gen/totable.py`. It needs numpy, scipy and python-flint, takes about 5 minutes, and runs single-threaded.

## 4. Novelty (checked 2026-11-05)
- **arXiv math.NT listings.** 2026-11 (140 entries) and 2026-10 (854 entries). The only de Bruijn–Newman titles:
  - 2610.12803 (HB, still v1, 0.13963);
  - 2610.18214 (Planat: 0.148612, under extra premises he states himself);
  - 2610.27343 (Dirichlet analogue);
  - 2611.06245 (Planat–Raymond, 4 Nov 2026: unconditional Λ < 0.156821, unrefereed). This is new for TREE N2, which lists 0.2 as the unconditional record.
- **Other sources.** math.CV 2026-11: none. arXiv search: none. Zenodo API (161 hits; newest checked): no Λ claim below 0.13963.
- **QuasiRiemannTracker.** It has no Λ entries. Two later half-planes exist: 0.87338196661320751993 (Watson–Bebel, 31 Oct) and "the closest declarable fraction above (777−√1081)/852" (Bebel, 1 Nov). Both move θ by less than 10⁻¹⁹, and the table already covers them.
- **Verdict.** No published Λ below 0.13963 was found, so 0.1383621 would be the record. The reviewer verdict is in §5; the orchestrator's own novelty pass is still needed.

## 5. Status and caveats
- This is a new table for HB's theorem; the analysis is unchanged. Our three new parts are the time grid, the per-row optimisation, and the half-plane input.
- **Reviewer (verdict only, 2026-11-05): SOUND.** No hypothesis of Theorem derivative-step or Prop two-density depends on the grid, and `verify.py` covers every table hypothesis. The step to M(0) ≤ (2θ−1)² is correct in the paper's normalisation.
- Premise status:
  - Bebel's half-plane statement: Lean-kernel-checked per the tracker; its axioms were not inspected.
  - Polymath's barrier computation: published, not replayed (as in HB).
  - 7/8: Lean-checked in our QRH port.
- RH relevance: RH ⟺ Λ ≤ 0. This shrinks the open interval [0, Λ⁺] from 0.13963 to 0.1383621.
