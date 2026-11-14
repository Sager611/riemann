**AUDIT PASSED**: x₁ < exp(3 484 971 332 733 178 124) holds as written. Every step of Bi–Hou's criterion re-derives, every constant recomputes (E = 2.39·10⁻¹¹, M = 22.767065, δ = 0.182865), and K2(v) = 1.000305250386512 recomputes from Gallego's separately certified zeros (rigorous ±6.9·10⁻⁷ against a margin of 3.05·10⁻⁴). Only the published inputs (Pintz's Lemmas 4–5 and Theorem A, and Platt–Trudgian) are taken on trust.

Audit of hal-05838582v1 (Zhenyu Bi and Enze Hou, 26 Oct 2026) and of github.com/bzyzz/mertens_bound @ a0b267c. Agent: bihouaudit (Orca), 2026-11-14, 18:32–19:10 UTC. No Lean, no commits. Line-by-line chain with source line numbers: `FULLPROOF.md`.

## 1. Short answer
- **Criterion (Thm 4.8)**: no false step.
  - Lemmas 4.1–4.3, the residue identity, the contour of Prop 4.4 and the bounds B1–B6 of Lemma 4.5 all re-derive.
  - Off-line zeros are handled by the ¬Q(u) case. It matches Pintz's Theorem A and Corollary A with Y = e^u exactly. We checked this against the 1987 scan (numdam).
  - Reviewer verdict, first pass: GAP. My paste omitted the hypothesis that T0 is not a zero ordinate.
  - That hypothesis is in the paper (l. 397, 508). We certify it: δ = 0.182865248898 (Arb counts, mpmath). Reviewer verdict, second pass: CORRECT.
- **Constants**: the six B terms, E, all hypotheses, the window, Remark 7.2 and three Table 2 rows recompute exactly (`e_bound.log`).
  - Prop 5.1 was redone over all 6 185 378 cells. It gives min N_ν = −8.767065 at ν = 5 042 995, the same cell as theirs, so M = 22.767065 ≤ 22.7671 (`prop51.log`, Arb cross-check in `arb_checks.log`).
- **K2**: summed from Gallego's published zeros, not from Bi–Hou's code. 2 325 572 zeros lie below T0.
  - The main part (n ≤ 6·10⁵) is 1.000305250386499, theirs 1.000305250386499012. The total is 1.000305250386512, inside their ball [1.0003052504 ± 2.3·10⁻¹¹].
  - Gallego's stated radii give a rigorous enclosure of ±6.9·10⁻⁷ (r_γ = 5·10⁻²⁶ makes the phases accurate to 1.7·10⁻⁷).
  - In fact the γ data are accurate to 6·10⁻³³ (mpmath), so the float value is good to ~10⁻¹⁴.
- **Falsification**: no slip found.
  - The hit is robust: K2 = 1.0003031 at v ± 10⁻⁶ and 1.00005 to 1.00013 at v ± 10⁻⁵.
  - The point v is a true local maximum: dK2/dv = 1.4·10⁻⁸, so the optimal shift is 3·10⁻¹⁵.
  - The slack in M is 16.36, i.e. M ≤ 39.1 would still work, so even a log t error in Pintz's 11 log t constant would be absorbed.

## 2. Criterion chain: findings (details in FULLPROOF.md A–B)
| Item | Finding | Effect |
|---|---|---|
| Statement used | Pintz 1987 §4–5, refined as in StR 2014 and re-derived for general (k, T0, T1). Pintz's L4, L5 and Thm A match the scan word for word. L7 is replaced by an Arb check, which we reproduce: 0.455603. | none |
| Gaussian identity (L4.3) | 1/(sζ(s)) = ∫M(e^w)e^{−ws}dw and (1/2πi)∫e^{ks²+as}ds = e^{−a²/4k}/(2√(πk)) | correct |
| Shift v = u + k | k(s²−¼) + u(s−½) = k(s−½)² + (u+k)(s−½) | correct; numerically irrelevant (γk ≤ 2·10⁻⁴) |
| Window 2√(ku) vs KN's √y | 2√(ku) comes from Pintz's M1. KN quote the weaker √y and note 2√(ky) after their (2.7), via StR Lemma E. | correct; 2√(ku) = 45 727.27 |
| Smoothing tail (L4.2) | The ½ of the erfc bound is dropped. | harmless (bound twice too big) |
| Off-line zeros | ¬Q(u) ⇒ β0 > ½ + 5logγ0/u. This is Cor A with Y = e^u, whose side condition Y > γ0⁵ holds automatically. The "+1" in (4.1) is StR Cor D with ε = 1. Under Q(u), the L5 distance satisfies σ − β > 1/u since uT1 > 64. | correct |
| Sign case h < −1 | \|K2\| is used throughout, and L4.1 is symmetric. | correct |
| One-sided strengthening | With a positive kernel and K2 > 0, case Q(u) gives M(x) > +√x in the window. Gallego states this; Bi–Hou §9.3 claim only \|M\|. | understatement, not an error |
| Range of M used | Only the trivial \|M(x)\| ≤ x (L4.2) and the window hypothesis. x₁ = min{x > 1 : \|M(x)\| ≥ √x}; ⌊x⌋ ≥ 2 is argued. | correct |
| R462-type slips | Every constant re-derived (c0 = tan 1, e^{150}, 14, 16 = 2+14, 18 = 4+14, 11, 2, 15/π², 1/π, 2^{−11}, Table 3). No misprint changes a result. | none |

## 3. K2 at v, recomputed (`k2_gallego.py` → `k2_bihou.json`, 61 s, 109 MB)
- **Data**: Gallego's per-zero files (Zenodo 23700212, `cc_*.crho`), 17 data members (n ≤ 9.14·10⁶) of the 833 MB zip taken with HTTP Range requests; all SHA-256 match his SHA256SUMS.
  - Each row gives γ to 42 digits (r_γ = 5·10⁻²⁶) and c_ρ = 1/(ρζ′(ρ)) to 30 digits.
  - Bi–Hou publish no per-zero data: their 330 MB cache is regenerable and not tracked.
- **Method**: phases γv mod 2π in exact integer arithmetic (v exact, 2π to 120 digits); float64 terms summed with fsum.
- **No truncation T\***: all zeros below T0 are summed. The part n > 6·10⁵ (γ > 3.77·10⁵) has Σ 2e^{−kγ²}\|c_ρ\| = 7.385·10⁻¹², computed from certified \|c_ρ\|. This replaces an analytic bound on Σα, for which no unconditional one exists.
- **Validation** (`val_mpmath.py`, 253 rows, 0 failures; mpmath at 50 digits, an independent library):
  - The rows are 232 uniformly random indices, n = 1..15 (the heaviest weights), n = 98 625/6, n = 600 000/1, and n = 2 325 571/2 (the highest used).
  - Each γ was refined by Newton on Z(t), then ζ′(ρ) was recomputed.
  - Worst \|Δγ\| = 6.0·10⁻³³ (1.2·10⁻⁷ of the radius). Worst \|Δc\| = 0.16 of the radius, 1.8·10⁻²⁴ relative.
  - γ_{N+1} = 1300000.18286524889819717690167 (mpmath; `delta_check.log`).
- **Results vs Bi–Hou**:

| Quantity | Bi–Hou | Ours |
|---|---|---|
| N(T0) | 2 325 572 | 2 325 572 data rows; Arb N(1.3·10⁶) = 2 325 572 |
| δ | 0.182865248898197 | 0.182865248898197 |
| K2, n ≤ 6·10⁵ | 1.000305250386499012059 | 1.000305250386499 |
| tail Σ 2e^{−kγ²}\|c\| | 7.3856·10⁻¹² | 7.3852·10⁻¹² |
| K2 | [1.0003052504 ± 2.30·10⁻¹¹] | 1.000305250386512 (rigorous ±6.9·10⁻⁷) |
| K2, n ≤ 98 625 (KN's zeros) | — | 1.0005544 (zeros 98 626…6·10⁵ add −2.49·10⁻⁴) |

- **Independence**: Bi–Hou's zeros, Gallego's zeros and our sum are three pipelines. Both data sets use Arb, but our mpmath sample does not.
  - Agreement to 1.3·10⁻¹⁴ rules out any zero-set discrepancy that contributes more than ~10⁻¹⁴ at v; the margin is 3·10⁻⁴.
  - The index also agrees at height 2.65·10⁶: ν = 5 042 995 in both data sets.

## 4. Falsification attempts (all negative)
- **Neighbourhood of v**: K2(v ∓ 10⁻⁶) = 1.0003032 / 1.0003031; K2(v ∓ 10⁻⁵) = 1.000133 / 1.000054.
  - The excess over 1 survives a 2·10⁻⁵ interval. Peak curvature K2″ ≈ −4·10⁶.
- **Gallego's y at Bi–Hou's (k, T0)**: K2 = 0.99976 < 1. Each paper's point works only with its own k.
- **Bi–Hou's 25 sub-record candidates** (`rescore.log`), k = 1.5·10⁻¹⁰, zeros n ≤ 742 986:
  - The four below 3.485·10¹⁸ reach \|K2\| ≤ 0.9902. These are 2.789·10¹⁸ (0.9902), 8.90·10¹⁷ (0.9842), 4.66·10¹⁷ (0.9839) and 1.507·10¹⁸ (0.9784).
  - Their CSV agrees to ≤ 10⁻³ (polish conventions differ). No better point is hidden in their list.
  - One quadratic extrapolation (y = 1.404·10¹⁹, dy = −4·10⁻⁴, outside its validity) suggests \|K2\| ≈ 1.008. Even if real, it would not beat 3.485·10¹⁸.
- **Zero count**: Arb N(t) = 2 325 571 at 1299999.53 and 2 325 572 from 1299999.54 through 1300000.1828. It is 2 325 573 at 1300000.1829, so T0 is not an ordinate.
- **Slack**: B3 dominates E; E < 3.05·10⁻⁴ holds for any M ≤ 39.1, i.e. 16.36 of slack in M.
  - The worst cluster up to 10⁷ is N_ν = −8.77 (theirs and StR's). Gallego's worst on [8·10⁶, 2.25·10⁷] is −10.25.
  - Replacing 11 log t by 12 log t costs ≈ log T0 = 14.1 < 16.36.
- **Pintz L4 constant**: \|log\|ζ(s)\| − Σ_{\|γ−t\|<c0} log\|s−ρ\|\| ≤ 4.3 at 24 points (t = 5·10⁵ to 6.5·10⁵, σ = 0.5001 to 2). The bound 11 log t + 14 ≈ 160 is never approached: the worst ratio is 0.026 (`arb_checks.log`).
- **Lemma 3.3**: min over 8192 Arb balls is 0.45560319 (theirs 0.455603…).

## 5. Gallego (Zenodo 23700212) vs Bi–Hou (full table: FULLPROOF.md D)
- **Same lattice point**: y − v = 1.44·10⁻⁵. Bi–Hou's repo has listed y_base = …396.5383311641 since 15 Oct. Gallego cites neither Bi–Hou nor the source of y ("a given number").
- **Different certification**, so the two works are independent checks of one candidate, not two discoveries:
  - k = 5.7·10⁻¹², T0 = 7 999 035.172, 16 622 627 zeros, S = 1.00207.
  - T1 = 22 494 898, with M = 24.25 taken from Platt's LMFDB zeros.
  - E < 3.9·10⁻⁴⁶; Bi–Hou's form of E at his parameters gives 1.08·10⁻⁴⁹.
- **Gallego's bound is slightly stronger**: B = 3484971332733141310.42 < …178 124 (W = 8 913.9 vs 45 727.3). Both round to exp(3.485·10¹⁸).
- **Differing conventions**:
  - x₁: Gallego uses \|M(n)\| > √n and Bi–Hou use \|M(x)\| ≥ √x; the bound holds for both.
  - Off-line zeros: StR Cor D with ε = 1 (\|M\| > e√x) vs Pintz Cor A.
  - Smoothing tail: 2e^{−u/4} vs 2e^{−u/2+2√(ku)}.
- **No conflicting constant.**

## 6. What stays unchecked (none load-bearing for the margin)
- **Literature inputs**:
  - Pintz L4/L5: the proofs are sketched in the 1987 paper. Our numerical sanity check leaves a 38× margin, and the B3 slack absorbs a log t error.
  - Pintz Thm A: the statement and proof skeleton match the scan.
  - Platt–Trudgian: RH and simplicity up to 3·10¹².
  - Every record since 1987 rests on the same inputs.
- **Completeness of Gallego's list**: certified by his own row-by-row Arb run plus the Arb count. We re-verified 253 rows plus counts at T0, not all 16.7 M rows (that would take 42 thread-hours).
  - Bi–Hou's own Arb zero list agrees in N(T0), δ, the ν index and K2 to 10⁻¹⁴.
- **Bi–Hou's Lean development** (`lean/`): not run. It covers only the final deduction from four stated hypotheses.

## 7. Files (this directory) and rerun (≈ 15 min, one core, peak 215 MB)
- `FULLPROOF.md`: the chain with bh.txt line numbers; Pintz scan cross-reference; numerics table; Gallego comparison.
- `k2_gallego.py` → `k2_bihou.json`, `k2_bihou_sample.json`, `k2_gallego.log`.
- `e_bound.py` → `e_bound.log`.
- `val_mpmath.py` → `val_mpmath.log`.
- `prop51.py` → `prop51.log`.
- `arb_checks.py` → `arb_checks.log`.
- `delta_check.py` → `delta_check.log`.
- `rescore.py` → `rescore.log`, `rescore.json`.
- **Rerun**:
  1. Fetch the zip's central directory: curl -r over the last 2 MB of `Gallego_2026_Mertens_bound_package.zip` (Zenodo 23700212).
  2. Range-fetch the members `package/crho_format.py` and `package/data/cc_00000001…cc_08056619.crho` (zip bytes 0–455 381 407), inflate them into DATA, and put crho_format.py next to the scripts.
  3. Run each script with `uv run --no-project --with numpy [--with mpmath | --with python-flint==0.9.0] python3 <script> DATA …`.
- **Sources**: bh.pdf from hal.science/hal-05838582/document; Pintz from numdam AST_1987__147-148__325_0.pdf; KN TeX from arXiv 2502.21021.

RH relevance: N10 (explicit Mertens bounds; the record x₁ < exp(3.485·10¹⁸) is sound and can be cited in TREE N10 and NOGO as audited; it is not a step toward RH, since the Mertens conjecture is false).
