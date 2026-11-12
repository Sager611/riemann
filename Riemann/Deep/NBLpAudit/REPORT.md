AUDIT PASSED — every step of PROOF.md §3–§7 (Lemmas T, M, PW, Γ, H, A, P, Theorem B, Theorem A (i)⟹(ii)) has a complete proof in FULLPROOF.md and none fails; the endpoint criterion is true, but it is new only at the single exponent p = 1/Θ.

Object: discovery `nyman-beurling-lp` (Riemann/Deep/NBLp/PROOF.md). Referee: Orca agent `nblpaudit`, 12 Nov 2026,
3 h budget. No Lean. Line-by-line proofs: [FULLPROOF.md](FULLPROOF.md) (§A–§F proofs, §G findings, §H numerics,
§I sources, §J significance).

## (1) Line-by-line verdicts

| Step | Verdict | Note |
|---|---|---|
| Prop. C (S_p, codim ∞, D_m-invariance, χ ∈ W_p ⟺ W_p = S_p) | ✓ | literal density in L^p(0,1) is false for every p, as P says |
| Prop. N (necessity, all 1 < p < ∞) | ✓ | `x^{w_0−1} ∈ L^q ⟺ β_0 > 1/p` |
| Lemma T, Lemma M (a < 1/p), Lemma PW | ✓ | PW: added `h ∈ L¹ ⟹ liminf h = 0` |
| Lemma Γ (a), (b) | ✓ | constants crude, true maxima −0.92, −0.42 vs 4, 5 |
| Lemma H (Hadamard midline) | ✓ | sharp: equality on `σ = θ* − ε/2` when all `β = θ*` |
| Lemma A (Lindelöf in zero-free half-plane) | ✓ | exponent and constant `C_A(θ,δ,η')` as stated; ineffective (`C_6`); Ω star-shaped (added) |
| Lemma P (Riesz–Perron) | ✓ | rate `N^{−η}` for `η < σ − θ`; sharp numerically |
| Theorem B (§6) | ✓ | open range, power rate |
| Theorem A (§7): A0–A3 | ✓ | shift to `σ' = 1/p − ε/2` legitimate (Ψ_ε analytic on `1/p − ε < σ < 1`) |
| A4 residue bookkeeping, signs, `B(ε) ≤ K_pε` | ✓ | poles only at `1` and `2σ'`; `J` formula confirmed to 2e−9; reviewer CORRECT |
| A5 (`O_p(ε^{1−1/p})`), A6 (finite span), A7 | ✓ | all constants depend on `(p, ε)` only; `N → ∞` at fixed ε |
| Use of the line `Re s = 1/p` | none | only `β ≤ 1/p` (Lemma H) and zero-freeness on `Re > 1/p` (Lemmas A, P with `δ = ε/8`) |

Mellin–Plancherel is used only in L² (weighted `L²_σ`); no Hausdorff–Young anywhere. Every L^p(0,1) norm comes from
an `L²_a` norm via Lemma M at `a = σ_0, σ', σ''`, each `< 1/p`. At `a = 1/p` the weight is only weak-`L^r`, which is
why the proof never sits on the line `1/p`. Gaps of exposition, all closed in FULLPROOF §G: PW's choice of `T_n`;
Ω simply connected; `X ≥ 1` in three circles; typo `S` for `S_p`; `e²` vs `e^{5/4}`. None changes a statement.

## (2) Cited results

Beurling 1955 as quoted verbatim in Dimitrov–Oliveira Thm A ✓ (PNAS scan not fetched). Báez-Duarte 2003
(math/0202141): Thm 1.1 = [C10] after the routine L²(0,∞) ⟷ L²(0,1) equivalence (FULLPROOF §G.9, valid for all
p > 1); its §"Some comments" already asks for the `L_p` extension. Báez-Duarte 2005 (math/0505453) Remark: L_p
sufficiency for the strong criterion not done ✓. Dimitrov–Oliveira (1608.07887) Thms 1–2: twisted-κ L²(0,∞) /
weighted L² on `Re s = 1/p`, endpoint included, not L^p with `{1/(kx)}` ✓. Titchmarsh (2.1.5), Thm 4.11; Davenport
Ch. 12; B–C and three circles: as used ✓. One citation slip: Binet's formula is not DLMF 5.9.1 (that equation is
`Γ(ν/μ)/(μz^{ν/μ}) = ∫exp(−zt^μ)t^{ν−1}dt`); the formula is correct (DLMF §5.9(ii), Whittaker–Watson §12.31).
[C9] and the 0.87338… half-plane: Lean-checked elsewhere (QRH port, hpaudit 11 Nov; TREE.md N7).

## (3) Numerical falsification (FULLPROOF §H; peaks ≤ 25 MB, one process at a time)

* A4 identity `J = 2π[R(2σ')/(2σ') − 1/((2σ'−1)ζ(1+ε))]`: quadrature vs formula agree to `2.6e−9, 2.1e−9, 9.2e−10`
  at `(p,ε) = (1.3,0.2), (1.1,0.1), (8/7,0.05)`. A coarse first run showed `1.3e−6`; it was the pole `s = 1` at
  distance 0.14 from the line spoiling the quadrature near `t = 0`, gone with fine panels.
* `|R|/|Γr_ε| ≤ 0.999` on the lines `σ'`; `2πB ≤ π/σ' − 2J + I_c ≤ bound` holds numerically; `B/ε² = 2.6–4.4` with the
  true zeros (the proof's `O(ε)` is the worst case).
* Lemma A: `max|log|ζ(σ+it)||/log(2+t)` falls with height (`σ = 0.75`: 0.62, 0.30, 0.27 on the decades up to 6000).
* Lemma P: `|P_N − 1/ζ|` decays exactly like `N^{−(σ−1/2)}` at `t = 14.13` (`N ≤ 10⁵`).
* `‖1 + Σ_{k≤N}μ(k)(1−k/N)ρ_{1/k}‖_p`, `N ≤ 10⁴`, `10⁷` breakpoints: local exponents for `N ≥ 300` are `0.45–0.53`
  (`p = 1.1`, Theorem B under RH: any `η < 0.409`) and `0.28–0.31` (`p = 1.3`, `η < 0.269`); `p = 2` stalls.
  Matches run1e9.txt inside its rigorous intervals. No slower-than-predicted decay.

## (4) What survives, and what the endpoint adds

Everything survives. Caveat on significance (FULLPROOF §J): Prop. N and the open range `p < 1/Θ` are routine
(Theorem B, or Dimitrov–Oliveira on a line `θ' ∈ (Θ, 1/p)` plus Lemma M). The endpoint is needed only for
`p = 1/Θ`, which lies in `(1,2)` only if RH fails. Hence Corollary 3 (RH ⟺ criterion for all `p < 2`) and the
discovery's unconditional ranges `p ≤ 8/7` and `p ≤ 1.14497` already follow from the routine parts once the
Lean-checked half-plane `Re s > 0.87338196661320752` is used (`1/0.87338… = 1.1449744`). The endpoint adds the exact
equivalence `χ ∈ W_p ⟺ p ≤ 1/Θ` and, unconditionally, the single exponent `p = 1.1449744…`. Suggest the
DISCOVERIES line say "exact endpoint" rather than imply the unconditional ranges need it.

Side checks (FULLPROOF §K): the Dirichlet extension's changed steps (Lemmas Γ_χ, E, H_χ; A4 for non-principal χ
with complex `R(2σ')`, `|1 − R(2σ')| ≤ ε|ζ'(u_0)|ζ(u_0)`) are correct; arXiv math.NT/CA/FA 2026-11 titles to 12 Nov
show nothing on L^p Nyman–Beurling criteria.

Method note: the scheme (fixed-ε family `Σμ(k)k^{−ε}ρ_k`, `n → ∞` first) is Báez-Duarte's 2003 one; his p = 2
proof bounds `ζ(½−ε+iτ)/ζ(½+ε+iτ)` by the functional equation. For `p < 2` the reflected lines about `1/p` are not
related by the functional equation; Lemma H (Hadamard product under (i)) replaces it, and the exact residue
evaluation of `J` supplies the rate `B(ε) = O(ε)` needed to beat Hölder's `ε^{−1/r}`. That is the new idea.

Files: FULLPROOF.md, reviewer_A4.txt, numA.py/.txt, numB.py, numB_T4000.txt, numB_T1000_fine.txt,
numB_T4000_fine.txt, numC.py/.txt, lpnorm.c, lpnorm_M1e7.txt.

RH relevance: N5
