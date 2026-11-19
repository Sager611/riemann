PADHI HOLDS WITH CORRECTED INPUTS: with Bellotti's corrected (or deleted) constants Padhi's certified max q is unchanged, 0.991368 on log x ∈ [872, 877], because his far term 2√(2π)·2.4·I′ is valid as printed (the 2.4 covers both sides) and Bellotti's tables never attain the minimum in his N⁺ (≥ e^{7.5} above his own cells); still unchecked and load-bearing: his two-zero repulsion table (Thm 5.2; q ≥ 1.17 without it, ≥ 1.07 with it halved).

# padhi2: Padhi's Ramanujan-inequality certificate re-audited with Bellotti's corrected zero-density constants

Claim audited: Padhi, "Ramanujan's prime-counting inequality holds for all x ≥ 38,358,837,683" (Zenodo 10.5281/zenodo.23440770 v2,
17 Oct 2026, still the latest version on 19 Nov 2026; code github.com/jaideepsaipadhi/ramanujans-inequality, last push 17 Oct 2026). Its own part is
e^58 ≤ x ≤ e^3159 (Thm 1.1); the rest is Johnston 2022 and Axler 2022 (cited; for Axler's lower end see §3).

## Short answer
- Padhi's far term is valid as printed. His C = 1.0067(1 + 2δΛ) + 2√(2π)·2.4·I′ (Lemma 6.2(v)) is Bellotti's c₄ with I replaced by
  I′ = δ n(ω) Σ_{i≥0} g(1 + iω). Per representative it asserts Σ_far ≤ √(2π)·2.4·n(ω)Σ_i g(1+iω)·X^{2δ}; the two-sided truth is
  ≤ √(2π)·2·n(ω)Σ_i g(1+iω)·X^{2δ} (windows of height ω tile each side, ≤ n(ω) each, g decreasing). 2.4 ≥ 2: 20% to spare. Our
  2026-10-22 objection doubled on top of the 2.4 (as with Bellotti's (4.13)), so its q = 1.03–1.07 near log x = 875 is void.
- Bellotti's tables never enter. In Padhi's certified scan N3 (Table 2 rows, the uniform 1.62·10¹¹ T^{1.448(1−σ)}) is never the
  minimum: wherever Padhi's own cells N4 exist, N3 ≥ e^{7.54}·N4; where they do not, the region is certified zero-free and skipped.
  Corrected constants (×1.78–2.05) or deleting N3 leave every q unchanged: max q = 0.991368 on [872, 877], as published.
- c₅ cites Bellotti's invalid |G(½+it)| ≤ d₅W/log W, but in all 3320 cells d₅·wL ≥ 11.47 > 9.757, so the certified
  9.75707 W/log²W is smaller and c₅ stands.
- Not touched by the Bellotti flaws, but load-bearing and still not re-derived by us: the repulsion table d* (Thm 5.2; §3).

## 1. The zero-density step exactly (main.tex of the record; TeX line numbers)
- N⁺ = min(N1 zero count (HSW), N2 Kadiri–Lumley–Ng Table 1, N3 Bellotti Thm 1.1 + Table 2 rows l = 29…41, N4 own Thm 6.1),
  l.601–606; "In the binding range only (N4) matters."
- N4 (Thm 6.1, l.618; Lemma 6.2, l.631): N(1−δ, T) ≤ ñ (d₄d₅/w)·C/(c₁² − c₅)·e^{2xδ_k l}/(δ_k l), evaluated at cell corners.
  C is in units X^{2δ}/(2δ) per representative: diagonal 1.0067, near 1.0067·2δΛ (Lemma 6.3, both sides explicit), far 2√(2π)·2.4·I′.
- Bellotti (2405.12545v2 TeX l.858–871): per r, Σ_{1<|γ_r−γ_s|≤T}|Γ(2−ρ_r−ρ̄_s)| ≤ √(2π)·2.4·I/δ, I = ∫₁^T g; c₄ = 2√(2π)·2.4·d₅d₄I/w
  (the 2 is |X^w − Y^w| ≤ 2X^{Re w}). bellotti2: the two-sided need is 2 + 4δe^{1/6−π/2}/I ≤ 2.12 ≤ 2.4.
- Lemma 6.4 (l.736): any window of height ω centred at t ≤ T holds ≤ n(ω) representatives (Hadamard budget, or ⌈ω/δ⌉ + 1).
  Replacing I/δ by n(ω)Σ_i g(1 + iω) needs no boundary term; windows above T shift to centre T − ω/2. One- or two-sided: the
  printed 2.4 carries the two sides, exactly as in Bellotti; I′ itself is one side.
- Code = paper (certify/arb/cells_arb.py l.258–263; Iprime_arb takes the worse of δ = (k−1)/10⁴, k/10⁴; tail y > 60 by e^{−1.4y}).
- Bellotti's other slips are not inherited: (4.11) replaced by Lemma 6.3 (two variables per rectangle index), (4.9) replaced by
  X^y/y increasing on [2ν, 2δ] (2νxL ≥ 1), (3.10) replaced by Padhi's own B₁, d₄ by his Thm 6.5. Only (4.15) (c₅) and the
  constants 1.0067 (≥ Γ(1+x) on [0, 1]), 2.4 and d₅ (Lemma 4.3) come from Bellotti.
- Paper/code mismatch, inert: the paper applies N3 for σ ≥ 0.9927; the code uses each row's own α₀ (0.9907–0.9927).

## 2. Padhi's certified pipeline rerun (code/qlab.py; his certify/arb/scan_arb.q_arb, his kernels (c, f) from scan_arb_{0,1}.txt;
any admissible kernel gives a valid bound; a labelled copy of his logN, checked equal to his on 3000 calls per run, records
which bound attains the minimum; ball arithmetic, 100 bits, as in his code)
| run | intervals | max q on [57, ·) | where | N3 attains the min |
|---|---|---|---|---|
| (a) published constants, (4.13) as printed | 126, [43, 1100) | 0.991367809 (= published, max diff 1.1e-16) | [872, 877] | 0 of 3.85·10⁶ calls |
| (b1) Table 2 C 'both flaws' (uniform C ×2.046) | 27, [797, 1000) | 0.991367809 (identical) | [872, 877] | 0 |
| (b2) Table 2 C 'all fixes' (uniform C ×2.046) | 27, [797, 1000) | 0.991367809 (identical) | [872, 877] | 0 |
| (b3) N3 deleted (every Bellotti bound) | 27, [797, 1000) | 0.991367809 (identical) | [872, 877] | deleted |
| (b4) N3 deleted, every 6th interval of [1100, 3159) | 35, [1100, 3150) | 0.759782448 (identical) | [1100, 1110] | deleted |
Every interval of (b1)–(b4) has the published q to the last digit. The rigorous reason (code/n3bind.py, reviewer CORRECT):
- N3 is eligible only for log T > log(3·10¹²) = 28.7297, and the scan evaluates N⁺ only at heights log t ∈ (log H, log T₁],
  max log T₁ = 36.2152 ([3120, 3130]); so every evaluated height has Padhi's cells N4 (28 < log T ≤ 40, δ ≤ 0.03).
- N3 increases in δ and log T; at each cell's lower-left corner, over all 418 (cell, N3 bound) pairs, min log(N3/N4) = 7.5365
  (cell 56,38, row e^37 < T ≤ e^38, a height never reached; 7.863 over reached heights). The 1833 pairs without an N4 cell lie in the certified zero-free region.
- Per bound (log T ≤ 36.2152): row 1 (T ≤ e^29, δ ≤ 0.0073) lies wholly inside the zero-free region (edge δ 0.007374 at
  log T = 29); rows 2–9 are evaluated only for δ between the zero-free edge and 1 − α₀ ≤ 0.0086, with min log(N3/N4) between 7.863 (row 5, cell 66,33) and 7.901 (row 2); rows 10–13 are never
  evaluated; the uniform bound is evaluated on (28.73, 36.22] but min log(N3/N4) = 20.29 there.
So any correction of Bellotti's Table 2 or Theorem 1.1, and their deletion, leaves Padhi's q, coverage and transfer unchanged:
[57, 3159) certified, max q 0.991368 on [872, 877], transfer ratio 0.990514 at log y = 877 (his transfer_arb.log, unchanged inputs).

## 3. Falsification
- Ranges of the Bellotti rows: see §2; none ever attains the minimum; the nearest miss at reached heights is a factor e^{7.86} (e^{7.54} over all cells). Rows 10–13 and every
  height above e^{36.2152} are never reached; row 1 sits inside the certified zero-free region.
- Transition to the other bounds (argmin counts, 126 intervals): N4 91.2%, N1 4.9%, N2 (KLN) 3.9%, N3 0. N4 covers δ ≤ 0.03,
  28 < log T ≤ 40; beyond δ = 0.03 only N1/N2 compete. Sensitivity on [867, 887): KLN A, B doubled → max q 0.991916 on [872, 877] (+0.00055); KLN deleted → q 852–1165 (KLN carries σ < 0.97; its Table 1
  rows were checked in RamRepair).
- Rounding (zero-density step): 1.0067 ≥ Γ(1+x) on [0, 1] (≤ 1); 2.4 ≥ 2 (two sides); √(2π), M(ε, L), −ζ′/ζ(1+ε) in arb;
  10.45 > 4ζ(3/2) = 10.44950; 5.922 ≥ I_Γ = 5.92164 and 4.72 ≥ 4.71363 (our first audit; his gammaint_cert.log); Iprime takes the
  worse δ of the cell and bounds the windows beyond y = 60 by e^{−1.4y} (valid: g(y) ≤ e^{(1/6−π/2)y} = e^{−1.404y}); n(ω) is the
  floor of a certified upper bound; tol uses the lower end of ε_R, q the upper end. Bellotti's own rounding (C up, 2 decimals) is
  irrelevant here (N3 inert). Found nothing in the unsafe direction.
- Load-bearing inputs not re-derived (sensitivity, his kernels): with the repulsion table switched off (NO_REP=1, d* = 0 in
  Lemma 6.3; 421 cells lk ≤ 32, k ≤ 140 recomputed by his cells_arb.cell, log N⁺ up by ≤ 0.4285 (cell 69, 31)), q exceeds 1 on
  [827, 947) at least: 1.017 [827, 837], 1.143 [872, 877], 1.1735 [917, 927] (logs/qv_norep_797.log; a lower bound, since cells
  with k > 140 or lk > 32 keep d*). With d* halved (DSTAR_SCALE=0.5, same cells; log N⁺ up by ≤ 0.1936) q is 1.027–1.074 on
  [847, 947), max 1.0742 on [877, 882] (logs/qv_dhalf_847.log): more than half of the claimed d* is needed. Kernel re-optimisation gained ≤ 1.7% in our first audit, so Padhi's Thm 1.1 stands or falls
  with Thm 5.2 (two-zero repulsion via BTY's Lemma 6), which no audit of ours has re-derived (RamRepair avoided it with the residue
  factor (1 + (Y/X)^{2ν}) instead).
- Thm 4.6 (λ*) and Thm 6.5 (d₄): RamRepair's independent re-derivations came out on Padhi's safe side (λ(29) 0.21494 vs his
  0.21384; k₁ ≤ −0.59807, k₀ ≤ 0.37066 vs his −0.5853, 0.4214), not fully certified.
- Outside the zero-density step, from RamRepair/NOVELTY.md (23 Oct, not re-checked here): Axler 2022 Thm 1.1 (x ≥ e^{3158.442}),
  which Padhi's Thm 1.2 needs above e^{3159}, uses Dudek–Platt Lemma 2.1 in its faulty case; with the repaired lemma Axler's
  parameters give about e^{3160.3}. So (e^{3159}, e^{3160.3}) is not covered as written; his q_arb with his top kernel (c = 70,
  f = 0.40) gives q ≤ 0.400012 on [3159, 3161] and [3161, 3163] (logs/ext3161.log), so extending his cover (and rerunning his
  transfer step there, not done here) closes it; Padhi's paper stops at 3159.

## 4. Verdict
- With Bellotti's corrected constants: Padhi's Thm 1.1 certificate is unchanged (max q 0.991368 on [872, 877], transfer 0.990514
  at log y = 877), because his own estimate N4 is valid as printed and Bellotti's bounds are never used. The question of the
  brief is settled: PADHI HOLDS WITH CORRECTED INPUTS. The same holds with every Bellotti bound deleted.
- As written: the paper cites two Bellotti statements that are not proved as printed (Table 2 / Thm 1.1 for N3; the |G(½+it)|
  bound in c₅). Both are inert or dominated (N3 never binds; c₅'s bound exceeds the certified 9.757 W/log²W in every cell), so
  they are presentation defects, not gaps in the numbers. Our 2026-10-22 refutation (far term ×2) is withdrawn.
- What the claim still rests on, unchecked by us: Thm 5.2 (repulsion d*), which is load-bearing (q ≥ 1.17 without it, ≥ 1.07 with half of it), the
  certified d*, λ* and d₄ computations (his code, reproduced only where rerun here), and, for N_R = 38,358,837,683 itself, the
  range (e^{3159}, e^{3160.3}) where Axler 2022 needs the faulty case of Dudek–Platt Lemma 2.1. The claim N_R = 38,358,837,683
  is therefore not yet established by Padhi's paper alone; liqsweep (14 Oct 2026) claims it by another route.
- Reviewer (reviewer subagent, two passes, logs/reviewer1.md, reviewer2.md): pass 1 GAP (heights below log T = 28
  not addressed); closed (N3 needs log T > 28.7297; the scan evaluates only log t ∈ (log 3·10¹², 36.2152]); pass 2 CORRECT.

## Draft replacement for the Padhi line of REFUTATIONS.md (do not paste blindly: the entry is no longer a refutation)
- Padhi, "Ramanujan's prime-counting inequality holds for all x ≥ 38,358,837,683" (Zenodo 10.5281/zenodo.23440770 v2, 2026-10-17;
  unconditional claim): refutation withdrawn (padhi2, 14 Nov 2026). His far term 2√(2π)·2.4·I′ is valid as printed (the 2.4 covers
  both sides, as in Bellotti's (4.13); our 'factor 2', q = 1.03–1.07, double-counted), and Bellotti's Table 2 and Thm 1.1, whose
  C are understated ×1.8–2.0, never attain the minimum in his N⁺ (≥ e^{7.5} above his own cells), so his certified max q = 0.991368
  on log x ∈ [872, 877] stands with corrected constants and with them deleted. Harmless: c₅ cites Bellotti's invalid |G(½+it)| proof
  (the certified 9.757 W/log²W is smaller in every cell). Open: his two-zero repulsion table (Thm 5.2) is load-bearing (q ≥ 1.17
  without it, ≥ 1.07 with it halved) and unchecked, and (e^{3159}, e^{3160.3}) needs Axler 2022, which uses Dudek–Platt Lemma 2.1 in its faulty case.
  liqsweep (14 Oct 2026) claims the result by another route; our repair reaches q ≤ 0.928 without the repulsion table but was not
  completed ([RamRepair](Riemann/Deep/RamRepair/REPORT.md)). [re-audit](Riemann/Deep/Padhi2/REPORT.md) [first audit](Riemann/Deep/Ramanujan/REPORT.md)
- Consequence for Riemann/Deep/Ramanujan/REPORT.md: its first line (REFUTATION, far term ×2, q > 1) is wrong; its §5 repair
  numbers (q 0.8957) doubled the far term too and are conservative. RamRepair is unaffected (it doubled, so it has more slack).

## Files and runs (one process at a time, single-threaded; every run under .scratch/memcap.sh 290 and /usr/bin/time -l)
- code/common.py (imports Padhi's certify/arb/scan_arb.py from PADHI_REPO; variants published | both | all | none | kln2 | nokln;
  CELLS_OVERRIDE swaps in recomputed cells), code/qlab.py (rerun with argmin labels), code/qvar.py (plain rerun), code/n3bind.py
  (cellwise N3 vs N4, zero-free check, max log T₁, per-row ranges), code/cells_g.py (d₅·wL), code/dstar0.py (cells with d* = 0),
  code/ext3161.py (cover past 3159), code/runhalf.sh (d* halved), code/run.sh (watchdog wrapper), code/runvars*.sh (queues). NOTE.md: working notes.
- Follow-up (Thm 5.2): code/repcert_cells.py, code/m2corner.py, code/runcert.sh; logs/repcert_*.log, m2corner*.log, m2grid.log;
  NOTE_THM52.md. Inputs re-downloaded to .scratch/padhi2 (Zenodo source, GitHub tarball, arXiv 2603.21490 e-print), deleted after.
- logs/: qlab_pub_lo, qv_{both,all,none,norep}_797, qv_{kln2,nokln}_867, qv_none_hi, qv_dhalf_847, ext3161, n3bind, cells_g, dstar0_lo
  (+ norep_lo.json), dstar_half (+ dhalf_lo.json), q_pub_872, reviewer_brief{1,2}.md, reviewer1.md, reviewer2.md; .err files hold the memcap and /usr/bin/time peaks (max RSS ≤ 102 MB; memcap peaks ≤ 70 MB).
- Inputs (downloaded to .scratch/padhi2, deleted after the run): Zenodo 23440770 v2 (main.pdf, source tar with main.tex),
  github.com/jaideepsaipadhi/ramanujans-inequality main (codeload tarball, last push 2026-10-17), arXiv 2405.12545v2 e-print.
  To rerun: unpack the repo tarball to .scratch/padhi2/repo/ and run bash code/run.sh 1800 <name> code/qlab.py published 0 1100.

## Thm 5.2 audit (follow-up)

THM 5.2 HOLDS: Padhi's two-zero repulsion (Thm 5.2, the d* table) is a valid extension of BTY's Lemma 6, his code computes exactly the printed margin M₂, every d* entry used is backed by a ball-arithmetic certificate in his run, and our reruns and an independent evaluation found no failure; the entries are tight (the margin vanishes 3–5% above d* in the binding cells), and the inputs Thm 5.2 shares with his zero-free region (Thm 4.6) were not re-derived.

Source: BTY = Bellotti–Trudgian–Yang, "Zero-free regions inspired by work of Heath-Brown", arXiv 2603.21490 v1 (main.tex l.1182;
TeX e-print). Working notes: NOTE_THM52.md. Code: code/repcert_cells.py (his certify), code/m2corner.py (independent float M₂).

### 1. Statement vs source
- Thm 5.2 (l.580): for integers 29 ≤ l_k ≤ 37 and λ*(l_k) < k·10⁻⁴ l_k ≤ 0.0125(l_k − 1), two zeros with H₀ < γ_j ≤ e^{l_k} and
  1 − β_j ≤ k·10⁻⁴ satisfy |γ₀ − γ₁| ≥ d*(l_k, k); d* = 0 outside the range.
- What it takes from BTY: the identity of Lemma 6 (TeX l.507: Re Σ_ρ F(s₁ − ρ) ≥ F(σ − β₀) + F(σ − 1 + β₀) − 10⁻⁸ for one zero,
  via pairs {ρ, 1 − ρ̄} and Lemma 5's pair non-negativity), Lemma 5(ii) (tails, constant 44), Lemma 13 (poles, 10⁻¹⁰η), the
  definitions of w (Def. 1), κ_m ((18)), P (§4), K = 16, T₀ = 10¹⁰. All parsed from the TeX and equal to his values
  (Σc_k² = 14912370 so a₀ = 1, a₁ = 865534/497079, a = 2919857/828465, κ = 433/859, w₀ = 5.672787599).
- Hypotheses: BTY's Lemma 6 assumes η ≤ A₀/log t (A₀ = 1/4.8596, η₀ = 0.0071093, σ₀ = 0.99347); Padhi needs η ≤ 0.0145,
  λ₀ ≤ 0.35, σ ≥ 0.993. He re-proves what depends on this: the pair non-negativity (his Lemma 4.2, l.442; ε₀ =
  51.9·0.0145²/(5.67279·0.986²) = 0.0019786 < 1/500, the root-isolation step not rerun) and the tails (l.500; we recomputed
  44·η(log t + 1)·7.6·10⁻¹¹ ≤ 2.3·10⁻⁹ per k and 88·0.015·log T₀/T₀ = 3.0·10⁻⁹ for k = 0, total ≤ 6·10⁻⁸ < 2·10⁻⁷). Normalisation:
  F(z) = Σ κ_m W((z + (2σ−1)m)/η) is the Laplace transform of f (checked); the second pair uses the primary's η, as printed.

### 2. Proof, line by line (main.tex)
- l.555 + Lemma 5.1 (l.558–571): BTY's k = 1 zero sum, all pairs {ρ, 1 − ρ̄} with |γ| ≤ Kt + T₀ are ≥ 0 (pair non-negativity,
  0 ≤ σ − β ≤ 2σ − 1 by the zero-free hypothesis and symmetry); keeping a second pair replaces a dropped non-negative term by its
  value. The pairs of ρ₀ and ρ₁ are distinct (β₁ > ½ > 1 − β₀); a double zero counts twice; only |d| matters (F(z̄) = conj F(z)).
  s₁ − ρ₁ = (σ − 1 + η₁) − id and s₁ − (1 − ρ̄₁) = (σ − η₁) − id: the displayed term is exact. Valid.
- l.573–577: M₂ = M₁ + a₁ Re(F(σ − 1 + η₁ − id) + F(σ − η₁ − id)) > 0 for one ordering excludes the pair; F uses the primary's η.
  rep_cert.margin computes exactly this (r₁ + r₂ − lhs), and our independent float M₂ equals his ball value at every point tried.
- l.585: e_min = λ*(l_k)/l_k ≤ 1 − β (Thm 4.6, λ*/l non-increasing: true on all 61 grid values; code uses the grid value at or
  above l_k, smaller, safe). Lemma 5.1's σ − β₁ > 0 needs e_min > 1 − σ: 0.007374 > 0.00676 (l_k = 29), 0.005707 > 0.005425 (37).
- l.586–589: boxes [l − 1, l + 10⁻¹²] × [e_min, e_max]² × [0, d*], both orderings; d ≤ 0.15 moves log t by ≤ 5·10⁻¹⁴ < 10⁻¹².
  σ is a ball with A = λ*(log(Ke^{l₁} + T₀)), the minimum over [H₀, Kt + T₀] because λ* is non-increasing (checked): valid.
  Nit: in the swapped ordering the primary is the upper zero, whose σ can exceed the ball by ~10⁻¹⁵ at an interior box edge
  (effect on M₂ ~10⁻¹³; margins at the binding corners ≥ 1.0·10⁻³).
- Constants at η₀ = ⌈10⁴e_max⌉/10⁴ (up) and σ₀ = ⌊10⁵σ_lo⌋/10⁵ (down), valid for η ≤ η₀, σ ≥ σ₀. Rounding: E₂, E₃ × (1 + 10⁻⁶),
  c₁₀, G ÷ (1 + 10⁻⁶): each lowers M₂ (safe). These four are mpmath 20-digit values, not ball enclosures as l.416 says (minor;
  the factor moves M₂ by ~10⁻⁸). Values agree with the paper: c₁₀ = 0.0598 (l.843), G(σ₀) between 1.56839 and 1.56888 (l.497),
  E₂λ₀ ≈ 186 ≈ 185.4 at BTY's point (l.858).
- W (l.502–516): closed form = Laplace transform of BTY's w (400-point Gauss–Legendre agrees to 4·10⁻¹⁴); the second-order
  enclosure bounds sup|W''| by P(ν) or (2w₀ + Q(ν))/|z|³ with ν rounded down to a half-integer and capped at 60 (P, Q decrease
  in ν: safe); the radius is inflated by 10⁻¹² (safe).
- l.590–592: d* certified at 0.97 × the float value, lowered by ×0.9 until the bisection closes (repcert_run.py); the other k take
  the value of the next certified k above (box grows with e_max: valid). l.593–594: min over l ≤ l_k: all 507 table entries are
  ≤ that minimum, and each l's first certified k is ⌊10⁴λ*(l)/l⌋ + 1 (cells below are zero-free).

### 3. His certification recomputed (rep_cert.certify unchanged, python-flint balls, 80 bits, maxdepth 40, the stored d*)
| cell (l_k, k) | d* | binding? | result | boxes | time |
|---|---|---|---|---|---|
| (29, 76) | 0.061013 | yes (zero-free edge, l_k = 29) | certified | 2979 | 29 s |
| (32, 68) | 0.055368 | yes (log T ∈ (31, 32]) | certified | 5927 | 65 s |
| (31, 71) | 0.056284 | yes (l_k = 31) | certified | 7833 | 110 s |
| (30, 72) | 0.092150 | yes | not finished in 420 s (killed by our time cap) | – | – |
| (31, 69) | 0.129986 | yes (largest NO_REP effect) | not finished in 10 min (killed) | – | – |
(logs/repcert_*.log; memcap peaks ≤ 33 MB.) The two unfinished cells have the largest d* (wide d-range, margin ~10⁻³ at the
corner) and need far more boxes; his run certified them (repcert_run.py, Pool(2)); we did not reproduce those two.
- Binding cells near log x = 872–877 (T₁ ≈ e^{31.1}) use d*(l_k, k) for l_k = ⌈lk⌉ ∈ {29, 30, 31, 32} and k from the zero-free
  edge up to ~140 (k ≤ 120/121 have d* > 0). Every such entry equals the certified value of his run (dcert_29…37.json, 14/13/14/14
  directly certified k for l_k = 29/30/31/32) at the same or the next larger certified k, and is ≤ the minimum over lower l
  (checked for all 507 entries). So every d* entry used is certified by his run; we re-ran 3 directly certified binding cells to completion (all certified), 2 more did not finish within our caps.

### 4. Falsification (independent float M₂, logs/m2corner.log, m2corner2.log, m2grid.log)
The certificate needs, on each box, one ordering with M₂ > 0. On a 5 × 7 × 7 × 11 grid of every box below (log t, η₀, η₁, d),
the minimum over the grid of max(M₂, M₂ swapped) is at the top corner log t = l_k, η₀ = η₁ = e_max, d = d*, and is positive:
| cell (l_k, k) | d* | min over grid (at top corner) | his ball margin there | d where M₂ = 0 at that corner |
|---|---|---|---|---|
| (29, 74) | 0.142299 | 1.069·10⁻³ | 0.0010691 | 0.146756 = 1.031 d* |
| (29, 76) | 0.061013 | 1.247·10⁻² | 0.012472 | 0.063202 = 1.036 d* |
| (30, 72) | 0.092150 | 5.605·10⁻³ | 0.0056052 | 0.095889 = 1.041 d* |
| (31, 69) | 0.129986 | 1.370·10⁻³ | 0.0013702 | 0.134493 = 1.035 d* |
| (31, 71) | 0.056284 | 1.496·10⁻² | 0.014959 | 0.058354 = 1.037 d* |
| (32, 68) | 0.055368 | 1.161·10⁻² | 0.011608 | 0.058011 = 1.048 d* |
| (29, 120) extreme δ | 0.058174 | 0.4208 (corner) | 0.42083 | 0.067812 = 1.166 d* |
| (37, 120) extreme δ and height | 0.044945 | 0.6461 (corner) | 0.64613 | 0.056708 = 1.262 d* |
One ordering alone can fail earlier (primary η₀ = e_max, second η₁ = e_min: M₂ = 0 at 1.010 d* in (29, 76)), but there the
swapped ordering has M₁ > 0 (a zero with 1 − β = e_min at height e^{l_k} is excluded on its own), so the box still closes.
- Smallest configuration the method cannot exclude: two zeros at height e^{29} with 1 − β₀ = 1 − β₁ = 0.0074 (cell (29, 74), the
  largest d*) at distance 0.146756 = 1.031 d*; in the binding cells the same happens at 1.035–1.048 d*. Below d* every grid
  point and every corner is excluded, so no entry is contradicted; the entries are tight (1–5% from the method's limit), as
  Padhi's 0.97 × float recipe intends. (37, 121) has d* = 0 in the table (set by hand); (37, 120) is the extreme non-zero entry.
- Mechanism: at η₀ = e_max the single-zero margin M₁ is negative (−0.0346 at (29, 76), −0.707 at (29, 120)); the second pair's
  term (0.047, 1.13) restores M₂ > 0 and decreases with d, which is why d* is finite and the entries are tight.

### 5. Verdict and residual inputs
- Thm 5.2 HOLDS: the statement follows from BTY's Lemma 6 argument with the second pair kept (valid), its hypotheses outside
  BTY's range are re-proved (pair non-negativity for σ ≥ 0.993, η ≤ 0.0145; tails), his code computes the printed M₂ with a
  correct W, every d* entry is backed by a certified cell of his run (table logic verified on all 507 entries), and our reruns
  (3 binding cells) and the independent grid found no failure.
- Not re-derived here (shared with his Thm 4.6, so already load-bearing for the zero-free region): the λ* table, the constants
  E₂, E₃ (his Lemma 4.3), c₁₀ (Lemma 4.4), G (Lemma 4.5), BTY's explicit-formula identity, Lemma 5(ii) (44η/|Im z|²) and
  Lemma 13. RamRepair's independent BTY-type iteration (23 Oct) is consistent with his λ* (0.21494 vs 0.21384 at l = 29).
  Minor: those four constants are mpmath values with a 10⁻⁶ outward factor, not ball enclosures (paper l.416 says ball
  arithmetic); the σ-ball nit (§2). Neither can move a margin of ≥ 10⁻³.
- Consequence: the open point "repulsion table unchecked" in the Padhi line is closed; Padhi's Thm 1.1 (e^{58} ≤ x ≤ e^{3159})
  now rests on inputs that were either checked here or independently re-derived (not certified) in RamRepair. The full claim
  N_R = 38,358,837,683 still needs (e^{3159}, e^{3160.3}).

### Updated draft of the REFUTATIONS Padhi line (replace its "Open: …" sentence)
- Old: "Open: his two-zero repulsion table (Thm 5.2) is load-bearing (q ≥ 1.17 without it, ≥ 1.07 with it halved) and unchecked,
  and (e^{3159}, e^{3160.3}) needs Axler 2022, …"
- New: "His two-zero repulsion table (Thm 5.2), load-bearing (q ≥ 1.17 without it), holds (padhi2 follow-up, 19 Nov 2026: valid
  extension of BTY's Lemma 6, his margin reproduced independently, 3 binding cells re-certified, entries 3–5% below the method's
  limit; the λ* table and four constants it shares with his zero-free region were not re-derived, consistent with RamRepair).
  Open: (e^{3159}, e^{3160.3}) needs Axler 2022, which uses Dudek–Platt Lemma 2.1 in its faulty case (his own q_arb gives
  q ≤ 0.400012 on [3159, 3163]; transfer step not rerun)." The rest of the line stays.

RH relevance: N10
