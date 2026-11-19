# padhi2 NOTE (working notes; REPORT.md is the deliverable)

## 1. Padhi's zero-density step as printed (main.tex of Zenodo 23440770 v2; line numbers of src/main.tex)
- N+ = min(N1 zero count, N2 KLN Table 1, N3 Bellotti Thm 1.1 uniform + Table 2 rows l = 29..41, N4 own Thm 6.1) (l.601-606).
- N4 = Thm density (l.618): N(1-δ,T) <= ñ (d4 d5/w) C/(c1^2-c5) e^{2xδ_k l}/(δ_k l); Lemma D (l.631) (v):
  C = 1.0067(1 + 2δ_k Λ(κ0 δ_k)) + 2√(2π)·2.4·I'. Lemma crowd (l.736): I' = δ n(ω) Σ_{i>=0} g_α(1+iω), ω = 1/4,
  n(ω) = max reps in ANY window of height ω; replaces Bellotti's I = ∫_1^T g_α in c4 [(4.13)-(4.14)].
- Units: C is in units X^{2δ}/(2δ) per representative r (diagonal 1.0067 X^{2δ}/(2δ), Lemma pair l.708).
  Bellotti l.862-871: per r, Σ_{1<|γr-γs|<=T} |Γ(2-ρr-ρ̄s)| <= √(2π)·2.4·I/δ; c4 = 2√(2π)·2.4·d5 d4 I/w (the 2 = |X^w-Y^w| <= 2X^Re w).
  So Padhi's far term = Bellotti's c4 with I -> I'; per r it asserts Σ_far <= √(2π)·2.4·n(ω)Σ_i g(1+iω)·X^{2δ}.
- Needed (two-sided): each side's distances (1, T] are tiled by windows [1+iω, 1+(i+1)ω], <= n(ω) reps each, g decreasing
  on [1,∞) (Padhi's proof) -> Σ_far <= √(2π)·2·n(ω)Σ_i g(1+iω)·X^{2δ}. Printed 2.4 >= 2: VALID, 20% slack, no boundary term.
  (Windows centred above T: shift to centre T-ω/2, same reps; n(ω) needs only t <= T.)
- Bellotti's own (4.13): 2.4 >= 2 + 4δe^{1/6-π/2}/I covers both sides (bellotti2, reviewer CORRECT). Our 2026-10-22 audit
  doubled on top of 2.4 (4√(2π)·2.4·I'): that 'correction' was wrong; q = 1.03-1.07 is void.
- Code = paper: certify/arb/cells_arb.py l.258-263: coef = 1.0067(1+2δΛ) + 2√(2π)·2.4·Iprime_arb; Iprime_arb takes the worse
  of dn = k-1, k for n(ω) and the δ factor, g at δ_k (g increasing in δ for y >= 1): conservative.
- Near term: Lemma pair counts both sides explicitly ('two variables for each j'), so Bellotti's (4.11) one-sidedness (c3) is
  not inherited; X^y/y monotonicity handled by 2νxL >= 1 (Bellotti's (4.9) slip not inherited). B1 is Padhi's own (no (3.10) slip).
- c5 uses |G(1/2+it)| <= d5 W/log W (Bellotti (4.15), proof invalid). Certified 9.75707 W/log^2 W (BellottiFix) is smaller iff
  wL >= 9.75707/d5 (< 9.698): check min wL over all cells (section 2).

## 2. Does N3 (Bellotti's tables) ever bind? (code/n3bind.py, logs/n3bind.log)
- scan_arb.logN: N+ = upper_min(N1, N2 KLN (largest σ0 <= s, T > 3.061e10), N3u = 1.62e11 T^{1.448(1-s)} and N3e29 =
  37341.72 T^{14.160(1-s)} (s >= 0.9927, T > 3e12), N3r = Table 2 rows 1-13 (own α0 0.9907-0.9927, log T in (l-1, l], l <= 41),
  N4 = certified cell (k = ceil(1e4 δ), lk)). Blocks with s >= edge(log t) (zero-free) are skipped.
- Cellwise bound (rigorous for printed constants): N3 increases in δ and log T, so on cell (k, lk) N3 >= its lower-left corner.
  418 (cell, N3 candidate) pairs where both apply: min (N3_lower - N4_upper) = 7.5365 in log N (cell 56,38, row 10): N3 is
  >= e^7.5 ≈ 1870 × N4 wherever both exist. Corrected C (×1.78-2.05) only widens this.
- N3 candidates on cells without N4 (1833 pairs): all inside the certified zero-free region (δ_k < 1 - edge(lk)), never evaluated.
  Row 1 (T <= e^29, δ <= 0.0073) lies wholly inside the zero-free region (edge δ 0.007374 at log T = 29).
- Max log T1 over the published scan = 36.2152 ([3120, 3130]): N4 cells exist for every evaluated height; rows 10-13
  (log T > 37) are never evaluated.
- c5 G-step: all 3320 cells have d5·wL >= 11.469 (min at cell k = 300, lk = 40; wL in [10.04, 24.14]) >= 9.75707: the printed
  |G| <= d5 W/log W exceeds the certified 9.75707 W/log^2 W, so c5 as printed is a valid upper bound (code/cells_g.py).

## 3. Reruns (code/qlab.py: Padhi's certify/arb/scan_arb.q_arb with his published kernels (c, f); labelled copy of his logN,
checked equal to the original on the first 3000 calls of each run; logs/qlab_*.log, logs/qv_*.log)
- (a) published constants, 126 intervals [43, 1100): every q equals scan_arb_{0,1}.txt (max |diff| 1.1e-16); max q on [57, 1100)
  0.991367809 on [872, 877]; argmin counts N1 190026, N2 150903, N4 3513525, N3 0. 246 s, peak 70 MB (memcap), RSS 101 MB.
- (b) 'both flaws' Table 2 C (uniform C ×2.046, T<=e^29 bound = row 1): [797, 1000) 27 intervals, identical, max 0.9913678095 [872, 877].
- (b) 'all fixes': identical, max 0.9913678095 [872, 877]; N3 deleted ('none'): identical on [797, 1000) (27 intervals).
- KLN ×2: max q 0.991916 [872, 877] (+0.00055); KLN deleted: q 852-1165 (KLN essential for σ < 0.97).
- NO_REP (d* = 0; 421 cells lk <= 32, k <= 140 recomputed, max Δlog N+ 0.4285 at 69,31): q 0.9356 [797,807], 1.0173 [827,837],
  1.0875 [847,857], 1.1429 [872,877], 1.1723 [877,882], 1.1735 [917,927] (lower bounds: other cells keep d*). Thm 5.2 load-bearing.
- n3bind per row (reached heights log T <= 36.2152): row 1 and T<=e^29 inside the zero-free region; rows 2-9 min log(N3/N4)
  7.863 (row 5, cell 66,33) .. 7.901; rows 10-13 never reached; uniform min 20.29.
- High range, N3 deleted, every 6th interval of [1100, 3159): 35 intervals identical, max q 0.759782 [1100, 1110]; N3 never.
- Extension (published N3, his top kernel c = 70, f = 0.40): q_arb [3150, 3159] 0.400012653 (= published), [3159, 3161]
  0.400011895, [3161, 3163] 0.400011758 (logs/ext3161.log). Transfer step not rerun there.
- Thm 5.2 sketch (code/rep/rep.py): F(z, η, σ) = Σ_m k_m W((z + (2σ−1)m)/η), W = Laplace transform of BTY's weight on [0, U];
  the second zero enters a1 Re(F(σ−1+η1−id) + F(σ−η1−id)) scaled by the primary η, so d* scales like η. Not audited: next target.
- d* halved (DSTAR_SCALE=0.5; 421 cells lk <= 32, k <= 140; max Δlog N+ 0.1936 at 120,32): q 1.0324 [847,857], 1.0610 [872,877],
  1.0742 [877,882], 1.0728 [917,927], 1.0507 [937,947] (logs/qv_dhalf_847.log; lower bounds). Peaks 62/59 MB (memcap).
- 19 Nov 2026 (session resumed after a 5-day pause): Zenodo 23440770 still has versions v1.0 (16 Oct) and v2 (17 Oct) only.
