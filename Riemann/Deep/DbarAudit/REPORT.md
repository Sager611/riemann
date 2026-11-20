AUDIT PASSED — every step of the chain (dual form, grouping, Lemmas A, B, Λ, D, E, Theorems R″, R′, Corollary, §6 tilings and (3v), explore59's (3i), (3ii), §7, Claim F, Burnol's liminf) re-derived line by line with its inputs stated from the sources; no false step, no counterexample in RH-allowed worst cases; a fresh reviewer on the most delicate step (Lemma Λ + Lemma D): CORRECT.

# DbarAudit — dbaraudit (Orca, 20 Nov 2026, 11:08–11:42 UTC)

## Claim audited
RH ⇒ lim_{λ→0} D(λ)² log(1/λ) = Σ_ρ m_ρ²/|ρ|², D(λ) = dist_{L²(0,∞)}(1_{(0,1)}, span{{θ/x} : λ ≤ θ ≤ 1}) (Riemann/Deep/DbarNB/,
built on Riemann/Deep/ExploreDbar/). Chain with line numbers and my checks: FULLPROOF.md §0–§13.

## Result by brief item
- (a) Group bound: groups are runs of consecutive ordinates with gaps < s₀; the first ⌊M*(t)⌋+1 of a group lie in an interval of length
  ≤ 7A₁ log L·L^{−3/4} < 1/log log t (t ≤ L⁶+4), so (S) (Titchmarsh 14.13, S(t) ≪ log t/log log t) gives M_g ≤ A₁ log t/log log t,
  uniformly for L ≥ L₀. Turán–Nazarov: Friedland–Yomdin 1107.0039 Thm 1.1 (TeX l.72–78): m+1 terms ↦ exponent m, absolute c, factor
  e^{|B|max|Re λ|} = 1 for real frequencies. Jets: not stated in FY; follow by continuity (divided differences of real exponentials →
  jets uniformly on compacts), exponent = number of terms with multiplicity − 1. ✓
- (b) Lemma B: only W_g/γ_g enters (W_g ≤ 2M_gL^{−3/4}); γ² gain uniform in the width. Numerics: endpoint ratio max at the coalescent jet
  (1.42 … 2.86e4 for M = 1 … 8, ≈ 3.4^M), finite γℓ never exceeds the limit. ✓
- (c) Lemma Λ/D: every leakage bound is relative to the group's own ‖q_g‖_χ (TN by-product handles ill-conditioned clusters), at most two
  windows per distance shell, so nothing depends on the number of groups (~L⁶log L). Reviewer: CORRECT. Numerics: cos of the principal
  angle between neighbouring groups decays ≈ e^{−√(s₀ℓ)} with an M-dependent prefactor; K3 → 1 (1.0004–1.18 at s₀ℓ = 100). ✓
- (d) R′: S_hi ≤ ε(H) (tail of a convergent series) uniformly in L; θ_L = C(H)L^{−3/4+o(1)} for fixed H. R′ only for fixed H (as stated). ✓
- (e) §6: residue sums depend only on enclosed zeros, so narrow/wide rectangles coincide and the two height families differ by one contour
  integral → 0 (no 1/ζ′(ρ) bound); divided-difference form valid for s ≠ 1, R_high(1) directly (1 outside all rectangles); trivial-zero
  residues r_k; growth in all four regions, exponent L′ + 2/T₁ ≤ L − ½; (3v) unconditional. ✓
- (f) Burnol math/0103058 TeX l.169–193: K = L²(]0,∞[), finite combinations of ρ(θ/t), λ ≤ θ ≤ 1, no constraint; l.203–206 liminf D√log(1/λ)
  ≥ √Σ m_ρ²/|ρ|², unconditional (l.197–201: +∞ if RH fails). Same D(λ). zbMATH 1008.11032 (Notes III) confirms the setting. ✓
- (g) (3i) window (spectrum [0, 2/T₁], |Φ| ≤ 1, ‖1−Φ‖² ≪ L⁻²), (3ii) trivial zeros, §7 Bochner closure, Claim F (Legendre endpoint
  Christoffel number m²): re-derived. ✓

## Falsification (numerics; FULLPROOF §12; logs/)
- audit_num.py B, D (mpmath Newton–Hermite bases, double QR, composite Gauss quadrature, h-converged to 1e-14), audit_chain.py
  (s₀-chains m = 6, 8; dense 4-clusters at every 2nd zero): L·U bounded, cost/decoupled-sum → 1. dbarnb_falsify.py dual and clusters
  (m = 5, h = 0.3) rerun from a copy (rerun/): identical to DbarNB's logs. audit_mpcheck.py: mpmath (dps 30) = double to 12 digits
  for cost ratios near 1 (0.99966, 0.99986: real, groups help slightly). Peaks ≤ 227 MB (one D run killed at 520 MB by memcap:
  leggauss(4800) companion matrix; redone with 48-node panels).
- Lemma D is invisible at computable L: for 4-fold and 8-fold jets the groups are numerically dependent (K3 ≥ 1e10) until s₀ℓ ≈ 30–100,
  i.e. L ≈ 10⁶–10⁸. No contradiction (the lemma is for L ≥ L₀); it quantifies the source's "ineffective" caveat.

## Novelty spot-check (brief (4))
- Balazard–Saias, Notes 4 (Adv. Math. 188 (2004) 69–86): zbMATH Open review Zbl 1096.11032 (A. Ivić, via api.zbmath.org): Báez-Duarte's
  sequence-space criterion d²(1, E_n) in ℓ²(ℕ*, 1/(k²+k)) (discrete), the conjecture d² ~ (2+γ−log 4π)/log n, an exact constant A for
  (I−S), and a limit in another space. No continuous D(λ), no RH upper bound: "Before" unchanged.
- Also Zbl 1008.11032 (BBLS, Notes III: liminf D√log(1/λ) > 0 and the conjectured value) and Zbl 1173.11343 (BBLS 2005, A(λ) only).
  GXL (1 call): 0 supports among 254 retrieved. Note: scout bs4b asserted that Notes 4 proves D(λ) = O((log 1/λ)^{−1/2}) under RH but
  found no citing sentence; the zbMATH review refutes it (scout claim discarded).
- With the theorem, RH + simple zeros ⇒ the BBLS conjecture for D(λ) (lim D²log(1/λ) = 2 + γ − log 4π); with Burnol's m² bound,
  that conjecture ⇔ RH + all zeros simple.

## Points the referees did not record (none breaks the chain; FULLPROOF §13)
1. Turán–Nazarov's jet case rests on a continuity remark not in FY; numbering "Nazarov Thm 1.3/1.4" unchecked. 2. R″/R′ are proved for
Z(L⁶) but used for Z(T₂), T₂ ≤ L⁶+4 (same proof). 3. Lemma B's 8/ℓ can be 4/ℓ. 4. L₀ ≥ 10⁸ in practice (above).

## Black boxes (not re-proved)
Nazarov's lemma (absolute constant); Titchmarsh Thm 9.7 (book statement; no open verbatim quote found); RH short-interval count
(Littlewood; verbatim in Goldston–Gonek math/0511092 TeX l.107–122, which even gives the constant ½ + o(1)); Burnol's
proof of the liminf (Adv. Math. 170 (2002)); Paley–Wiener–Schwartz. Ineffective: no L₀, no rate.

## Files (Riemann/Deep/DbarAudit/)
REPORT.md; FULLPROOF.md (chain with source line numbers, §12 numerics, §13 remarks); audit_num.py (modes B, D); audit_chain.py; audit_mpcheck.py;
rerun/ (copies of dbarnb_falsify.py, dbar_U.py, import path pointed at the copy); logs/ (B, B_v1_singlepanel, B_v1_vs_v2, D,
D_killed520MB, chain, mpcheck, rerun_dual, rerun_clusters_m5_h0.3, with memcap and /usr/bin/time footprints); evidence/ (zbMATH reviews of
Notes III, IV, BBLS 2005, Balazard–Saias 2000; Burnol and Friedland–Yomdin TeX extracts with line numbers; reviewer verdict on Λ + D).

RH relevance: N5
