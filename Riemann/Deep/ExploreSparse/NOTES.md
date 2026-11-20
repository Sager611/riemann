# explore61 working notes (appended as found; REPORT.md is the deliverable)
- W1 := {f ∈ L²(0,∞): f = a_m + b/x on I_m = (1/(m+1), 1/m) (m ≥ 1), f = b/x on (1,∞); one b}. Closed.
  ρ_{1/k} ∈ W1: on I_m, {1/(kx)} = (1/k)(1/x) − ⌊m/k⌋. χ ∈ W1. Hence V_N ⊂ W1 for all N.
  RH ⇒ V_∞ = W1 (χ ∈ V_∞ and contractions f(kx) keep V_∞; 1_(0,1/k) = χ(kx) span the piecewise constants; ρ_1 gives b).
  Conversely V_∞ = W1 ⇒ χ ∈ V_∞ ⇒ RH (Báez-Duarte). (Bagchi 2006's sequence-space form, presumably.)
- ρ_θ, θ ∉ 1/ℕ, jumps at θ/n inside I_m: ω(θ) := dist(ρ_θ, W1) > 0, and dist(ρ_θ, V_N) ≥ ω(θ) for every N (sampling defect).
- Split (exact, Pythagoras, V_N ⊂ W1): P = d_N² − D² = R + Q, R = dist(h*, W1)², Q = dist(P_W1 h*, V_N)².
- Reduction: d_N ≤ ‖χ − h‖ + dist(P_W1 h, V_N) for every h (P_W1 is a contraction toward χ ∈ W1): roughness is free.
- Lemma B: ⟨P_W r*, P_W h* − h_N⟩ = R (r* ⊥ C ⊇ V_N), so R² ≤ (D² − R)·Q (Cauchy–Schwarz); P ≤ ε(1+ε)D², ε = Q/R.
- s3_rate: dist(P_W1 ρ_θ, V_N)² = κ(θ)/log N to ±2 % for N = 10..2000 (θ = 2/3: κ = .0033; 2/5: .026; 2/9: .018).
- s1_split (head θ ∈ (1/5,1), grid Q/k): R/P = .72–.85; L·Q/R = 1.4–1.7 flat in N; cos(P_W r*, P_W h* − h_N) ≈ .95.
- s5_zero (200 zeros + tail): zero formula κ(θ) = Σ_ρ |M[P_W1⊥ ρ_θ](ρ)|² holds within 0–7 % for all 14 θ (e.g. 2/3: .00338 vs .00351;
  5/6: .00165 vs .00165; 2/5: .0263 vs .0269; 2/9: .0183 vs .0171). In-W1 cost = Burnol's zero sum for the roughness.
- Lemma B sharpened: R(1+ε) ≤ εD² ⇒ P = R + Q ≤ ε D², ε = Q/R. I.e. d_N² ≤ (1 + Q/R)·D(1/N)² (exact, unconditional).
  N=2000,Q=16: P/D² = .1935 ≤ Q/R = .208.
- s6_optzero (300 zeros + tail): along the hybrid optimizer Q·L = Σ_ρ |M[P_W1⊥ h](ρ)|² to 0.1–2.6 % (N=160, 640, 2000; grids 8, 16):
  L³Q measured/predicted .0534/.0534, .0589/.0582, .0527/.0513, .0750/.0739. Zero visibility Z = Σ|M g(ρ)|²/‖g‖² = 1.36–1.56 = L·Q/R.
- Mechanism: Q ≈ Burnol's zero bound for the roughness g = P⊥h* (lower-bound side expected from Burnol's finite-λ zero vectors, not checked; sharp numerically);
  with Lemma B (cos ≈ .95): P·L² → cos²·Z·C ≈ .92·1.5·.0462 ≈ .064 = the flat integer-dilation penalty.
- s5 (3000 zeros + tail) vs s3: κ(N=2000)/K0 = .96–1.06 over 14 θ (s10_cmp.log). DH (s8, 304 on-line zeros < 420): ratio .87–1.09;
  DH off-line floor slope Σ(2β−1)|G(ρ0)|² = 1e-5–1e-4 per unit L: invisible (≤ .4 % at N = 1000). Mechanism RH-blind at N ≤ 10³.
- s9 (Arb, 3000 intervals): dist(ρ_θ, V_N)² ≥ .13256 (2/3), .17548 (3/4), .04876 (2/5), .11340 (3/5), .01517 (2/9) for every N.
- Novelty: Bagchi math/0607733 (TeX read): piecewise-constant space M on (1/(n+1),1/n], g_l = {1/(lx)} − {1/x}/l ∈ M, RH ⟺ γ_l total:
  so W1 and RH ⟺ V_∞ = W1 are Bagchi's (known). arXiv API (6 queries), listings Oct–Nov 2026 (17 title hits, none on NB),
  Zenodo (74 NB records; 24382036, 24054105 abstracts read), GXL (2 claims, 0 supports): no P = R + Q split, Lemma B or zero formula for W1 targets.
- s13_tail (tail midpoints 2/(2k+1), k ≤ 300): R/P ≈ .81–.83, c_e = 1.56 flat in N, cos ≈ .89–.91; tail = 81–95 % of ‖e‖₁, +4 % of P.
- Reviewer (verdict only): W1 facts, defect, split, Theorem C CORRECT; Lemma B GAP (R > 0) fixed as R·P ≤ D²·Q.
