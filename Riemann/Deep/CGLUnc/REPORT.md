OBSTRUCTION: Chirre–Gonçalves–de Laat's semidefinite gain cannot enter zeta23's unconditional frame: both its inputs hold without RH (F ≥ 0 and F ~ |α| on [−1,1] for complex γ, Baluyot–Goldston–Suriajaya–Turnage-Butterbaugh 2306.04799 Thm 1), but its diagonal step S(K) ≥ K(0)Σm² fails and cannot be repaired in the configuration-wise frame: any CGL-type kernel (band-limited, or CGL's polynomial × Gaussian) with K̂ < 0 beyond |α| = 1 gives an off-line pair at depth y the self-energy 2K(0) + 2K(2iy) → −∞ faster than every bandwidth-one cost, so even one deep pair (or N^θ shallower ones) makes CGL's inequality vacuous (for band-limited K̂ only a zero-density slope > 2 would stop this, known 8/7; for CGL's own Gaussian-type kernels no density slope does, only thin-box hypotheses). BARRIER: the arena's TEN class (every certificate that yields the arena's universal affine-defect estimate — any block size, weights, pressures, floors — on any window cos(√2 s) + Σ_{j≤J}c_j cos(2πjs), J ≤ 96, their mixtures or the best of several) is capped at 153/227 = 0.6740088 (Arb dual), the record's own window at 31/46 = 0.6739130 (float limits 0.67394 and 0.67384), against the record 0.6736393; every bandwidth-one certificate, CGL-type or not, stays ≤ 0.6818287 (Alpöge–Furman's Lean ceiling).

# CGLUnc — report (cglunc, Orca, 12 Nov 2026, 14:12–15:50 UTC)

Files: DERIVATION.md (steps D1–D7), PROOF.md (Props U, O, Lemma P, Prop D), code: tenwin.py, periodic.py, necklace.py,
winopt.py (+winopt.log), dual_build.py, certify_dual.py (+certify_dual_TEN.log), cover.py, verify_cover.py
(+verify_cover_J16.log, verify_cover_J48.log), coverJ.py, crosscheck*.py, marks.py (+marks_*.log), cover_marks.py; data: dual_TEN.json, cover_*.json, winopt_best.json.
All runs single-threaded under memcap 290 MB (peaks ≤ 127 MB; cover_marks.py stopped by its 1200 s watchdog, result in D7).

## Answer
1. zeta23's 2/3 is one inequality: N0s ≥ rank P₁ ≥ 4 tr G̃ − 2N − ‖G̃‖²_F = (2 − R(ψ))N for a Gabor–Gram matrix of Weil's
 form; TEN adds the configuration-wise Gram defect D∘ of the simple on-line block: κ = (H_W − B)/(1 − c).
2. CGL's two inputs are unconditional (BGSTB 2306.04799 Thm 1); its diagonal lower bound fails as soon as one deep
 off-line pair is present, and no configuration-wise variant survives (Prop O, Lemma P): CGL-type information adds
 nothing to bandwidth-one certificates unconditionally.
3. The TEN class has a hard cap 0.6740088 (certified; ≈ 0.67394 float): records of this type can add at most 3.7e−4;
 any one-window refinement (Q-side defects, cross terms) ≤ 0.6777 on the record's window, ≤ 0.6765 on MT (float; no
 window-independent adversary found, so a family-wide one-window cap stays open); every bandwidth-one certificate ≤ 0.6818287 (AF).

## 1. zeta23 (Alpöge–Furman 2608.13637) as mathematics
Window ψ ≥ 0 on [−1/2,1/2]; φ ∈ C_c², supp φ = [−L/2, L/2], φ² ≈ ψ(u/L), L = λ log(T/2π) (λ ≤ 1, X = e^L ≤ T).
Grid α_k = T + 2πk/L (0 ≤ k < d ≈ N(T,2T)); v_ρ = (φ̂(γ_ρ − α_k))_k, γ_ρ = (ρ − 1/2)/i;
G̃ = (aL²)⁻¹Σ m_ρ v_ρ v_ρᵀ, real symmetric because v_{1−ρ̄} = conj v_ρ. Poisson–Gabor (exact at density h = 2π/L):
Σ_k φ̂(z−α_k)φ̂(z'−α_k) = L(φ²)^(z − z') for complex z, z'.
- Zero side (block structure): on-line zeros give rank-one PSD terms (P ⪰ 0, rank P ≤ N0*); an off-line pair gives
 2m(aaᵀ − bbᵀ), the pull-back of m·diag(1, −1), so n₊(Q) ≤ #pairs (Sylvester). Off-line zeros are never estimated one
 by one: their huge entries (|φ̂(x+iy)| ≤ L e^{L|y|/2}) only enter through n₊(Q) and the global tr, ‖·‖².
- Trace (EF + RvM + N(t,t+1) ≪ log t): tr G̃ = N + o(N).
- Frobenius (EF → prime side; Chebyshev–Mertens Σ_{n≤X}Λ²/n for n = n', Montgomery–Vaughan for n ≠ n', X ≤ T):
 ‖G̃‖² = Σ m m' K(Δ) = (R(ψ) + o(1))N, K̂ = ψ⋆ψ/(∫ψ)² on [−1,1], R(ψ) = (∫ψ² + ∬|u−v|ψψ)/(∫ψ)².
- Rank–trace (von Neumann; matrix form of (x−1)² ≥ 0, x² ≥ 0, (x−2)² ≥ 0): r ≥ 2trP + 4trQ − 4b − ‖P+Q‖².
 With P₁ = simple on-line part and doubles/pairs at flat charge 4: N0s ≥ 4trG̃ − 2N − ‖G̃‖² = (2 − R)N.
 ψ₀ = 1: R = 4/3 → 2/3; R_λ(ψ₀) = 1/λ + λ/3 gives the profile H(λ) = 2 − 1/λ − λ/3, maximal at λ = 1;
 ψ_MT = cos(√2 s): R = c_MT⁻¹ = 1.3274992, κ₀ = 0.6725007, optimal for 2 − R(ψ) (CCLM 1406.5462 Cor. 14).

## 2. TEN (arena record knausgard-ten-win10-lp)
Exact identity: for the retained simple on-line zeros, r = 2trΓ − ‖Γ‖² + D∘ with D∘ = ‖Γ − I‖² ≈ Σ_{i≠j}k_W²(x_i − x_j),
k_W(x) = ψ̂(2πx)/∫ψ. MT/CCLM discard D∘; TEN bounds it configuration-wise: a 10-point (9-gap) local certificate with
45 pair weights a_ij, gap pressures b_r (→ B·span), a floor c per point and a reflection half-space row proves
D∘ ≥ c·r − B·span − O(1) for every ZeroConfig (Lean ActualDefectEstimate, l.16111). Window
ψ = cos(√2 s) + Σ_{j≤16}c_j cos(2πjs) trades H_W = H_MT − Σ w_jc_j² (here 0.672063, −4.4e−4) for larger k_W at the
near-integer distances the adversary uses. Record: c = 0.0093733, B = 0.0047380, κ = (0.6720630 − B)/(1 − c) = 0.6736393.
Side finding: the Lean proof passes through (κ − ε)N ≤ N0simple for ζ (actual_dyadic, l.16237–16241, then
N0simple ≤ N0star); so the arena constant also holds for simple zeros on the line, beating 2610.08965's 0.6735382 there
too. TREE N8's "not for simple zeros on the line" should be corrected.

## 3. Can a CGL-type function enter? (exact failing step)
| CGL step (1810.08843) | RH-world input | Unconditional replacement | Status |
|---|---|---|---|
| F(α) ~ T^{−2α}log T + α on [0,1] (l.474–478) | Montgomery, Goldston–Montgomery | BGSTB 2306.04799 Thm 1 with x^{ρ−ρ'} (= AF's Thm HS, polarised over windows) | holds |
| F(α) ≥ 0 for |α| > 1 (l.446) | Mueller, Heath-Brown | BGSTB Thm 1; reason: ρ ↦ 1−ρ̄ gives Σm m'e^{iξ(γ_ρ−γ_ρ')} = |Σ m e^{iξγ_ρ}|² | holds |
| upper bound S(K) ≤ N(K̂(0) + ∫_{−1}^{1}K̂|α|) for K̂ ≤ 0 off [−1,1] | the two lines above | the same | holds (Prop U) |
| lower bound S(K) ≥ K(0)Σ_γ m_γ² (l.~500, eq:oneineq) | real ordinates, K ≥ 0 | none | FAILS (Prop O) |
| rank–trace with ‖·‖²_F in place of S(K) | — | only positive-definite pair kernels | F ≥ 0 gives the wrong direction (Lemma P) |
Prop O (K̂ ≤ 0 on |α| ≥ 1, ≢ 0 beyond 1; class C1: band-limited C^k_c([−A,A]), class C2: CGL's p(x)e^{−πx²}; both have
the exponential moments S(K) = N∫K̂F needs for complex zeros): an off-line pair at depth y = (β−1/2)L/2π has self-energy
2K(0) + 2K(2iy), K(2iy) = ∫K̂cosh(4παy) ≤ ‖K̂‖₁cosh 4πy − I₂(y), I₂ ≥ c e^{4πAy} (Gaussian-type CGL kernels: e^{cy²});
n_p ≍ N^θ such pairs (any θ ∈ [0, 1/2]; in C1 even ONE pair at β ≈ 1/2 + 1/(2A)), added to any shallow (e.g. on-line)
configuration, lower S(K) by ≥ MN while moving N, N0*, N0s, the trace and every bandwidth-one second moment by o(N)
(padding block E: ‖E‖²_F ≪ N^{θ+(1−θ)/(A−ε)} in C1, N^{θ+o(1)} in C2, then Cauchy–Schwarz for 2⟨G̃, E⟩). So CGL's lower bound is false for padded configurations, and its upper bound cannot
lower any bound coming from shallow adversaries (AF's on-line law included). Scope: configuration-wise arguments on
bandwidth-one data and counts, i.e. the zeta23/TEN frame; not arguments with new arithmetic input.
Zero density: in C1 the padding sits at β − 1/2 ≈ (1−θ)/(2A), allowed by N(σ,T) ≪ T^{1−a(σ−1/2)}log T whenever a < 2A;
a slope a > 2A (> 2) would exclude it (known 8/7); in C2 it sits at |β − 1/2| ≍ (log T)^{−1/2} and is compatible with
every slope: only a box hypothesis without exceptions (as in BGSTB and Goldston–Suriajaya, |β − 1/2| ≤ b/log T) stops it.
Freedoms tested (lesson 191): analytic weights η (|η(γ)|² multiplies the pair term; log|η| subharmonic, so no
real-analytic η is small on Im z = ±y and ≈ 1 on ℝ), complex shifts K(Δ + iσ) (K̂ → K̂e^{−2πσα}, same tail), complex or
sign-changing windows and several windows (still positive-definite HS kernels), Gaussian vs compact tails (both fail),
more points/defects (bandwidth-one, covered by AF's ceiling). Under GRH CGL add F ≥ 3/2 − |α| on [1, 3/2): not available.
The inequality one CAN prove is Prop U; its best unconditional κ is the bandwidth-one value itself (CGL adds 0).

## 4. Ceilings (dual certificates)
Dual lemma (Prop D): a periodic configuration X of the retained points with density ν and Gram-defect energy
e_W(X) ≤ ν − H_W forces every TEN-type certificate on window W (admissible c < 1) to κ ≤ ν — X is realised as an
on-line simple ZeroConfig Z_X (ordinates with (γ/2π)(log(γ/2π) − 1) ∈ SZ + X), for which the arena's
ActualDefectEstimate (proved in Lean for every ZeroConfig) must hold; e_W is an exact finite sum (Poisson,
K̂ supported in [−1,1]). The extremal configurations are near-integer patterns (.xx at ν = 2/3, e = 1.697e−3;
.xx.xx.xxx, .xx.xx.xx.xxx), i.e. points sitting near the near-integer zeros of k_W.
| Class | Bound on κ | Type | Room above 0.6736393 |
|---|---|---|---|
| TEN class on the record's window | 31/46 = 0.6739130 (float hull 0.673840) | Arb, one configuration, slack 7.18e−5 | ≤ 2.7e−4 (2.0e−4) |
| TEN class on the MT window | 0.673522 | float, necklaces S ≤ 12 | none (window change needed) |
| TEN class (universal affine-defect estimate), windows cos(√2 s) + Σ_{j≤J}c_j cos(2πjs) with J ≤ 96, mixtures, best-of-several (complex coefficients only via Hermitian G, |k|²) | 153/227 = 0.6740088 (float optimum 0.67394–0.67397) | Arb LDLᵀ, 12 configurations, J = 16, 48 here (min pivot 3.85e−4), J = 96 orchestrator rerun (3.8468e−4); one common periodic configuration | ≤ 3.7e−4 (3.0e−4) |
| (b) CGL-type class | = value of bandwidth-one data alone | Prop O | — |
| any bandwidth-one certificate for N0s (TEN included) | 0.6818287 + 2.55e−6(|r′(1)| + ∫|r″|) | Alpöge–Furman, Lean PairCeiling (given EnclOK) | ≤ 8.2e−3 |
| reference, under RH | simple ≥ 0.6792 (CGL) | | |
N0* (the arena's count): TEN certifies N0s, so the same caps apply. For general N0* certificates AF's law with every
even mark replaced by near-line off-line pairs gives ≤ 0.6818287 + (fraction of odd marks ≥ 3), not computable here
(the law cert_N256 is not public).
Escape routes from the TEN cap (honest status):
| Route | What it adds | Status |
|---|---|---|
| any further configuration-wise lemma on ONE window's data (Q-side defects Σ(q_j − 2)², cross term 2tr(P₁Q₊), von Neumann slack) | penalises doubles/pairs next to simple zeros | capped by on-line double adversaries matching that window's HS value (marks.py, patterns S ≤ 9–10, LP mixtures, float): ≤ 0.676483 (MT), ≤ 0.677662 (record window), ≤ 0.677965 (float-optimal window); room ≤ 4.0e−3 there; no single mixture works for all windows (cover_marks.py, S ≤ 8 library: LMI margin −0.018), so other windows stay open |
| bandwidth-one data of several windows used as constraints on the pair statistics (global LP over laws) | excludes near-lattice adversaries (their form factor is lattice-like) | AF's LP; ceiling 0.6818 |
| higher moments tr G̃^k | | dead: X^k ≤ T^{2−ε} (AF §7(e)) |
| CGL positivity beyond |α| = 1 | | dead unconditionally (Prop O) |

## 5. Verification and novelty
Own checks: H_W closed form = quadrature (3 random windows, 1e−14); quadratic-form slack = direct slack (1e−14);
Poisson energy = real-space sum (2e−7, tail); H_W of the record window = the Lean constant to 25 digits.
Reviewers (fresh, verdict-only, decisive steps pasted inline): see §6. Novelty: axq ('pair correlation' ∧ zeros ∧ critical line; semidefinite ∧ zeta; simple zeros ∧
unconditional; Montgomery ∧ optimal), lst ('pair correlation|simple zeros|distinct zeros|critical line', 2026-10/11, Zenodo
"pair correlation"): no ceiling for the arena/TEN class and no unconditional CGL import anywhere. Known and used: BGSTB Thm 1
(F ≥ 0 and the asymptotic unconditionally, 2023; they need thin boxes |β − 1/2| < 1/(2log T) for 61.7% simple, GS 2501.14545
b = 0.3185 for 2/3), AF's bandwidth-one ceiling 0.6818287. R453's 'no unconditional import of CGL' stands only for the
lower bound: the inputs were already unconditional in BGSTB.


## 6. Review
Repairs (orchestrator referee, Prop D/D2: GAP repairable): (1) class = the arena's ActualDefectEstimate for every
ZeroConfig (Solution.lean l.16111–16115, l.16204–16224; ZeroConfig fields zeta23 Defs.lean l.120–130), realised by an
on-line simple ZeroConfig Z_X, O(1) = (n−1)c uniform, limits T → ∞ then δ → 0; (2) tr(MG) > 0 only for the augmented
normalised Gram G ⪰ 0, G₀₀ = 1, with certified rational λ_i ≥ 0, Σλ_i = 1; (3) one common configuration X̄_B via
‖Q_{X̄_B} − Σλ_iQ_{X_i}‖_op ≤ (J+1)C_J log(227B)/(227B) and the PD margin; (4) complex coefficients only via Hermitian G
and |k|² energies; (5) family narrowed to J ≤ 96 (J = 96 by the orchestrator), uncovered freedoms listed in PROOF.md.
- Round 1 (both claims): Claim B (TEN-class dual lemma, Poisson normalisation S⁻², S-procedure, Gram/complex
 coefficients) CORRECT for the stated class (J ≤ 48, admissible c < 1; not unrestricted harmonics). Claim A GAP:
 exponential moments needed for K(2iy); padding block needs ‖E‖²_F = o(N) plus Cauchy–Schwarz; scope.
- Round 2 (repaired A): GAP — integration-by-parts constant; zero-density "iff" false for Gaussian kernels.
- Round 3: GAP — C1 constants (e^{o(y)}, two-sided pinning of y); "all but ≪ √N" boxes do not block.
- Round 4 (final, PROOF.md Prop O verbatim): CORRECT ("both bounds on I₂ pin y; bounded-pair padding survives a box
 allowing a pair of exceptions"); one harmless rounding (= → ≤) applied.
Not reviewed: AF's ceiling (theirs, Lean modulo EnclOK); the float numbers (window optimum 0.67394–0.67397, one-window
caps 0.6765/0.6777, hull limit 0.67384) are evidence, not certificates.

## 7. Notes for the orchestrator
- TREE N8 says the arena bound "beats 0.6735382 for zeros on the line, not for simple zeros on the line": the accepted
 proof establishes (κ − ε)·N(T,2T) ≤ N0simple(T,2T) for ζ as a kernel-checked intermediate (`hs` in
 OpethN.ActualDefect.actual_dyadic, Solution.lean l.16231–16243, via defect_to_dyadic l.16116 with zetaZeroConfig and
 paperInputs_zeta), then N0simple ≤ N0star. So 0.673639288 also holds for simple zeros on the line (audit before
 recording; the definition N0simple is identified with zeta23's N0s by zetaZeroConfig_N0s).
- Records of TEN type cannot pass 0.6740088; a record beyond needs information the class discards: Q-side defects or
 cross terms on one window (≤ 0.6777, float), or the bandwidth-one data of several windows used as constraints on the
 pair statistics (≤ 0.6818287). Not a DISCOVERIES candidate by itself (method ceiling); NOGO.

## Suggested NOGO line
- cglunc (Orca, N8, 12 Nov 2026): Chirre–Gonçalves–de Laat's semidefinite pair-correlation gain cannot enter zeta23's unconditional frame — F ≥ 0 and F ~ |α| on [−1,1] hold without RH (BGSTB 2306.04799 Thm 1), but every CGL-type kernel (K̂ < 0 beyond |α| = 1, band-limited or polynomial × Gaussian) gives an off-line pair self-energy → −∞, so one deep pair (or N^θ shallower ones) voids the diagonal bound at o(N) bandwidth-one cost (band-limited: only a density slope > 2A would stop it; Gaussian: only exception-free thin boxes); the arena's TEN method (universal affine-defect estimates) is capped at 153/227 = 0.6740088 over all windows cos(√2 s) + Σ_{j≤96}c_j cos(2πjs), mixtures and best-of-several, at 31/46 on the record's window (Arb duals; record 0.6736393), one-window refinements on the record's window ≤ 0.6777 (float) [CGLUnc](Riemann/Deep/CGLUnc/REPORT.md)

RH relevance: N8
