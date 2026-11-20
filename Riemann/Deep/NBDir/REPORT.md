THEOREM: for every real primitive χ mod q ≥ 3 (and, §8, every complex primitive χ), GRH for L(s,χ) ⇒ lim_{λ→0} D_χ(λ)² log(1/λ) = Σ_ρ m_ρ²/|ρ|²; so lim = 2Re L′/L(1,χ) + log(q/π) + ψ((1+a)/2) ⟺ GRH for χ and simple zeros. Lower half known (de Roton 2006); upper half a routine transfer of our ζ proof (corollary by the gate).

## Statement
χ primitive mod q ≥ 3 (nonprincipal, conductor q), parity a (χ(−1) = (−1)^a), A_χ(y) = Σ_{n≤y}χ(n), L = log(1/λ),
D_χ(λ) = dist_{L²(0,∞)}(1_{(0,1)}, span{A_χ(θ/x) : λ ≤ θ ≤ 1}) (Mellin of A_χ(1/x): L(s,χ)/s, entire).
(i) Unconditional (Theorem Bχ = de Roton, Bull. SMF 134 (2006) 417–445, Thm 2, Selberg class): liminf L·D_χ(λ)² ≥ Σ_{Re ρ = ½} m_ρ²/|ρ|²;
    a nontrivial zero off the line forces D_χ(λ) ≥ c_χ > 0.
(ii) GRH for L(s,χ) ⇒ lim L·D_χ(λ)² = Σ_ρ m_ρ²/|ρ|² (distinct zeros, multiplicities m_ρ). No simplicity, no condition on L′(ρ,χ); L(½,χ) = 0 allowed.
(iii) Equivalent form: lim L·D_χ(λ)² = C_χ := 2Re L′/L(1,χ) + log(q/π) + ψ((1+a)/2) ⟺ GRH for L(s,χ) and all its zeros simple.
(iv) Complex χ: (i)–(iii) hold verbatim (phase W(χ)^{−1}(s/(1−s))γ_a(s); no (χ, χ̄)-pairing; two-sided good heights).
Constants: C_{χ₋₄} = 0.155567980 (brief: 0.15556798), C_{χ₅} = 0.156556954. Not uniform in q (L₀ depends on χ).

## Proof
DERIVATION.md (§1 phase operator V_χ = W^{−1}(M − 1)D_qF_±I; §2 pointwise kernel lemma Bχ; §3 vectors Y; §4 orthogonality, Gram, target;
§5 Theorem Bχ; §6 upper bound: DbarNB §1–7 + ExploreDbar (3i), (3ii), §7, Claim F with ζ → L(s,χ), changes (c1)–(c6); §7 theorem, constant,
equivalence, uniformity; §8 complex χ).
New ideas beyond our ζ proof (the gate): none needed.
- Lower half: Burnol's construction needs the phase conj(â)/â explicitly; for χ it is W^{−1}(s/(1−s))γ_a(s) by the functional equation:
  even χ: V_χ = D_q·U with U Báez-Duarte's operator itself (Burnol/BBLS verbatim, dilated); odd χ: the sine analogue, kernel
  ∂_v[(1 − cos(2πt/v))/(πt/v)], whose pointwise lemma is the same computation (∫_0^∞x^{μ−1}(1 − cos x)dx = −Γ(μ)cos(πμ/2)). One (M − 1)
  instead of Burnol's three (â_χ = L(s,χ)/s is already in H²). And de Roton had done the whole Selberg class in 2006.
- Upper half: verbatim; the non-mechanical points are (c2) the trivial zero s = 0 of even χ (the existing R_tr kernel e^{−(L−1)(s−s_k)}/(s−s_k)
  gives ‖L R_tr‖² ≪ e^{−(√L−2)}), (c6) a possible central zero (Lemma E needs the full N-Gram, nondegenerate at γ = 0), and for complex χ
  heights good at ±t simultaneously and actual residues in the lower half-plane.
Uses of GRH (only in the upper half): (a) real ordinates: Turán–Nazarov without growth factor, Gevrey decay at real arguments (Lemmas A,
B, Λ, D, E; Claim F needs real nodes); (b) S(t,χ) ≪ log qt/log log qt (cluster bound, §3(iii)); in the assembly, (c) all nontrivial zeros
on the line, so (i) gives the full sum. The reduction (3v)χ is unconditional.
Cited inputs, exact hypotheses: Friedland–Yomdin arXiv:1107.0039 Thm 1.1 (Turán–Nazarov; real frequencies; jets by the continuity argument
of DbarAudit §0); Titchmarsh Thm 9.7 (unconditional) and 14.13 (RH) — their proofs for L(s,χ), fixed χ, using Davenport ch. 16
(L′/L = Σ_{|t−γ|<1} m_ρ/(s−ρ) + O(log q(|t|+2)) on −1 ≤ σ ≤ 2; N⁺(T,χ) = (T/2π)log(qT/2πe) + S(T,χ) + O(1)); Davenport ch. 9
(functional equation, W(χ) = τ(χ)/(i^a√q)) and ch. 12 (Hadamard product, Re B(χ) = −Σ m_ρ Re(1/ρ)); convexity |L(½+it,χ)| ≪ (1+|t|)^{¼+ε};
Paley–Wiener–Schwartz; Burnol math/0103058 (framework, l.296–672); de Roton 2006 Thm 2 (for (i), as a second source).

## Evidence (falsification first; all single-threaded, memcap peaks ≤ 82 MB, logs/)
Upper bound D_h² (θ-grid θ_j = e^{−jh}; exact Gram: K(r) = ∫_1^∞A(rz)A(z)z^{−2}dz = c₀L(1,χ) + r^{−1}Σ_{n≤N}χ(n)V(rn)/n², checked against brute-force
piecewise integration to 15 digits (r = 1, 1.37, e, e²); N = 32768 vs 131072: max |ΔĜ| = 4e-9; Toeplitz solve = Cholesky to 1e-9; Gram
eigenvalues in [4e-3, 66], coefficient norms ≈ 1: well conditioned). L·D_h², h = 0.1 → 0.0016 (each halving nests the grid, so values decrease):
| χ | L | h = 0.1 | 0.025 | 0.00625 | 0.0031 | 0.0016 | Burnol LB (zeros |γ| < 60 / 45) | C_χ |
|---|---|---|---|---|---|---|---|---|
| χ₋₄ | 4 | 0.1359 | 0.1177 | 0.1098 | 0.1067 | 0.1053 | 0.0917 | 0.15557 |
| χ₋₄ | 8 | 0.2282 | 0.1639 | 0.1359 | 0.1272 | 0.1229 | 0.1037 | |
| χ₋₄ | 12 | 0.3293 | 0.2129 | 0.1587 | 0.1424 | 0.1345 | 0.1105 | |
| χ₋₄ | 16 | 0.4351 | 0.2670 | 0.1839 | 0.1589 | 0.1469 | 0.1186 | |
| χ₅ | 4 | 0.2180 | 0.1855 | 0.1640 | 0.1615 | 0.1590 | 0.0828 | 0.15656 |
| χ₅ | 8 | 0.2961 | 0.2013 | 0.1433 | 0.1360 | 0.1291 | 0.0989 | |
| χ₅ | 12 | 0.4335 | 0.2693 | 0.1676 | 0.1537 | 0.1415 | 0.1057 | |
| χ₅ | 16 | 0.5756 | 0.3389 | 0.1908 | 0.1695 | 0.1510 | 0.1077 | |
Grid convergence: halving h at h = 0.0031 changes L·D_h² by 0.0015 / 0.0043 / 0.0079 / 0.0120 (χ₋₄, L = 4, 8, 12, 16) and 0.0025 / 0.0069 /
0.0122 / 0.0185 (χ₅): not converged (error ≈ h^{0.5–1}: the grid resolves ordinates only up to ≈ π/h, D(λ) uses them up to ≈ e^L/q).
Aitken on the last three χ₋₄ values (geometric differences) gives 0.104, 0.119, 0.127, 0.136 at L = 4, 8, 12, 16 (indicative only; χ₅'s
differences are irregular, no extrapolation).
Lower bound (rigorous at each finite λ up to quadrature): L·‖P_R 1‖² for Burnol's vectors X^λ_{ρ,0} of the zeros in R (Lemma Bχ closed forms
for Φ_ρ, Gauss–Legendre on (λ/q, 1] and [1, 40] plus analytic tails): χ₋₄, R = 50 zeros (|γ| < 60): 0.0917, 0.1037, 0.1105, 0.1186, 0.1194,
0.1240, 0.1280, 0.13050, 0.130732 at L = 4, 8, 12, 16, 24, 40, 100, 10³, 10⁵ → Σ_R 1/|ρ|² = 0.130739 (like Σ_R − 0.27/L); χ₅, R = 36 zeros
(|γ| < 45): 0.0828 … 0.123819 (L = 10⁵) → 0.123823; |√L(1, X_ρ)|² → 1/|ρ|² (0.027396, 0.022496) as §4(T) predicts. χ₋₄ with R = 100 zeros
(|γ| < 100): 0.0964, 0.1094, 0.1169, 0.1252, 0.1265, 0.1317, 0.1360, 0.13882, 0.139089 → 0.139101. So L·D_{χ₋₄}(λ)² ∈ [0.1252, 0.1469] at L = 16
([0.1169, 0.1345] at L = 12), and the Aitken values 0.104, 0.119, 0.127, 0.136 lie inside the brackets.
Verdict: no counterexample. At every L tested LB < UB, and both lie below C_χ, increasing in L: consistent with an approach to C_χ from
below as in the ζ toy models (Σ − c/L). The numerics cannot pin the limit at L ≤ 16 (as for ζ: L₀ is astronomically large).

## Checks
Complex χ (§8): χ mod 5 of order 4 (χ(2) = i, odd), s = ½ + 3i: W q^{−½}∫_1^∞A_χ̄(t)t^{−1}conj φ_s(t/q)dt = 0.15138123 − 0.62666322i
vs L(s,χ)/s = 0.15138123 − 0.62666322i (2e-9); dropping W gives −0.2007 − 0.6127i, using A_χ instead of A_χ̄ gives 0.1975 − 0.1961i: the
phase W^{−1}(s/(1−s))γ_a(s) is right and no pairing is needed (logs/complex_chi5_order4.log; FE check 2e-21).
Other checks: γ_± closed forms and L(1−s,χ)/L(s,χ) = γ_a(s) to 1e-21 (χ₋₄, χ₅); Lemma Bχ's closed form vs direct oscillatory quadrature
(w = 0.4 + 5i, t = 0.3, 1.2, 3, 7.5) to 12 digits for both parities; the orthogonality identity (D₁a_χ, Y_s) = L(s,χ)/s at s = ½ + 3i:
χ₅ to 1e-9, χ₋₄ to 1e-9 (after adding the κ/x tail term, logs/burnol_m4_C_v2.log; 4e-4 without it); zeros: χ₋₄ 6.0209489, 10.2437703, …; χ₅ 6.6484533, 9.8314443, …; Σ_{|γ|<60} (χ₋₄) = 0.130739 +
N(T)-tail 0.0246 = 0.15537 ≈ C_{χ₋₄}; Σ_{|γ|<45} (χ₅) = 0.123823 + 0.0324 = 0.15621 ≈ C_{χ₅}.

## Review (reviewer subagents, the exact DERIVATION text pasted inline)
- Block 1, Burnol's liminf for χ (§1–5): GAP → fixed: "entire in w" → holomorphic for Re w < 1 (L(0,χ₋₄) = ½ gives a pole at w = 1); the
  off-line claim holds for nontrivial zeros only (L(−1,χ₋₄) = 0). Every other step OK (kernels, signs, Lemma Bχ, Y-existence, orthogonality,
  Gram, target 1/ρ, projection). ROUTINE: "the sine-kernel extension requires routine estimates, not a new idea".
- Block 2, upper-bound transfer (§6): GAP → fixed: scope (primitive, nonprincipal, conductor q: an induced χ₅ mod 10 has zeros on Re s = 0);
  one-sided zero count with χ̄ for negative ordinates; critical-line convexity (for L(s,χ)/s ∈ L² and ‖L R_tr‖); for complex χ, heights good at
  ±t and actual lower-half residues; R_tr meromorphic, PWS and maximum principle only for the entire G. (c1)–(c6) OK; no simplicity, central
  non-vanishing or L′(ρ) hypothesis needed. ROUTINE.
- Block 3, assembly, constant, equivalence, complex χ (§7–8): GAP → fixed: multiplicities in the Hadamard/Davenport formulas; explicit (iii),
  (iv) and combined R_high for complex χ. Equivalence (both directions), Re(1/ρ) = 1/(2|ρ|²), central zeros (4m² vs 4m), fixed-χ only: OK. ROUTINE.
- Second pass on the repaired §6–8: (a) scope, (b) counting, (c) convexity, (e) PWS/maximum principle CLOSED; (d) complex-χ growth bound in
  both half-planes and (f) m_ρ in L′/L OPEN → written out (DERIVATION §8 (iv), §6 (Cχ)); third pass: (d), (f) CLOSED, verdict CORRECT.
- Orchestrator's referees: §6 CORRECT ROUTINE; §7 CORRECT ROUTINE (C_{χ₋₄} = 0.155567979923586 recomputed); §8 GAP ((iv) dropped the opposite
  half's contour sum) → repaired: G = −Σ I^{w,+} − Σ I^{w,−} − R_tr, opposite-half terms bounded with separation |t| + t_m (DERIVATION §8); ROUTINE.
Real χ (the brief's theorem) does not depend on (d): for real χ, conj L(s̄,χ) = L(s,χ) and DbarNB's mirror tiling applies verbatim.

## Novelty
- Lower half (Theorem Bχ): KNOWN. De Roton, "Une approche hilbertienne de l'hypothèse de Riemann généralisée", Bull. SMF 134 (2006) 417–445,
  Thm 2 (read from numdam, open access): for F in the Selberg class with Ψ_F^{(1)} ∈ L²(0,∞), liminf D(λ)² log(1/λ) ≥ Σ_ρ m_ρ²/|ρ|² over the
  critical zeros, for D(λ) = dist_{L²(0,∞)}(1_{[0,1]}, span{Ψ_F(α/t) : λ ≤ α ≤ 1}); for F = L(s,χ), Ψ_F = −A_χ: our D_χ. Found via the GXL
  evidence API (claim 3: Dimitrov–Oliveira 1608.07887, TeX l.291: "De Roton extended the latter lower bound for L-functions in the Selberg
  class"). Our §§1–5 is an independent explicit proof for Dirichlet χ (dilated cosine/sine kernels).
- Upper half and the exact limit: not found. arXiv API (abs:Nyman AND abs:Dirichlet; abs:Nyman AND abs:"L-functions"; abs:Beurling AND Selberg
  class; Burnol): 1101.1199 (criterion, zero-free discs), 1608.07887 (L^p criteria, finite-zero extremal problem d_{n,m,p}² ~ Σ_{j≤m}/log n),
  1705.09918, math/0001013 (Burnol, adelic causality: criterion only); GXL claim 2 (exact asymptotic under GRH for Dirichlet L): 1 support,
  BCF 1211.5191, for ζ only; de Roton 2006 abstract and intro: lower bound only ("permet de majorer la vitesse de convergence"); arXiv listings
  2026-10/11 math.NT, math.CA (titles matching Nyman|Beurling|Báez-Duarte|Burnol|Dirichlet: 32 hits, none on the Nyman–Beurling distance);
  Zenodo "Nyman-Beurling" (60 records: ζ only); our files (TREE, NOGO, REPORTs, COROLLARIES: R520's readiness check; CeilAudit6's χ₋₄
  d_N numerics; NBLpK, which read de Roton 2007 for L^p): no continuous Dirichlet rate. Not read: de Roton JNT 129 (2009) (sequential version).
- Gate: the exact limit is new as a statement but needs no new idea beyond our ζ theorem (three ROUTINE verdicts): a corollary
  (COROLLARIES.md, the orchestrator's call), not a DISCOVERIES entry; the equivalent form (iii) likewise.

## Caveat
Ineffective (L₀ astronomically large, as for ζ; Lemma D invisible below L ≈ 10⁶–10⁸); not uniform in q. The upper half inherits DbarNB's black
boxes (Turán–Nazarov with the jet continuity remark, Titchmarsh 9.7/14.13 now for L(s,χ) "by the same proof", Claim F, explore59's (3i),
(3ii), §7). Numerics at L ≤ 16 are consistent but far from the limit (grid not converged; LB uses finitely many zeros). Complex χ (§8) rests on
the third-pass review of the two-sided tiling. The result does not move toward GRH itself.

RH relevance: N5 (the continuous Nyman–Beurling rate, now for every primitive Dirichlet L-function under its GRH; sibling of our ζ theorem).

## Files
Riemann/Deep/NBDir/: REPORT.md, DERIVATION.md (proof), nbdir_D.py (upper bound D_h², exact Gram), nbdir_burnol.py (identities, Lemma Bχ,
orthogonality, zeros, C_χ, finite-λ Burnol lower bound), nbdir_complex.py (§8 check), logs/ (D_m4_*, D_5_*, burnol_m4_AB,
burnol_m4_C_v2, burnol_m4_CDE(_v2), burnol_5_ABCDE, burnol_m4_DE_T100, complex_chi5_order4). Scratch: .scratch/nbdir/ (STATE.md; downloads deleted).
