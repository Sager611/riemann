AUDIT PASSED — every step of Theorem E1 (Steps 1–5 and all "other branch" claims) re-derived from Santibáñez-Leal's TeX, no error found; E1 holds exactly as written, given the same external inputs as the audited SL Thm 1.1 (Wang 2609.07918 Thm 2.2, Young's main-term evaluation, SL's parity product).

# E1Audit — referee of Theorem E1 (simple critical zeros in every window (T, T+T^θ], θ > 1/2)
e1audit (Orca), 12–14 Nov 2026. Object: Riemann/Deep/SantiAudit/E1.md. Source text: SL v0.03 main.tex (CAOS zip, Zenodo 23135354; line numbers in FULLPROOF.md), Wang v2 TeX.

## Short answer
E1-M (SL Thm 1.2 for every ν < (3θ−1)/4) and E1 (liminf S/N ≥ h(θ;κ1) ≥ h(1/2;κ1) = 0.0055353042386 on (1/2,1)) are correct. The new input is an absolute
x-space bound of SL's oscillatory dual term (14): |𝒯(n)| ≪ (T/H²)^j H (pq/(nT))^{1/4}(1+√(n/N0))^{−K}, total ≪ T^{1/2+ε}H^{−1/2}y². It is the classical
"trivial Voronoi" size (Steuding 2002 line; ν < 1/2 at θ = 1), now with complex shifts, general Q and SL's exact Gaussian transformation.

## What was checked (MECE)
| block | method | result |
|---|---|---|
| SL (10)–(14) as used by E1 | re-derived: Taylor profiles, Estermann FE (q = 1 reduction, standard form), u = 1−s−z/2 substitution, prefactor 4π/(gpq), branch phases | ✓ |
| Step 1: Mellin of J_d, Y_d, factor −π, contour, Parseval | DLMF 10.22.43–44 + reflection; strip (|Re d|/2, 1/4); per-n contour move, no poles, no horizontal terms | ✓ |
| Step 1–(12)–(14) end to end | voronoi_check: primal = residues + full dual, complex shifts, rel. error 1.7e−19 | ✓ (new test) |
| Step 2: Hankel form, symbol bounds for complex d | J, Y via H¹, H²; |(1±it/2X)^{d−1/2}| ≤ e^{π|Im d|/2}; X^k∂^k a_± = O_k(1) for all X > 0 | ✓ |
| Step 3: profile bounds uniform in H, j < J | ψ_j Schwartz uniformly (H^{2−m} ≤ 1); derivative costs 1/H; phase X = c√ξ non-stationary | ✓ |
| Step 4: integration by parts | each 𝒟 costs (T/X0)/H; T/(HX0) = (2√(2π))^{−1}√(N0/n) | ✓ |
| Step 5: n, g, p, q sums, σ_{−α,β}, both N0 regimes | K ≥ 3; pq powers cancel; Σ_g y²/g³ | ✓ |
| per-term scaling, j = 0, 1, 2, complex shifts, T to 1e6, θ = 0.505, 0.55, d = 0 | perterm_check(2): ratio ∈ [0.19, 0.66], tails ≤ 1.1e−3 at 8N0 | ✓ (E1 had j = 0, T ≤ 25600) |
| other branches: ε = +1, K-Bessel, residues, Taylor J, off-band, compact windows, shift derivatives | against SL l. 383–405, 444–447, 537–542, 707–758 | ✓ |
| Thm 1.3 (Lemma 4.2, Littlewood), (Emean), Thm 5.1 (Wang, 0 < λ < θ < 1) | SL l. 795–963, Wang TeX l. 229–240 | no θ ≥ 0.527, ν_A or H ≥ T^{1/2+c} requirement ✓ |
| constants κ1, c(1/2), k*(1/2), h(1/2;κ1), dh/dθ, h at 6 θ's, ν = 1249/10000 case | independent route (Chebyshev-form Q, Gauss–Legendre, Arb) | all digits agree ✓ |
| reviewer (Steps 4–5 pasted inline) | verdict-only | CORRECT |

## Falsification (smallest cases)
- θ = 0.501, ν = 3/25: dual exponent −0.0115 (−0.01075 at θ − η = 0.5005); needs fixed J ≥ 621, D ≥ 621, η < 0.001; all other exponents negative. No failure.
- Genuine failure points: ν ≥ (3θ−1)/4 and θ ≤ 1/2 — both excluded by E1's statement. d → 0 and N0 < 1 are harmless (strip nonempty; tail bound covers every n).
- Numerics never exceeded the claimed scaling; constants grow with j (as stated: liminf for fixed θ).

## Notes (no effect on the result)
1. The inverse-Mellin formula for 𝒢 on Re u < 1/4 converges only conditionally; the proof correctly rests on Parseval (G ξ^{σ−1} ∈ L¹, F̃ ∈ L¹(line)).
2. Step 3's "X ≥ 1 so Step 2 applies": unnecessary — the symbol bounds hold for all X > 0.
3. Purely asymptotic: actual shifts are −R/L, so |d| ≤ 1/4 needs L > 8R ≈ 71; j ≥ 1 residues at θ = 0.501 need T^{−0.001}L^5 → 0 (L ≳ 7·10⁴).
4. Not run (time): near-coalescent residue case of voronoi_check, T = 1e7 rows (cut by the 1200 s cap under background QoS); the K-branch is ≈ e^{−56} in the
   Voronoi test (SL's unchanged bound, re-derived by hand).
5. Unchanged external inputs (audited earlier in Riemann/Deep/SantiAudit, not re-audited here): Wang's Thm 2.2, the parity product [SL ref 4, Prop 2.5], Young's Lemma 7.

## What survives
Everything: E1 and E1-M exactly as written; least θ: every θ > 1/2; constants as in E1 §4 (κ1 = 0.0844779079697, h(1/2;κ1) = 0.0055353042386).

## Files
FULLPROOF.md (line-by-line notes); voronoi_check.py/.log; perterm_check.py/.log, perterm_check2.py/.log; constants_check.py/.log.

RH relevance: N8
